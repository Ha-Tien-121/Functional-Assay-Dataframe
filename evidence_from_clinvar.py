"""
Hybrid Evidence Strength Calibration using ClinVar as truth set.

This module calculates LR+ and LR- (odds of pathogenicity) for functional assays
using two methods:
1. THRESHOLD-BASED: Uses MaveDB-curated score thresholds when available
2. CLASSIFICATION-BASED: Falls back to author-reported classifications

Output:
- Per-variant evidence columns added to master dataframe
- Per-assay summary with calibration method, LR values, and performance metrics
"""

import argparse
import pathlib
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

PSEUDOCOUNT = 0.5

# Tavtigian / Brnich likelihood ratio thresholds
LR_PATHOGENIC_THRESHOLDS = {
    "very_strong": 350,
    "strong": 18.7,
    "moderate": 4.3,
    "supporting": 2.1,
}

LR_BENIGN_THRESHOLDS = {
    "strong": 0.053,
    "moderate": 0.23,
    "supporting": 0.48,
}

# Token lists used to normalize assay class values
ABNORMAL_TOKENS = [
    "abnormal",
    "loss",
    "lof",
    "deleterious",
    "damaging",
    "nonfunctional",
    "non-functional",
    "non functional",
    "deficient",
    "inactive",
    "pathogenic",
    "LOF",
    "likely pathogenic",
    "high_risk"
    "reduced", "compromised", "severe", "strong functional effect",
]
NORMAL_TOKENS = [
    "normal",
    "wt",
    "wild",
    "functional",
    "neutral",
    "no effect",
    "like wt",
    "wildtype",
    "wild type",
    "FUNC",
    "likely not pathogenic",
    "not pathogenic",
    "likely neutral"
    "benign", "no functional effect",
]
INDETERMINATE_TOKENS = [
    "indeterminate",
    "intermediate",
    "uncertain",
    "ambiguous",
    "partial",
    "mixed",
    "INT",
    "not cleard",
    "low_risk"
    "not specified",
    "unclear",
]

# Mapping of MaveDB assay names to score columns in master dataframe
# Format: "MaveDB Assay Name" -> "score_column_name"
MAVEDB_TO_SCORE_COLUMN = {
    "Findlay et al 2018 BRCA1": "BRCA1_Findlay_auth_func_score",
    "Adamovich et al 2022 BRCA1 HDR": "BRCA1_Adamovich_HDR_auth_func_score",
    "Adamovich et al 2022 BRCA1 Cisplatin": "BRCA1_Adamovich_Cisplatin_auth_func_score",
    "Fernandes et al 2019 BRCA1": "BRCA1_Fernandes_eta",
    "Bouwman et al 2020 BRCA1 HDR": "BRCA1_Bouwman_2020_nGFP+b",
    "Bouwman et al 2020 BRCA1 Cisplatin Sensitivty": "BRCA1_Bouwman_2020_Cisplatin_nIC50b",
    "Bouwman et al 2020 BRCA1 Olaparib Sensitivity": "BRCA1_Bouwman_2020_Olaparib_nIC50b",
    # Add more mappings as needed for other genes
}

# Mapping of MaveDB assay names to classification columns (fallback)
MAVEDB_TO_CLASS_COLUMN = {
    "Findlay et al 2018 BRCA1": "BRCA1_Findlay_reported_functional_class",
    "Adamovich et al 2022 BRCA1 HDR": "BRCA1_Adamovich_HDR_auth_reported_functional_class",
    "Adamovich et al 2022 BRCA1 Cisplatin": "BRCA1_Adamovich_Cisplatin_auth_reported_functional_class",
    "Fernandes et al 2019 BRCA1": "BRCA1_Fernandes_fClass_Category",
    "Bouwman et al 2020 BRCA1 HDR": "BRCA1_Bouwman_2020_DR-GFP_prediction",
    "Bouwman et al 2020 BRCA1 Cisplatin Sensitivty": "BRCA1_Bouwman_2020_Cisplatin_prediction",
    "Bouwman et al 2020 BRCA1 Olaparib Sensitivity": "BRCA1_Bouwman_2020_Olaparib_prediction",
    # Add more mappings as needed
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ThresholdInterval:
    """Represents a single threshold interval from MaveDB."""
    name: str
    lower: float  # -inf if unbounded
    upper: float  # inf if unbounded
    lower_inclusive: bool
    upper_inclusive: bool
    mavedb_class: str  # 'Normal', 'Abnormal', 'Not specified', etc.
    
    def contains(self, value: float) -> bool:
        """Check if value falls within this interval."""
        if pd.isna(value):
            return False
        lower_ok = (value > self.lower) if not self.lower_inclusive else (value >= self.lower)
        upper_ok = (value < self.upper) if not self.upper_inclusive else (value <= self.upper)
        return lower_ok and upper_ok


@dataclass
class AssayThresholds:
    """Collection of threshold intervals for an assay."""
    assay_name: str
    intervals: List[ThresholdInterval] = field(default_factory=list)
    score_column: Optional[str] = None
    class_column: Optional[str] = None
    
    def get_normal_threshold(self) -> Optional[float]:
        """Get the boundary for 'Normal' classification."""
        for iv in self.intervals:
            if iv.mavedb_class.lower() == 'normal':
                # Return the more restrictive bound (the one closer to abnormal)
                if iv.lower > -np.inf:
                    return iv.lower
                if iv.upper < np.inf:
                    return iv.upper
        return None
    
    def get_abnormal_threshold(self) -> Optional[float]:
        """Get the boundary for 'Abnormal' classification."""
        for iv in self.intervals:
            if iv.mavedb_class.lower() == 'abnormal':
                if iv.upper < np.inf:
                    return iv.upper
                if iv.lower > -np.inf:
                    return iv.lower
        return None
    
    def classify_score(self, score: float) -> Optional[str]:
        """Classify a score value using the thresholds."""
        if pd.isna(score):
            return None
        for iv in self.intervals:
            if iv.contains(score):
                mclass = iv.mavedb_class.lower().strip()
                if mclass in ['normal']:
                    return 'normal'
                elif mclass in ['abnormal']:
                    return 'abnormal'
                elif mclass in ['not specified', 'intermediate']:
                    return 'indeterminate'
        return 'indeterminate'  # Default if no interval matches


@dataclass 
class AssayCalibrationResult:
    """Results of calibrating a single assay."""
    assay_name: str
    calibration_method: str  # 'threshold-based' or 'classification-based'
    score_column: Optional[str]
    class_column: Optional[str]
    n_pathogenic: int
    n_benign: int
    n_abn_P: int
    n_abn_B: int
    n_norm_P: int
    n_norm_B: int
    n_indeterminate: int
    lr_plus: float
    lr_minus: float
    evidence_strength_path: str
    evidence_strength_benign: str
    sensitivity: Optional[float]
    specificity: Optional[float]
    precision: Optional[float]
    npv: Optional[float]
    indeterminate_rate: Optional[float]
    thresholds_used: Optional[str] = None  # Description of thresholds


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_master_df(path: str) -> pd.DataFrame:
    """Load a master dataframe (CSV/CSV.GZ/PARQUET)."""
    file_path = pathlib.Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Master dataframe not found: {file_path}")
    
    suffix = file_path.suffix.lower()
    name = file_path.name.lower()
    if suffix in {".parquet"} or name.endswith(".parquet.gz"):
        return pd.read_parquet(file_path)
    return pd.read_csv(file_path, low_memory=False)


def detect_clinvar_column(df: pd.DataFrame, provided: Optional[str] = None) -> str:
    """Resolve which ClinVar significance column to use."""
    if provided:
        if provided not in df.columns:
            raise KeyError(f"Provided ClinVar column '{provided}' not found.")
        return provided
    
    preferred = ["clinvar_sig_2025", "clinvar_sig_2018"]
    for col in preferred:
        if col in df.columns:
            return col
    
    candidates = [
        c for c in df.columns
        if "clinvar" in c.lower() and ("sig" in c.lower() or "significance" in c.lower())
    ]
    if candidates:
        return candidates[0]
    raise KeyError("No ClinVar significance column found.")


def classify_clinvar_sig(value: object) -> Optional[str]:
    """Map a ClinVar significance string to {'pathogenic','benign',None}."""
    if pd.isna(value):
        return None
    
    val = str(value).lower()
    if not val.strip():
        return None
    
    # Exclusions
    if "conflict" in val:
        return None
    if "uncertain" in val or "vus" in val:
        return None
    if "risk" in val or "drug" in val or "association" in val:
        return None
    
    has_path = "pathogenic" in val
    has_benign = "benign" in val
    
    if has_path and has_benign:
        return None  # ambiguous
    if has_path:
        return "pathogenic"
    if has_benign:
        return "benign"
    return None


def normalize_assay_class(value: object) -> Optional[str]:
    """Normalize an assay functional class value to {'abnormal','normal','indeterminate',None}."""
    if pd.isna(value):
        return None
    
    val = str(value).strip().lower()
    if not val or val in {"na", "nan", "none", ""}:
        return None
    
    # Check abnormal first (more specific tokens)
    if any(tok in val for tok in ABNORMAL_TOKENS):
        return "abnormal"
    if any(tok in val for tok in NORMAL_TOKENS):
        return "normal"
    if any(tok in val for tok in INDETERMINATE_TOKENS):
        return "indeterminate"
    
    # Special cases
    if val in {"lof", "func", "int"}:
        if val == "lof":
            return "abnormal"
        elif val == "func":
            return "normal"
        elif val == "int":
            return "indeterminate"
    
    if val in {"+", "-"}:
        # Context-dependent, treat as indeterminate
        return "indeterminate"
    
    return None


def evidence_from_lr_path(lr_plus: float) -> str:
    """Evidence strength derived from LR+ (pathogenic evidence)."""
    if pd.isna(lr_plus):
        return "Indeterminate"
    if lr_plus > LR_PATHOGENIC_THRESHOLDS["very_strong"]:
        return "Pathogenic Very Strong"
    if lr_plus > LR_PATHOGENIC_THRESHOLDS["strong"]:
        return "Pathogenic Strong"
    if lr_plus > LR_PATHOGENIC_THRESHOLDS["moderate"]:
        return "Pathogenic Moderate"
    if lr_plus > LR_PATHOGENIC_THRESHOLDS["supporting"]:
        return "Pathogenic Supporting"
    return "Indeterminate"


def evidence_from_lr_benign(lr_minus: float) -> str:
    """Evidence strength derived from LR- (benign evidence)."""
    if pd.isna(lr_minus):
        return "Indeterminate"
    if lr_minus < LR_BENIGN_THRESHOLDS["strong"]:
        return "Benign Strong"
    if lr_minus < LR_BENIGN_THRESHOLDS["moderate"]:
        return "Benign Moderate"
    if lr_minus < LR_BENIGN_THRESHOLDS["supporting"]:
        return "Benign Supporting"
    return "Indeterminate"


def _smoothed_rate(count: int, total: int) -> float:
    """Apply Haldane-Anscombe-like smoothing: (count + 0.5) / (total + 1.0)"""
    return (count + PSEUDOCOUNT) / (total + 1.0)


def _safe_div(num: float, den: float) -> Optional[float]:
    """Safe division returning None if denominator is zero."""
    return num / den if den else None


# =============================================================================
# THRESHOLD PARSING
# =============================================================================

def parse_interval_range(range_str: str) -> Tuple[float, float, bool, bool]:
    """
    Parse a MaveDB interval range string like "[-0.748, Inf)" or "(0.05, 0.95]".
    
    Returns: (lower, upper, lower_inclusive, upper_inclusive)
    """
    if pd.isna(range_str) or not range_str.strip():
        return (-np.inf, np.inf, False, False)
    
    s = str(range_str).strip()
    
    # Determine inclusivity from brackets
    lower_inclusive = s.startswith('[')
    upper_inclusive = s.endswith(']')
    
    # Remove brackets
    s = s.strip('[]()').strip()
    
    # Split by comma
    parts = [p.strip() for p in s.split(',')]
    if len(parts) != 2:
        return (-np.inf, np.inf, False, False)
    
    # Parse lower bound
    lower_str = parts[0].lower().replace('inf', 'inf').replace('-inf', '-inf')
    if 'inf' in lower_str and '-' in lower_str:
        lower = -np.inf
    elif 'inf' in lower_str:
        lower = np.inf
    else:
        try:
            lower = float(lower_str)
        except ValueError:
            lower = -np.inf
    
    # Parse upper bound
    upper_str = parts[1].lower().replace('inf', 'inf')
    if 'inf' in upper_str and '-' not in upper_str:
        upper = np.inf
    elif '-inf' in upper_str:
        upper = -np.inf
    else:
        try:
            upper = float(upper_str)
        except ValueError:
            upper = np.inf
    
    return (lower, upper, lower_inclusive, upper_inclusive)


def extract_mavedb_assays(df: pd.DataFrame) -> List[str]:
    """Extract unique MaveDB assay names from interval columns."""
    assays = set()
    pattern = re.compile(r'^Interval \d+ (name|range|MaveDB class) (.+)$')
    
    for col in df.columns:
        match = pattern.match(col)
        if match:
            assay_name = match.group(2)
            # Filter out the "class X" duplicates
            if not assay_name.startswith('class '):
                assays.add(assay_name)
    
    return sorted(assays)


def parse_assay_thresholds(df: pd.DataFrame, assay_name: str) -> AssayThresholds:
    """Parse all threshold intervals for a given assay from the dataframe."""
    thresholds = AssayThresholds(assay_name=assay_name)
    
    # Try to find score and class columns
    thresholds.score_column = MAVEDB_TO_SCORE_COLUMN.get(assay_name)
    thresholds.class_column = MAVEDB_TO_CLASS_COLUMN.get(assay_name)
    
    # Parse intervals 1-6
    for i in range(1, 7):
        name_col = f"Interval {i} name {assay_name}"
        range_col = f"Interval {i} range {assay_name}"
        class_col = f"Interval {i} MaveDB class {assay_name}"
        
        if name_col not in df.columns:
            continue
        
        # Get first non-null values (they should be constant across rows)
        name_val = df[name_col].dropna().iloc[0] if df[name_col].notna().any() else None
        range_val = df[range_col].dropna().iloc[0] if range_col in df.columns and df[range_col].notna().any() else None
        class_val = df[class_col].dropna().iloc[0] if class_col in df.columns and df[class_col].notna().any() else None
        
        if name_val is None:
            continue
        
        lower, upper, lower_inc, upper_inc = parse_interval_range(range_val)
        
        interval = ThresholdInterval(
            name=str(name_val).strip(),
            lower=lower,
            upper=upper,
            lower_inclusive=lower_inc,
            upper_inclusive=upper_inc,
            mavedb_class=str(class_val).strip() if class_val else "Not specified"
        )
        thresholds.intervals.append(interval)
    
    return thresholds


# =============================================================================
# CALIBRATION LOGIC
# =============================================================================

def calibrate_assay_threshold_based(
    df: pd.DataFrame,
    thresholds: AssayThresholds,
    clinvar_truth_col: str = "clinvar_truth"
) -> Optional[AssayCalibrationResult]:
    """
    Calibrate an assay using score thresholds.
    Returns None if threshold-based calibration is not possible.
    """
    score_col = thresholds.score_column
    
    # Check if we can do threshold-based calibration
    if score_col is None or score_col not in df.columns:
        return None
    if not thresholds.intervals:
        return None
    
    # Check if we have any Normal/Abnormal intervals defined
    has_normal = any(iv.mavedb_class.lower() == 'normal' for iv in thresholds.intervals)
    has_abnormal = any(iv.mavedb_class.lower() == 'abnormal' for iv in thresholds.intervals)
    if not (has_normal or has_abnormal):
        return None
    
    # Create standardized classification column
    std_col = f"{thresholds.assay_name}__std_threshold"
    df[std_col] = df[score_col].apply(thresholds.classify_score)
    
    # Filter to rows with both truth and classification
    mask = df[clinvar_truth_col].notna() & df[std_col].notna()
    subset = df.loc[mask]
    
    if len(subset) < 5:  # Minimum sample requirement
        return None
    
    # Compute counts
    nP = (subset[clinvar_truth_col] == "pathogenic").sum()
    nB = (subset[clinvar_truth_col] == "benign").sum()
    
    n_abn_P = ((subset[std_col] == "abnormal") & (subset[clinvar_truth_col] == "pathogenic")).sum()
    n_abn_B = ((subset[std_col] == "abnormal") & (subset[clinvar_truth_col] == "benign")).sum()
    n_norm_P = ((subset[std_col] == "normal") & (subset[clinvar_truth_col] == "pathogenic")).sum()
    n_norm_B = ((subset[std_col] == "normal") & (subset[clinvar_truth_col] == "benign")).sum()
    n_indet = (subset[std_col] == "indeterminate").sum()
    
    # Compute LR+ and LR- with smoothing
    p_abn_path = _smoothed_rate(n_abn_P, nP)
    p_abn_ben = _smoothed_rate(n_abn_B, nB)
    p_norm_path = _smoothed_rate(n_norm_P, nP)
    p_norm_ben = _smoothed_rate(n_norm_B, nB)
    
    lr_plus = p_abn_path / p_abn_ben if p_abn_ben != 0 else np.nan
    lr_minus = p_norm_path / p_norm_ben if p_norm_ben != 0 else np.nan
    
    # Classification metrics
    tp, fp = n_abn_P, n_abn_B
    tn, fn = n_norm_B, n_norm_P
    
    # Build threshold description
    thresh_desc_parts = []
    for iv in thresholds.intervals:
        if iv.mavedb_class.lower() in ['normal', 'abnormal']:
            thresh_desc_parts.append(f"{iv.mavedb_class}: [{iv.lower}, {iv.upper}]")
    thresh_desc = "; ".join(thresh_desc_parts)
    
    return AssayCalibrationResult(
        assay_name=thresholds.assay_name,
        calibration_method="threshold-based",
        score_column=score_col,
        class_column=None,
        n_pathogenic=int(nP),
        n_benign=int(nB),
        n_abn_P=int(n_abn_P),
        n_abn_B=int(n_abn_B),
        n_norm_P=int(n_norm_P),
        n_norm_B=int(n_norm_B),
        n_indeterminate=int(n_indet),
        lr_plus=lr_plus,
        lr_minus=lr_minus,
        evidence_strength_path=evidence_from_lr_path(lr_plus),
        evidence_strength_benign=evidence_from_lr_benign(lr_minus),
        sensitivity=_safe_div(tp, tp + fn),
        specificity=_safe_div(tn, tn + fp),
        precision=_safe_div(tp, tp + fp),
        npv=_safe_div(tn, tn + fn),
        indeterminate_rate=_safe_div(n_indet, len(subset)),
        thresholds_used=thresh_desc
    )


def calibrate_assay_classification_based(
    df: pd.DataFrame,
    thresholds: AssayThresholds,
    clinvar_truth_col: str = "clinvar_truth"
) -> Optional[AssayCalibrationResult]:
    """
    Calibrate an assay using author-reported classifications.
    Returns None if classification-based calibration is not possible.
    """
    class_col = thresholds.class_column
    
    if class_col is None or class_col not in df.columns:
        return None
    
    # Create standardized classification column
    std_col = f"{thresholds.assay_name}__std_class"
    df[std_col] = df[class_col].apply(normalize_assay_class)
    
    # Filter to rows with both truth and classification
    mask = df[clinvar_truth_col].notna() & df[std_col].notna()
    subset = df.loc[mask]
    
    if len(subset) < 5:
        return None
    
    # Compute counts
    nP = (subset[clinvar_truth_col] == "pathogenic").sum()
    nB = (subset[clinvar_truth_col] == "benign").sum()
    
    n_abn_P = ((subset[std_col] == "abnormal") & (subset[clinvar_truth_col] == "pathogenic")).sum()
    n_abn_B = ((subset[std_col] == "abnormal") & (subset[clinvar_truth_col] == "benign")).sum()
    n_norm_P = ((subset[std_col] == "normal") & (subset[clinvar_truth_col] == "pathogenic")).sum()
    n_norm_B = ((subset[std_col] == "normal") & (subset[clinvar_truth_col] == "benign")).sum()
    n_indet = (subset[std_col] == "indeterminate").sum()
    
    # Compute LR+ and LR-
    p_abn_path = _smoothed_rate(n_abn_P, nP)
    p_abn_ben = _smoothed_rate(n_abn_B, nB)
    p_norm_path = _smoothed_rate(n_norm_P, nP)
    p_norm_ben = _smoothed_rate(n_norm_B, nB)
    
    lr_plus = p_abn_path / p_abn_ben if p_abn_ben != 0 else np.nan
    lr_minus = p_norm_path / p_norm_ben if p_norm_ben != 0 else np.nan
    
    # Classification metrics
    tp, fp = n_abn_P, n_abn_B
    tn, fn = n_norm_B, n_norm_P
    
    return AssayCalibrationResult(
        assay_name=thresholds.assay_name,
        calibration_method="classification-based",
        score_column=None,
        class_column=class_col,
        n_pathogenic=int(nP),
        n_benign=int(nB),
        n_abn_P=int(n_abn_P),
        n_abn_B=int(n_abn_B),
        n_norm_P=int(n_norm_P),
        n_norm_B=int(n_norm_B),
        n_indeterminate=int(n_indet),
        lr_plus=lr_plus,
        lr_minus=lr_minus,
        evidence_strength_path=evidence_from_lr_path(lr_plus),
        evidence_strength_benign=evidence_from_lr_benign(lr_minus),
        sensitivity=_safe_div(tp, tp + fn),
        specificity=_safe_div(tn, tn + fp),
        precision=_safe_div(tp, tp + fp),
        npv=_safe_div(tn, tn + fn),
        indeterminate_rate=_safe_div(n_indet, len(subset)),
        thresholds_used=None
    )


def calibrate_assay_hybrid(
    df: pd.DataFrame,
    thresholds: AssayThresholds,
    clinvar_truth_col: str = "clinvar_truth"
) -> Optional[AssayCalibrationResult]:
    """
    Hybrid calibration: try threshold-based first, fall back to classification-based.
    """
    # Try threshold-based first
    result = calibrate_assay_threshold_based(df, thresholds, clinvar_truth_col)
    if result is not None:
        return result
    
    # Fall back to classification-based
    result = calibrate_assay_classification_based(df, thresholds, clinvar_truth_col)
    return result


# =============================================================================
# AUTO-DETECTION OF ADDITIONAL ASSAY COLUMNS
# =============================================================================

def auto_detect_class_columns(df: pd.DataFrame) -> Dict[str, str]:
    """
    Auto-detect assay classification columns not covered by MaveDB mappings.
    Returns dict of {assay_name: column_name}.
    """
    detected = {}
    mavedb_class_cols = set(MAVEDB_TO_CLASS_COLUMN.values())
    
    for col in df.columns:
        name = col.lower()
        # Skip ClinVar and already-mapped columns
        if "clinvar" in name:
            continue
        if col in mavedb_class_cols:
            continue
        
        # Check if it looks like a classification column
        if any(x in name for x in ["func_class", "_class", "functional", "classification", "prediction", "category"]):
            # Extract assay name (heuristic: everything before the classification suffix)
            assay_name = col
            for suffix in ["_func_class", "_class", "_Functional_classification", "_prediction", "_fClass_Category"]:
                if col.endswith(suffix):
                    assay_name = col[:-len(suffix)]
                    break
            detected[assay_name] = col
    
    return detected


# =============================================================================
# MAIN CALIBRATION PIPELINE
# =============================================================================

def build_calibration_summary(
    df: pd.DataFrame,
    clinvar_col: Optional[str] = None,
    include_auto_detected: bool = True
) -> Tuple[pd.DataFrame, Dict[str, AssayCalibrationResult]]:
    """
    Build comprehensive calibration summary for all detected assays.
    
    Returns:
        - DataFrame with per-assay summary
        - Dict mapping assay_name -> AssayCalibrationResult
    """
    # Detect ClinVar column and add truth column
    clinvar_col = detect_clinvar_column(df, clinvar_col)
    df = df.copy()
    df["clinvar_truth"] = df[clinvar_col].apply(classify_clinvar_sig)
    
    print(f"Using ClinVar column: {clinvar_col}")
    print(f"  Pathogenic variants: {(df['clinvar_truth'] == 'pathogenic').sum()}")
    print(f"  Benign variants: {(df['clinvar_truth'] == 'benign').sum()}")
    
    results: Dict[str, AssayCalibrationResult] = {}
    
    # 1. Calibrate MaveDB assays (with thresholds)
    mavedb_assays = extract_mavedb_assays(df)
    print(f"\nFound {len(mavedb_assays)} MaveDB assays with threshold intervals")
    
    for assay_name in mavedb_assays:
        thresholds = parse_assay_thresholds(df, assay_name)
        result = calibrate_assay_hybrid(df, thresholds, "clinvar_truth")
        
        if result:
            results[assay_name] = result
            print(f"  [OK] {assay_name}: {result.calibration_method} (nP={result.n_pathogenic}, nB={result.n_benign})")
        else:
            print(f"  [--] {assay_name}: could not calibrate (insufficient data or missing columns)")
    
    # 2. Auto-detect additional classification columns
    if include_auto_detected:
        auto_detected = auto_detect_class_columns(df)
        print(f"\nAuto-detected {len(auto_detected)} additional classification columns")
        
        for assay_name, class_col in auto_detected.items():
            if assay_name in results:
                continue  # Already calibrated
            
            # Create minimal thresholds object for classification-based
            thresholds = AssayThresholds(
                assay_name=assay_name,
                class_column=class_col
            )
            result = calibrate_assay_classification_based(df, thresholds, "clinvar_truth")
            
            if result:
                results[assay_name] = result
                print(f"  [OK] {assay_name}: {result.calibration_method} (nP={result.n_pathogenic}, nB={result.n_benign})")
    
    # Convert results to DataFrame
    records = []
    for name, res in results.items():
        records.append({
            "assay_name": res.assay_name,
            "calibration_method": res.calibration_method,
            "score_column": res.score_column,
            "class_column": res.class_column,
            "n_pathogenic": res.n_pathogenic,
            "n_benign": res.n_benign,
            "n_abn_P": res.n_abn_P,
            "n_abn_B": res.n_abn_B,
            "n_norm_P": res.n_norm_P,
            "n_norm_B": res.n_norm_B,
            "n_indeterminate": res.n_indeterminate,
            "LR_plus": res.lr_plus,
            "LR_minus": res.lr_minus,
            "evidence_strength_path": res.evidence_strength_path,
            "evidence_strength_benign": res.evidence_strength_benign,
            "sensitivity": res.sensitivity,
            "specificity": res.specificity,
            "precision": res.precision,
            "npv": res.npv,
            "indeterminate_rate": res.indeterminate_rate,
            "thresholds_used": res.thresholds_used,
        })
    
    summary_df = pd.DataFrame.from_records(records)
    return summary_df, results


def add_evidence_columns(
    df: pd.DataFrame,
    calibration_results: Dict[str, AssayCalibrationResult],
    clinvar_col: Optional[str] = None
) -> pd.DataFrame:
    """
    Add per-variant evidence columns to the master dataframe.
    
    For each calibrated assay, adds:
    - {assay}__LR_value
    - {assay}__evidence_strength
    """
    df = df.copy()
    clinvar_col = detect_clinvar_column(df, clinvar_col)
    df["clinvar_truth"] = df[clinvar_col].apply(classify_clinvar_sig)
    
    evidence_cols_added = []
    
    for assay_name, result in calibration_results.items():
        lr_plus = result.lr_plus
        lr_minus = result.lr_minus
        
        # Determine which column to use for classification
        if result.calibration_method == "threshold-based" and result.score_column:
            # Re-parse thresholds and classify scores
            thresholds = parse_assay_thresholds(df, assay_name)
            std_col = f"{assay_name}__std"
            df[std_col] = df[result.score_column].apply(thresholds.classify_score)
        elif result.class_column:
            std_col = f"{assay_name}__std"
            df[std_col] = df[result.class_column].apply(normalize_assay_class)
        else:
            continue
        
        # Create LR value column
        lr_col = f"{assay_name}__LR_value"
        ev_col = f"{assay_name}__evidence_strength"
        
        def compute_lr(val):
            if val == "abnormal":
                return lr_plus
            elif val == "normal":
                return lr_minus
            return 1.0  # Neutral
        
        def compute_evidence(val):
            if val == "abnormal":
                return evidence_from_lr_path(lr_plus)
            elif val == "normal":
                return evidence_from_lr_benign(lr_minus)
            return "Indeterminate"
        
        df[lr_col] = df[std_col].apply(compute_lr)
        df[ev_col] = df[std_col].apply(compute_evidence)
        evidence_cols_added.append(ev_col)
    
    # Create summary column of applied evidence
    def gather_sources(row):
        sources = []
        for col in evidence_cols_added:
            if col in row and row[col] not in ["Indeterminate", None, np.nan]:
                assay = col.replace("__evidence_strength", "")
                sources.append(assay)
        return ";".join(sources) if sources else ""
    
    df["evidence_sources_applied"] = df.apply(gather_sources, axis=1)
    
    return df


# =============================================================================
# CLI ENTRYPOINT
# =============================================================================

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hybrid evidence strength calibration using ClinVar as truth set."
    )
    parser.add_argument(
        "--master-path", type=str, required=True,
        help="Path to master dataframe (csv/csv.gz/parquet)"
    )
    parser.add_argument(
        "--clinvar-col", type=str, default=None,
        help="ClinVar significance column (auto-detect if omitted)"
    )
    parser.add_argument(
        "--output-master", type=str, default=None,
        help="Output path for updated master dataframe"
    )
    parser.add_argument(
        "--summary-output", type=str, default="assay_calibration_summary.csv",
        help="Output CSV for per-assay calibration summary"
    )
    parser.add_argument(
        "--no-auto-detect", action="store_true",
        help="Disable auto-detection of additional classification columns"
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    master_path = pathlib.Path(args.master_path)
    
    print(f"Loading master dataframe from: {master_path}")
    df = load_master_df(master_path)
    print(f"Loaded {len(df)} variants, {len(df.columns)} columns")
    
    # Build calibration summary
    summary_df, results = build_calibration_summary(
        df,
        clinvar_col=args.clinvar_col,
        include_auto_detected=not args.no_auto_detect
    )
    
    print(f"\n=== CALIBRATION SUMMARY ===")
    print(f"Total assays calibrated: {len(results)}")
    threshold_based = sum(1 for r in results.values() if r.calibration_method == "threshold-based")
    class_based = sum(1 for r in results.values() if r.calibration_method == "classification-based")
    print(f"  Threshold-based: {threshold_based}")
    print(f"  Classification-based: {class_based}")
    
    # Add evidence columns
    print("\nAdding per-variant evidence columns...")
    updated_df = add_evidence_columns(df, results, args.clinvar_col)
    
    # Save outputs
    output_master = args.output_master
    if output_master is None:
        stem = master_path.stem
        if stem.endswith(".csv"):
            stem = stem[:-4]
        output_master = master_path.parent / f"{stem}_with_evidence.csv.gz"
    else:
        output_master = pathlib.Path(output_master)
    
    compression = "gzip" if str(output_master).endswith(".gz") else None
    updated_df.to_csv(output_master, index=False, compression=compression)
    print(f"Saved updated master dataframe: {output_master}")
    
    summary_df.to_csv(args.summary_output, index=False)
    print(f"Saved calibration summary: {args.summary_output}")
    
    # Print summary table
    print("\n=== PER-ASSAY RESULTS ===")
    display_cols = ["assay_name", "calibration_method", "n_pathogenic", "n_benign", 
                    "LR_plus", "LR_minus", "evidence_strength_path", "evidence_strength_benign"]
    print(summary_df[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
