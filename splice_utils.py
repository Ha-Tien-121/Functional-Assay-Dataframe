import numpy as np
import pandas as pd
from typing import Iterable, Set

# Columns to check; accept both dotted and upper-case variants.
SPLICE_COL_CANDIDATES = [
    "spliceai.dp_ag",
    "spliceai.dp_al",
    "spliceai.dp_dg",
    "spliceai.dp_dl",
    
]


def compute_spliceai_max(df: pd.DataFrame) -> pd.Series:
    """
    Row-wise max across available SpliceAI columns; returns NaN if none present.
    """
    present = [c for c in SPLICE_COL_CANDIDATES if c in df.columns]
    if not present:
        return pd.Series([np.nan] * len(df), index=df.index, name="spliceai_max")
    return df[present].astype(float).max(axis=1, skipna=True)


def get_spliceai_exclusion_mask(df: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
    """
    True where spliceai_max > threshold. Computes spliceai_max if missing.
    """
    if "spliceai_max" not in df.columns:
        splice_max = compute_spliceai_max(df)
    else:
        splice_max = df["spliceai_max"]
    return splice_max > threshold


def null_out_assay_for_splice(dataset_df: pd.DataFrame, exclude_keys: Iterable, join_key: str = "join_key") -> int:
    """
    Set all assay columns (except join_key) to NaN for rows whose join_key is in exclude_keys.
    Returns number of rows nulled.
    """
    if dataset_df is None or dataset_df.empty or join_key not in dataset_df.columns:
        return 0
    mask = dataset_df[join_key].isin(exclude_keys)
    count = int(mask.sum())
    if count:
        cols = [c for c in dataset_df.columns if c != join_key]
        dataset_df.loc[mask, cols] = np.nan
    return count

