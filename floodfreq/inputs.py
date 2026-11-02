"""Reading annual maxima from pasted text or CSV files."""
import re
import numpy as np
import pandas as pd

# 1,200 / 12,500 / 1,234,567.5 are thousands separators; 520,610,480 is a list.
_THOUSANDS = re.compile(r"^\d{1,2}(,\d{3})+(\.\d+)?$")
YEAR_NAMES = {"year", "years", "yr", "water_year", "wateryear", "wy", "date", "time"}


def parse_numbers(text):
    """Parse numbers separated by spaces, newlines, semicolons or commas.

    Returns (array, notes). Tokens like '1,200' are read as 1200 (and reported in notes).
    """
    vals, notes = [], []
    for tok in re.split(r"[\s;]+", text.strip()):
        if not tok:
            continue
        if _THOUSANDS.match(tok):
            notes.append(f"read '{tok}' as {tok.replace(',', '')}")
            parts = [tok.replace(",", "")]
        else:
            parts = [p for p in tok.split(",") if p]
        for p in parts:
            try:
                vals.append(float(p))
            except ValueError:
                raise ValueError(f"Could not read '{p}' as a number.") from None
    return np.array(vals), notes


def _is_year_like(name, s):
    if str(name).strip().lower() in YEAR_NAMES:
        return True
    if pd.api.types.is_numeric_dtype(s) and s.notna().all() and len(s) >= 5:
        v = s.to_numpy(float)
        return bool(np.all(v == np.round(v)) and v.min() >= 1800 and v.max() <= 2200
                    and (np.all(np.diff(v) > 0) or np.all(np.diff(v) < 0)))
    return False


def read_csv_values(path, column=None):
    """Read annual maxima from a CSV. Year-like columns are detected by name or content, used to
    put the records in chronological order, and never analysed as values.

    Returns (values, years or None, notes, value_column_name).
    """
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("The CSV file is empty.")
    year_col = next((c for c in df.columns if _is_year_like(c, df[c])), None)
    notes = []
    if column is not None:
        if column not in df.columns:
            raise ValueError(f"Column '{column}' not found. Columns: {list(df.columns)}")
        value_col = column
        if year_col == column:
            year_col = None
    else:
        cands = [c for c in df.columns if c != year_col and pd.api.types.is_numeric_dtype(df[c])]
        if not cands:
            raise ValueError("No numeric value column found in the CSV.")
        value_col = cands[-1]
        if len(cands) > 1:
            notes.append(f"several numeric columns {cands}; using '{value_col}' (choose with --column)")
    if year_col is not None:
        key = df[year_col]
        if not pd.api.types.is_numeric_dtype(key):
            parsed = pd.to_datetime(key, errors="coerce")
            key = parsed if parsed.notna().all() else None
        if key is not None:
            df = df.assign(_k=key).sort_values("_k")
        else:
            year_col = None
    sub = df.dropna(subset=[value_col])
    values = pd.to_numeric(sub[value_col], errors="raise").to_numpy(float)
    years = sub[year_col].tolist() if year_col is not None else None
    if year_col is not None:
        notes.append(f"year column '{year_col}' used for chronological order")
    return values, years, notes, value_col
