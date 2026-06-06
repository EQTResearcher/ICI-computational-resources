"""Critical slowing down diagnostics for R(t) time series."""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


def rolling_variance(series: np.ndarray, window: int = 20) -> np.ndarray:
    s = pd.Series(series)
    return s.rolling(window=window, min_periods=window).var().to_numpy()


def rolling_ar1(series: np.ndarray, window: int = 20) -> np.ndarray:
    n = len(series)
    ar1 = np.full(n, np.nan)
    for i in range(window - 1, n):
        sub = series[i - window + 1 : i + 1]
        x = sub[:-1] - sub[:-1].mean()
        y = sub[1:] - sub[1:].mean()
        denom = np.dot(x, x)
        if denom > 1e-10:
            ar1[i] = np.dot(x, y) / denom
    return ar1


def _hurst_rs(series: np.ndarray) -> Optional[float]:
    n = len(series)
    if n < 8:
        return None
    rs_vals = []
    for lag in range(4, n // 2 + 1):
        chunks = [series[i : i + lag] for i in range(0, n - lag + 1, lag)]
        chunk_vals = []
        for chunk in chunks:
            devs = np.cumsum(chunk - np.mean(chunk))
            R = np.max(devs) - np.min(devs)
            S = np.std(chunk, ddof=1)
            if S > 1e-10:
                chunk_vals.append(R / S)
        if chunk_vals:
            rs_vals.append((lag, float(np.mean(chunk_vals))))
    if len(rs_vals) < 3:
        return None
    lags_log = np.log([r[0] for r in rs_vals])
    rs_log = np.log([r[1] for r in rs_vals])
    slope, _, r, _, _ = stats.linregress(lags_log, rs_log)
    return float(slope) if abs(r) > 0.7 else None


def rolling_hurst(series: np.ndarray, window: int = 40) -> np.ndarray:
    n = len(series)
    vals = np.full(n, np.nan)
    for i in range(window - 1, n):
        h = _hurst_rs(series[i - window + 1 : i + 1])
        if h is not None:
            vals[i] = h
    return vals


def kendall_tau_trend(series: np.ndarray) -> Tuple[float, float]:
    clean = np.asarray(series)
    clean = clean[~np.isnan(clean)]
    if len(clean) < 5:
        return 0.0, 1.0
    tau, p = stats.kendalltau(np.arange(len(clean)), clean)
    return float(tau), float(p)


def detect_hurst_transition(hurst_series: np.ndarray, threshold: float = 0.5, min_delta: float = 0.2) -> bool:
    valid = hurst_series[~np.isnan(hurst_series)]
    if len(valid) < 10:
        return False
    mid = len(valid) // 2
    mean_first = np.mean(valid[:mid])
    mean_last = np.mean(valid[mid:])
    return bool(mean_first < threshold and mean_last > threshold and (mean_last - mean_first) > min_delta)


def compute_csd_signals(
    Rt_series: np.ndarray,
    window_var: int = 20,
    window_ar1: int = 20,
    window_hurst: int = 40,
) -> Dict:
    Rt_series = np.asarray(Rt_series, dtype=float)
    var_s = rolling_variance(Rt_series, window=window_var)
    ar1_s = rolling_ar1(Rt_series, window=window_ar1)
    hurst_s = rolling_hurst(Rt_series, window=window_hurst)

    tau_var, p_var = kendall_tau_trend(var_s)
    tau_ar1, p_ar1 = kendall_tau_trend(ar1_s)
    hurst_trans = detect_hurst_transition(hurst_s)

    score = 0
    Rt_recent = float(np.nanmean(Rt_series[-5:])) if len(Rt_series) >= 5 else float(Rt_series[-1])
    if Rt_recent < -8.5:
        score += 2
    if Rt_recent < -9.5:
        score += 2
    if tau_var > 0.4 and p_var < 0.05:
        score += 1
    if tau_ar1 > 0.4 and p_ar1 < 0.05:
        score += 1
    if hurst_trans:
        score += 2

    if score <= 1:
        alert = "green"
    elif score <= 3:
        alert = "yellow"
    elif score <= 5:
        alert = "orange"
    else:
        alert = "red"

    return {
        "variance": var_s,
        "ar1": ar1_s,
        "hurst": hurst_s,
        "tau_var": tau_var,
        "p_var": p_var,
        "tau_ar1": tau_ar1,
        "p_ar1": p_ar1,
        "hurst_trans": hurst_trans,
        "alert_score": score,
        "alert_level": alert,
    }
