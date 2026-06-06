"""Cross-domain validation tools: mutual information and Procrustes alignment."""
from __future__ import annotations

from typing import Optional
import warnings

import numpy as np
from scipy.spatial import procrustes

try:
    from npeet.entropy_estimators import mi as ksg_mi
    HAS_NPEET = True
except Exception:  # pragma: no cover
    HAS_NPEET = False
    ksg_mi = None


def _approx_mi(x: np.ndarray, y: np.ndarray, bins: int = 10) -> float:
    if x.ndim > 1:
        x = x[:, 0]
    if y.ndim > 1:
        y = y[:, 0]
    hist_xy, _, _ = np.histogram2d(x, y, bins=bins)
    if hist_xy.sum() == 0:
        return 0.0
    pxy = hist_xy / hist_xy.sum()
    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)
    mi = 0.0
    for i in range(pxy.shape[0]):
        for j in range(pxy.shape[1]):
            if pxy[i, j] > 1e-12 and px[i] > 1e-12 and py[j] > 1e-12:
                mi += pxy[i, j] * np.log2(pxy[i, j] / (px[i] * py[j]))
    return float(mi)


def mutual_information(x: np.ndarray, y: np.ndarray, k: int = 5) -> float:
    x = np.asarray(x)
    y = np.asarray(y)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    if HAS_NPEET:
        return float(ksg_mi(x.tolist(), y.tolist(), k=k))
    warnings.warn("npeet is not installed; using a lower-precision histogram MI approximation.")
    return _approx_mi(x, y)


def functional_equivalence_test(
    source_params: np.ndarray,
    target_params: np.ndarray,
    param_names: Optional[list[str]] = None,
    threshold_bits: float = 0.5,
    n_permutations: int = 200,
) -> dict:
    if param_names is None:
        param_names = ["D", "C", "S", "F", "W", "M"]
    source_params = np.asarray(source_params, dtype=float)
    target_params = np.asarray(target_params, dtype=float)
    if source_params.shape[1] != 6 or target_params.shape[1] != 6:
        raise ValueError("source_params and target_params must each have six columns.")

    rng = np.random.default_rng(42)
    results = {}
    passed_count = 0
    for i, name in enumerate(param_names):
        src = np.log10(np.clip(source_params[:, i], 1e-300, None))
        tgt = np.log10(np.clip(target_params[:, i], 1e-300, None))
        mi_val = mutual_information(src, tgt)
        null_mis = [mutual_information(src, rng.permutation(tgt)) for _ in range(n_permutations)]
        p_val = float(np.mean(np.array(null_mis) >= mi_val))
        passed = mi_val >= threshold_bits and p_val < 0.05
        results[name] = {
            "mutual_information": mi_val,
            "p_value": p_val,
            "passed": bool(passed),
            "threshold": threshold_bits,
        }
        passed_count += int(passed)
    results["overall"] = {
        "n_passed": passed_count,
        "n_total": 6,
        "verdict": "PASS" if passed_count >= 4 else "FAIL",
        "note": f"{passed_count}/6 parameters passed functional equivalence test.",
    }
    return results


def topological_alignment_test(
    source_params: np.ndarray,
    target_params: np.ndarray,
    threshold_r2: float = 0.75,
    n_permutations: int = 500,
) -> dict:
    source_params = np.asarray(source_params, dtype=float)
    target_params = np.asarray(target_params, dtype=float)

    def standardize(x: np.ndarray) -> np.ndarray:
        mu = x.mean(axis=0)
        sd = x.std(axis=0)
        sd[sd < 1e-10] = 1.0
        return (x - mu) / sd

    log_src = standardize(np.log10(np.clip(source_params, 1e-300, None)))
    log_tgt = standardize(np.log10(np.clip(target_params, 1e-300, None)))
    _, _, disparity = procrustes(log_src, log_tgt)
    r2 = 1.0 - disparity

    rng = np.random.default_rng(42)
    null_r2 = []
    for _ in range(n_permutations):
        _, _, d_null = procrustes(log_src, log_tgt[rng.permutation(len(log_tgt))])
        null_r2.append(1.0 - d_null)
    p_val = float(np.mean(np.array(null_r2) >= r2))
    passed = r2 >= threshold_r2 and p_val < 0.05
    return {
        "procrustes_r2": float(r2),
        "disparity": float(disparity),
        "p_value": p_val,
        "threshold": threshold_r2,
        "passed": bool(passed),
        "verdict": "PASS" if passed else "FAIL",
    }


def full_crossdomain_validation(source_params: np.ndarray, target_params: np.ndarray) -> dict:
    feq = functional_equivalence_test(source_params, target_params)
    topo = topological_alignment_test(source_params, target_params)
    both = feq["overall"]["verdict"] == "PASS" and topo["verdict"] == "PASS"
    return {
        "functional_equivalence": feq,
        "topological_alignment": topo,
        "combined_verdict": "VALID" if both else "INVALID",
    }
