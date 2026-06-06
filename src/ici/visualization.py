"""Visualization utilities for ICI and R(t)."""
from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

ICI_COLORS = {
    "green": "#2E7D32",
    "yellow": "#F9A825",
    "orange": "#E65100",
    "red": "#B71C1C",
    "blue": "#1565C0",
    "gray": "#546E7A",
}


def plot_rt_timeseries(
    time_points: list,
    Rt_values: np.ndarray,
    Rt_lower: Optional[np.ndarray] = None,
    Rt_upper: Optional[np.ndarray] = None,
    system_name: str = "",
    events: Optional[dict] = None,
    figsize: tuple = (12, 5),
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize)
    ax.axhspan(-9.0, 5.0, alpha=0.08, color=ICI_COLORS["green"], label="Resilient state (R > -9)")
    ax.axhspan(-10.0, -9.0, alpha=0.12, color=ICI_COLORS["yellow"], label="Metastable band (-10 <= R < -9)")
    ax.axhspan(-15.0, -10.0, alpha=0.12, color=ICI_COLORS["red"], label="Collapse zone (R < -10)")
    ax.axhline(-9.0, color=ICI_COLORS["yellow"], lw=1.5, ls="--", alpha=0.7)
    ax.axhline(-10.0, color=ICI_COLORS["red"], lw=2.0, ls="--", alpha=0.9, label="Critical threshold R = -10")

    if Rt_lower is not None and Rt_upper is not None:
        ax.fill_between(time_points, Rt_lower, Rt_upper, alpha=0.25, color=ICI_COLORS["blue"], label="68% CI")

    for i in range(len(Rt_values) - 1):
        mid = (Rt_values[i] + Rt_values[i + 1]) / 2
        color = ICI_COLORS["green"] if mid > -9 else ICI_COLORS["yellow"] if mid > -10 else ICI_COLORS["red"]
        ax.plot(time_points[i : i + 2], Rt_values[i : i + 2], color=color, lw=2.5, solid_capstyle="round")

    for t, r in zip(time_points, Rt_values):
        color = ICI_COLORS["green"] if r > -9 else ICI_COLORS["yellow"] if r > -10 else ICI_COLORS["red"]
        ax.scatter(t, r, color=color, s=60, zorder=5, edgecolors="white", linewidths=1.0)

    if events:
        for t, label in events.items():
            ax.axvline(t, color="black", lw=1.0, ls=":", alpha=0.6)
            ax.annotate(label, xy=(t, ax.get_ylim()[1] * 0.95), xytext=(5, 0), textcoords="offset points", fontsize=9, rotation=45)

    ax.set_xlabel("Time")
    ax.set_ylabel("R(t)")
    ax.set_title(f"R(t) trajectory: {system_name}")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9, ncol=2)
    ax.set_ylim(min(Rt_values) - 1.5, 2.0)
    ax.grid(True, alpha=0.3, ls="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return fig


def plot_ici_spectrum(systems: list[str], ici_values: list[float], ici_errors: Optional[list[float]] = None, figsize: tuple = (10, 7)) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize)
    order = np.argsort(ici_values)
    sys_sorted = [systems[i] for i in order]
    ici_sorted = [ici_values[i] for i in order]
    err_sorted = [ici_errors[i] for i in order] if ici_errors else None

    def ici_color(val: float) -> str:
        if val < 50:
            return ICI_COLORS["gray"]
        if val < 5000:
            return ICI_COLORS["blue"]
        if val < 8000:
            return "#7B1FA2"
        return ICI_COLORS["red"]

    y_pos = np.arange(len(sys_sorted))
    ax.barh(y_pos, ici_sorted, color=[ici_color(v) for v in ici_sorted], alpha=0.8, height=0.7)
    if err_sorted:
        ax.errorbar(ici_sorted, y_pos, xerr=err_sorted, fmt="none", color="black", capsize=4, lw=1.5)
    ax.axvspan(0, 50, alpha=0.04, color="gray", label="Primary layer (0-50)")
    ax.axvspan(50, 5000, alpha=0.04, color=ICI_COLORS["blue"], label="Intermediate layer (50-5000)")
    ax.axvspan(5000, 8000, alpha=0.04, color="purple", label="Advanced layer (5000-8000)")
    ax.axvspan(8000, 12000, alpha=0.04, color="red", label="Top layer (8000+)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sys_sorted, fontsize=10)
    ax.set_xscale("log")
    ax.set_xlabel("ICI value, log scale")
    ax.set_title("System complexity spectrum: ICI comparison")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, axis="x", alpha=0.3, which="both")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return fig
