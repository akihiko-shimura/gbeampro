from __future__ import annotations

import numpy as np

from .beam import GaussBeam
from .elements import (
    CurvedMirrorSag,
    CurvedMirrorTan,
    Element,
    Interface,
    InterfaceCurved,
    ThinLens,
)


def plot_caustic(trajectory: list[GaussBeam], ax, **kwargs):
    """caustic曲線 (z, w) をプロット。"""
    zs = np.array([b.z_mm for b in trajectory])
    ws = np.array([b.w_mm * 1e3 for b in trajectory])
    ax.plot(zs, ws, **kwargs)
    ax.set_xlabel("$z$ (mm)")
    ax.set_ylabel("$w$ (µm)")


def plot_system(system, trajectory: list[GaussBeam], ax, beam_kw: dict | None = None, label: str = ""):
    """caustic + 各光学素子のシンボルを重ねてプロット。"""
    zs = np.array([b.z_mm for b in trajectory])
    ws = np.array([b.w_mm * 1e3 for b in trajectory])
    kw = beam_kw or {}

    kw_no_label = {k: v for k, v in kw.items() if k != "label"}
    (line,) = ax.plot(zs, ws, label=label or "beam", **kw)
    color = line.get_color()
    ax.fill_between(zs, ws, -ws, alpha=0.2, color=color, **{k: v for k, v in kw_no_label.items() if k != "color"})
    ax.plot(zs, -ws, color=color, **{k: v for k, v in kw_no_label.items() if k != "color"})

    h = float(np.max(ws)) * 1.2
    z_cur = trajectory[0].z_mm
    for el in system.elements:
        z_cur += el.dz
        _draw_element(ax, el, z_cur, h, color)

    ax.axhline(0, color="k", lw=0.5, ls="--")
    ax.set_xlabel("$z$ (mm)")
    ax.set_ylabel("$w$ (µm)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=10, frameon=False)


def _draw_element(ax, el: Element, z: float, h: float, color: str):
    if isinstance(el, ThinLens):
        ax.plot([z, z], [-h, h], color=color, lw=1.5, zorder=5)
        ax.plot(z, h,  marker="^", color=color, ms=7, zorder=5)
        ax.plot(z, -h, marker="v", color=color, ms=7, zorder=5)
        ax.text(z, h * 1.05, f"f={el.f:.0f}", ha="center", va="bottom", fontsize=7, color=color)
    elif isinstance(el, (Interface, InterfaceCurved)):
        ax.axvline(z, color=color, lw=1.5, ls="--", alpha=0.8)
    elif isinstance(el, (CurvedMirrorTan, CurvedMirrorSag)):
        ax.axvline(z, color=color, lw=2.0, alpha=0.7)
