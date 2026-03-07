from __future__ import annotations
import numpy as np
from .beam import GaussBeam


def find_waists(trajectory: list[GaussBeam]) -> list[GaussBeam]:
    """軌跡中のビームウエスト（R が負→正に転じる点）を返す。

    ウエストでは R が -inf から +inf へ遷移する（収束→発散）。
    レンズ等による正→負の遷移はウエストではないため除外する。
    """
    Rs = np.array([b.R_mm for b in trajectory])
    signs = np.sign(Rs)
    flips = np.where(np.diff(signs) > 0)[0] + 1  # -1 -> +1 のみ
    return [trajectory[int(i)] for i in flips]


def rayleigh_range(beam: GaussBeam) -> float:
    """Rayleigh長 z_R (mm)。"""
    return np.pi * beam.n * beam.w_mm**2 / (beam.wl_um * 1e-3)


def confocal_parameter(beam: GaussBeam) -> float:
    """共焦点パラメータ 2*z_R (mm)。"""
    return 2.0 * rayleigh_range(beam)
