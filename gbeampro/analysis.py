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


def beam_at(trajectory: list[GaussBeam], z_mm: float) -> GaussBeam:
    """軌跡上の位置 z_mm におけるビームを返す。

    z_mm 以下で最後の軌跡点から自由伝搬させて求める（軌跡終端より先は自由空間として外挿）。
    """
    from .elements import Propagation

    if z_mm < trajectory[0].z_mm:
        raise ValueError(f"z_mm ({z_mm}) は軌跡の開始位置 ({trajectory[0].z_mm}) より前です")
    b = [t for t in trajectory if t.z_mm <= z_mm][-1]
    d = z_mm - b.z_mm
    return Propagation(d).apply(b) if d > 0 else b


def aperture_loss(beam: GaussBeam, r_mm: float, beam_y: GaussBeam | None = None) -> float:
    """半径 r_mm の円形アパーチャ（光軸中心）で遮られるパワーの割合 (0–1)。

    beam_y を与えると x 方向半径 beam.w_mm，y 方向半径 beam_y.w_mm の楕円ビームとして扱う。
    """
    wx = beam.w_mm
    wy = wx if beam_y is None else beam_y.w_mm
    if np.isclose(wx, wy):
        return float(np.exp(-2.0 * r_mm**2 / wx**2))
    # 極座標で動径方向を解析積分し，角度方向は周期関数の台形則で数値積分する
    phi = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    a = np.cos(phi)**2 / wx**2 + np.sin(phi)**2 / wy**2
    T = np.mean(-np.expm1(-2.0 * r_mm**2 * a) / a) / (wx * wy)
    return float(1.0 - T)
