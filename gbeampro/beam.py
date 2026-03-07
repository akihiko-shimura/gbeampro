from __future__ import annotations
from dataclasses import dataclass
import numpy as np

_pi = np.pi


@dataclass(frozen=True)
class GaussBeam:
    wl_um: float
    n: float = 1.0
    z_mm: float = 0.0
    R_mm: float = np.inf
    w_mm: float = 1.0

    @property
    def q(self) -> complex:
        return 1.0 / (1j * self.wl_um * 1e-3 / (_pi * self.n * self.w_mm**2) + 1.0 / self.R_mm)

    @property
    def theta(self) -> float:
        return np.arctan(self.wl_um * 1e-3 / (_pi * self.n * self.w_mm))

    @classmethod
    def from_waist(cls, wl_um: float, w0_mm: float, z_mm: float = 0.0, n: float = 1.0) -> GaussBeam:
        return cls(wl_um=wl_um, n=n, z_mm=z_mm, R_mm=np.inf, w_mm=w0_mm)

    @classmethod
    def from_q(cls, wl_um: float, n: float, q: complex, z_mm: float = 0.0) -> GaussBeam:
        qinv = 1.0 / q
        qinv_real = float(np.real(qinv))
        qinv_imag = float(np.imag(qinv))
        R_mm = float(1.0 / qinv_real) if qinv_real != 0.0 else np.inf
        w_mm = float(np.sqrt(wl_um * 1e-3 / (_pi * n * abs(qinv_imag))))
        return cls(wl_um=wl_um, n=n, z_mm=z_mm, R_mm=R_mm, w_mm=w_mm)

    def __repr__(self) -> str:
        R_str = "inf" if not np.isfinite(self.R_mm) else f"{self.R_mm:.5e}"
        return (
            f"GaussBeam(wl_um={self.wl_um}, n={self.n}, z_mm={self.z_mm:.5f}, "
            f"R_mm={R_str}, w_mm={self.w_mm:.5f})\n"
            f"  q = {self.q.real:.5e}{self.q.imag:+.5e}j  "
            f"theta = {self.theta * 1e3:.4f} mrad"
        )
