from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from .beam import GaussBeam


class Element(ABC):
    @property
    def label(self) -> str:
        return self.__class__.__name__

    @property
    def param_str(self) -> str:
        return ""

    @property
    def dz(self) -> float:
        return 0.0

    @property
    @abstractmethod
    def matrix(self) -> np.ndarray:
        ...

    def _q_transform(self, q: complex) -> complex:
        A, B, C, D = self.matrix.ravel()
        return (A * q + B) / (C * q + D)

    def apply(self, beam: GaussBeam) -> GaussBeam:
        q_new = self._q_transform(beam.q)
        return GaussBeam.from_q(beam.wl_um, beam.n, q_new, z_mm=beam.z_mm + self.dz)


class Propagation(Element):
    def __init__(self, d_mm: float):
        self.d = float(d_mm)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, self.d], [0.0, 1.0]])

    @property
    def label(self) -> str:
        return "Propagation"

    @property
    def param_str(self) -> str:
        return f"d = {self.d:>8.3f} mm"

    @property
    def dz(self) -> float:
        return self.d


class ThinLens(Element):
    def __init__(self, f_mm: float):
        self.f = float(f_mm)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, 0.0], [-1.0 / self.f, 1.0]])

    @property
    def label(self) -> str:
        return "ThinLens"

    @property
    def param_str(self) -> str:
        return f"f = {self.f:>8.3f} mm"


class Interface(Element):
    def __init__(self, n1: float, n2: float):
        self.n1 = float(n1)
        self.n2 = float(n2)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, 0.0], [0.0, self.n1 / self.n2]])

    @property
    def label(self) -> str:
        return "Interface"

    @property
    def param_str(self) -> str:
        return f"n {self.n1:.4f} -> {self.n2:.4f}"

    def apply(self, beam: GaussBeam) -> GaussBeam:
        q_new = self._q_transform(beam.q)
        return GaussBeam.from_q(beam.wl_um, self.n2, q_new, z_mm=beam.z_mm)


class InterfaceCurved(Element):
    def __init__(self, n1: float, n2: float, r_mm: float):
        self.n1 = float(n1)
        self.n2 = float(n2)
        self.r = float(r_mm)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, 0.0],
                         [(self.n2 - self.n1) / (self.n2 * self.r), self.n1 / self.n2]])

    @property
    def label(self) -> str:
        return "InterfaceCurved"

    @property
    def param_str(self) -> str:
        return f"r={self.r:.3f} mm, n {self.n1:.4f}->{self.n2:.4f}"

    def apply(self, beam: GaussBeam) -> GaussBeam:
        q_new = self._q_transform(beam.q)
        return GaussBeam.from_q(beam.wl_um, self.n2, q_new, z_mm=beam.z_mm)


class CurvedMirrorTan(Element):
    def __init__(self, r_mm: float, theta_deg: float):
        self.r = float(r_mm)
        self.theta_deg = float(theta_deg)

    @property
    def _r_eff(self) -> float:
        return self.r * np.cos(self.theta_deg * np.pi / 180)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, 0.0], [-2.0 / self._r_eff, 1.0]])

    @property
    def label(self) -> str:
        return "CurvedMirrorTan"

    @property
    def param_str(self) -> str:
        return f"r={self.r:.3f} mm, th={self.theta_deg:.1f} deg"


class CurvedMirrorSag(Element):
    def __init__(self, r_mm: float, theta_deg: float):
        self.r = float(r_mm)
        self.theta_deg = float(theta_deg)

    @property
    def _r_eff(self) -> float:
        return self.r / np.cos(self.theta_deg * np.pi / 180)

    @property
    def matrix(self) -> np.ndarray:
        return np.array([[1.0, 0.0], [-2.0 / self._r_eff, 1.0]])

    @property
    def label(self) -> str:
        return "CurvedMirrorSag"

    @property
    def param_str(self) -> str:
        return f"r={self.r:.3f} mm, th={self.theta_deg:.1f} deg"
