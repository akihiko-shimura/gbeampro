from __future__ import annotations
import numpy as np
from .beam import GaussBeam
from .elements import Element, Propagation


class OpticalSystem:
    def __init__(self):
        self.elements: list[Element] = []

    def add(self, element: Element) -> OpticalSystem:
        self.elements.append(element)
        return self

    def trace(self, beam: GaussBeam, dz: float = 0.01) -> list[GaussBeam]:
        """全軌跡を返す。Propagation素子はdz刻みで中間状態も含む。"""
        traj = [beam]
        current = beam
        for el in self.elements:
            if isinstance(el, Propagation) and el.d > 0:
                steps = max(1, round(el.d / dz))
                step_el = Propagation(el.d / steps)
                for _ in range(steps):
                    current = step_el.apply(current)
                    traj.append(current)
            else:
                current = el.apply(current)
                traj.append(current)
        return traj

    def __str__(self) -> str:
        W = 57
        lines = [
            "OpticalSystem",
            "=" * W,
            f" {'#':>3}  {'Type':<18}  {'Parameters':<22}  {'z (mm)':>8}",
            "-" * W,
        ]
        z = 0.0
        lines.append(f" {'0':>3}  {'--- input ---':<18}  {'':22}  {z:>8.3f}")
        for i, el in enumerate(self.elements, 1):
            z += el.dz
            lines.append(f" {i:>3}  {el.label:<18}  {el.param_str:<22}  {z:>8.3f}")
        lines.append("=" * W)
        lines.append(f"Total length: {z:.3f} mm  |  {len(self.elements)} elements")
        return "\n".join(lines)

    def summary(self, beam: GaussBeam, dz: float = 0.01) -> str:
        """各素子の境界でのビーム状態表とウエスト情報を返す。"""
        from .analysis import find_waists

        checkpoints: list[tuple[str, GaussBeam]] = [("--- input ---", beam)]
        full_traj: list[GaussBeam] = [beam]
        current = beam

        for el in self.elements:
            if isinstance(el, Propagation) and el.d > 0:
                steps = max(1, round(el.d / dz))
                step_el = Propagation(el.d / steps)
                for _ in range(steps):
                    current = step_el.apply(current)
                    full_traj.append(current)
                checkpoints.append((el.label, current))
            else:
                current = el.apply(current)
                full_traj.append(current)
                checkpoints.append((el.label, current))

        W = 72
        lines = [
            f"OpticalSystem trace  [wl={beam.wl_um:.4g} um]",
            "=" * W,
            f" {'#':>3}  {'Type':<18}  {'z (mm)':>9}  {'w (um)':>9}  {'R (mm)':>10}  {'th (urad)':>10}",
            "-" * W,
        ]
        for i, (name, b) in enumerate(checkpoints):
            R_str = "inf" if not np.isfinite(b.R_mm) else f"{b.R_mm:.3e}"
            lines.append(
                f" {i:>3}  {name:<18}  {b.z_mm:>9.3f}"
                f"  {b.w_mm * 1e3:>9.2f}  {R_str:>10}  {b.theta * 1e6:>10.2f}"
            )
        lines.append("=" * W)

        waists = find_waists(full_traj)
        if waists:
            parts = [f"z={w.z_mm:.3f} mm (2w0={w.w_mm * 2e3:.1f} um)" for w in waists]
            lines.append("Beam waists:  " + "  |  ".join(parts))

        return "\n".join(lines)
