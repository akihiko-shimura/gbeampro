from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .beam import GaussBeam
from .elements import Propagation, ThinLens
from .system import OpticalSystem

# ── Operand ──────────────────────────────────────────────────────────────────

_OPERAND_TYPES = {'wx', 'wy', 'cvx', 'cvy', 'thx', 'thy'}


@dataclass
class Operand:
    """Single term in the merit function.

    type:
        wx  / wy  : beam radius (mm) at z_mm for x / y axis
        cvx / cvy : wavefront curvature 1/R (mm⁻¹) at z_mm; target=0 → waist
        thx / thy : half-divergence angle (mrad) at z_mm
    target:
        desired value in the unit implied by type.
    weight:
        relative importance; merit contribution ∝ weight × residual².
    """
    type: str
    z_mm: float
    target: float
    weight: float = 1.0

    def __post_init__(self):
        if self.type not in _OPERAND_TYPES:
            raise ValueError(f"operand type '{self.type}' must be one of {_OPERAND_TYPES}")


# ── Result ────────────────────────────────────────────────────────────────────

@dataclass
class OptResult:
    specs: list[dict]
    merit: float
    operand_values: list[float]
    success: bool


# ── Internal helpers ──────────────────────────────────────────────────────────

LensSpec = dict  # {'type': str, 'z_mm': float, 'f_mm': float}


def _build_xy_systems(
    beam: GaussBeam,
    lens_specs: list[LensSpec],
    z_target_mm: float,
) -> tuple[OpticalSystem, OpticalSystem]:
    specs = sorted(lens_specs, key=lambda s: s["z_mm"])
    sys_x, sys_y = OpticalSystem(), OpticalSystem()
    z_prev = beam.z_mm

    for s in specs:
        d = float(s["z_mm"]) - z_prev
        if d > 0:
            sys_x.add(Propagation(d))
            sys_y.add(Propagation(d))
        t, f = s["type"], float(s["f_mm"])
        if t == "spherical":
            sys_x.add(ThinLens(f))
            sys_y.add(ThinLens(f))
        elif t == "cyl_x":
            sys_x.add(ThinLens(f))
        elif t == "cyl_y":
            sys_y.add(ThinLens(f))
        z_prev = float(s["z_mm"])

    d_final = z_target_mm - z_prev
    if d_final > 0:
        sys_x.add(Propagation(d_final))
        sys_y.add(Propagation(d_final))

    return sys_x, sys_y


def _eval_operand(op: Operand, bx: GaussBeam, by: GaussBeam) -> float:
    if op.type == 'wx':
        return bx.w_mm
    elif op.type == 'wy':
        return by.w_mm
    elif op.type == 'cvx':
        return 1.0 / bx.R_mm if np.isfinite(bx.R_mm) else 0.0
    elif op.type == 'cvy':
        return 1.0 / by.R_mm if np.isfinite(by.R_mm) else 0.0
    elif op.type == 'thx':
        return bx.theta * 1e3  # rad → mrad
    elif op.type == 'thy':
        return by.theta * 1e3


def _residual_scale(op: Operand) -> float:
    """Normalization factor so residuals are dimensionless and O(1) near the target."""
    if op.type in ('wx', 'wy'):
        return abs(op.target) if op.target != 0 else 1.0
    elif op.type in ('cvx', 'cvy'):
        # 1/R → 0 at waist; scale by 1/z_mm so residual ≈ z/R (dimensionless)
        return 1.0 / max(op.z_mm, 1.0)
    elif op.type in ('thx', 'thy'):
        return abs(op.target) if op.target != 0 else 1.0
    return 1.0


# ── Main optimizer ────────────────────────────────────────────────────────────

def optimize_astigmatic(
    beam: GaussBeam,
    lens_types: list[str],
    operands: list[Operand],
    z_bounds: tuple[float, float] | None = None,
    f_bounds: tuple[float, float] | None = None,
    f_abs_bounds: tuple[float, float] | None = None,
    min_lens_sep_mm: float = 0.0,
    algorithm: str = 'de',
    seed: int = 42,
    maxiter: int = 1000,
    popsize: int = 15,
    tol: float = 1e-12,
) -> OptResult:
    """Optimize an astigmatic Gaussian beam system against a merit function.

    Parameters
    ----------
    beam : GaussBeam
        Input beam (circular; shared for x and y axes).
    lens_types : list of str
        Lens types: 'spherical', 'cyl_x', or 'cyl_y'.
    operands : list of Operand
        Merit function terms. Use ``waist_operands()`` for the common case.
    z_bounds : (float, float), optional
        Search range for lens z positions (mm). Default: (beam.z_mm, max operand z).
    f_bounds : (float, float), optional
        Search range for focal lengths (mm). Default: (1.0, 10000.0).
    f_abs_bounds : (float, float), optional
        Constraint on |f| (mm). If given, lenses with |f| outside this range are
        penalized. f_bounds is automatically set to (-f_abs_bounds[1], f_abs_bounds[1])
        unless explicitly provided. Example: f_abs_bounds=(30, 1000) allows
        f ∈ [-1000, -30] ∪ [30, 1000].
    min_lens_sep_mm : float
        Minimum separation between adjacent lenses (mm).
    algorithm : str
        'de'     — differential evolution (global, robust)
        'lm'     — Levenberg-Marquardt / TRF (local, fast, gradient-based)
        'hammer' — DE global search → TRF polish on top candidates (like Zemax Hammer)
    seed, maxiter, popsize, tol
        Algorithm tuning parameters.

    Returns
    -------
    OptResult
    """
    try:
        from scipy.optimize import differential_evolution, least_squares
    except ImportError:
        raise ImportError("scipy が必要です: pip install gbeampro[optimize]")

    n = len(lens_types)
    z_lo, z_hi = z_bounds if z_bounds is not None else (
        beam.z_mm, max(op.z_mm for op in operands)
    )
    if f_abs_bounds is not None:
        f_abs_lo, f_abs_hi = f_abs_bounds
        f_lo, f_hi = f_bounds if f_bounds is not None else (-f_abs_hi, f_abs_hi)
    else:
        f_abs_lo = 0.0
        f_lo, f_hi = f_bounds if f_bounds is not None else (1.0, 10000.0)
    bounds_list = [(z_lo, z_hi), (f_lo, f_hi)] * n
    bounds_lo = np.array([b[0] for b in bounds_list])
    bounds_hi = np.array([b[1] for b in bounds_list])

    def _build_specs(params: np.ndarray) -> list[LensSpec]:
        zs = np.sort(params[0::2])
        fs = params[1::2]
        return [{"type": t, "z_mm": float(z), "f_mm": float(f)}
                for t, z, f in zip(lens_types, zs, fs)]

    def _sep_ok(params: np.ndarray) -> bool:
        if min_lens_sep_mm <= 0 or n < 2:
            return True
        zs = np.sort(params[0::2])
        return bool(np.all(np.diff(zs) >= min_lens_sep_mm))

    def _fabs_ok(params: np.ndarray) -> bool:
        if f_abs_lo <= 0:
            return True
        fs = params[1::2]
        return bool(np.all(np.abs(fs) >= f_abs_lo))

    def _residuals(params: np.ndarray) -> np.ndarray:
        if not _sep_ok(params) or not _fabs_ok(params):
            return np.full(len(operands), 1e3)
        specs = _build_specs(params)
        # Trace once per unique evaluation z
        cache: dict[float, tuple[GaussBeam, GaussBeam]] = {}
        res = []
        for op in operands:
            if op.z_mm not in cache:
                sx, sy = _build_xy_systems(beam, specs, op.z_mm)
                dz = max(1.0, op.z_mm / 100)
                cache[op.z_mm] = sx.trace(beam, dz=dz)[-1], sy.trace(beam, dz=dz)[-1]
            bx, by = cache[op.z_mm]
            val = _eval_operand(op, bx, by)
            scale = _residual_scale(op)
            res.append(np.sqrt(op.weight) * (val - op.target) / scale)
        return np.array(res)

    def _scalar(params: np.ndarray) -> float:
        r = _residuals(params)
        return float(np.dot(r, r))

    # ── Algorithm dispatch ────────────────────────────────────────────────

    if algorithm == 'de':
        de = differential_evolution(
            _scalar, bounds_list,
            tol=tol, seed=seed, maxiter=maxiter, popsize=popsize,
            polish=True, workers=1,
        )
        best_params = de.x
        success = bool(de.success)

    elif algorithm == 'lm':
        rng = np.random.default_rng(seed)
        x0 = rng.uniform(bounds_lo, bounds_hi)
        r = least_squares(
            _residuals, x0,
            bounds=(bounds_lo, bounds_hi),
            method='trf', max_nfev=maxiter * 100, ftol=tol,
        )
        best_params = r.x
        success = bool(r.success)

    elif algorithm == 'hammer':
        # Phase 1: DE global search (coarse)
        de = differential_evolution(
            _scalar, bounds_list,
            tol=tol * 1e3, seed=seed, maxiter=maxiter // 2, popsize=popsize,
            polish=False, workers=1, init='sobol',
        )
        # Phase 2: TRF polish on top candidates from DE population
        candidates = [de.x]
        top_idx = np.argsort(de.population_energies)[:min(5, len(de.population_energies))]
        for idx in top_idx:
            candidates.append(de.population[idx])

        best_fun, best_params = np.inf, de.x
        for x0 in candidates:
            x0 = np.clip(x0, bounds_lo, bounds_hi)
            r = least_squares(
                _residuals, x0,
                bounds=(bounds_lo, bounds_hi),
                method='trf', max_nfev=maxiter * 50, ftol=tol,
            )
            if r.cost < best_fun:
                best_fun, best_params = r.cost, r.x
        success = True

    else:
        raise ValueError(f"algorithm='{algorithm}' unknown. Use 'de', 'lm', or 'hammer'.")

    # ── Assemble result ───────────────────────────────────────────────────

    if not _sep_ok(best_params):
        return OptResult(specs=[], merit=np.inf, operand_values=[], success=False)

    specs = _build_specs(best_params)
    cache: dict[float, tuple[GaussBeam, GaussBeam]] = {}
    operand_values = []
    for op in operands:
        if op.z_mm not in cache:
            sx, sy = _build_xy_systems(beam, specs, op.z_mm)
            dz = max(1.0, op.z_mm / 100)
            cache[op.z_mm] = sx.trace(beam, dz=dz)[-1], sy.trace(beam, dz=dz)[-1]
        bx, by = cache[op.z_mm]
        operand_values.append(_eval_operand(op, bx, by))

    return OptResult(
        specs=specs,
        merit=_scalar(best_params),
        operand_values=operand_values,
        success=success,
    )


# ── Convenience ───────────────────────────────────────────────────────────────

def waist_operands(
    z_mm: float,
    wx_mm: float,
    wy_mm: float,
    size_weight: float = 1.0,
    waist_tol_x_mm: float | None = None,
    waist_tol_y_mm: float | None = None,
    curvature_weight_x: float | None = None,
    curvature_weight_y: float | None = None,
) -> list[Operand]:
    """Return operands targeting beam waists (wx, wy) at z_mm.

    Combines size operands (wx/wy) and waist condition operands (cvx/cvy = 0).

    Parameters
    ----------
    waist_tol_x_mm, waist_tol_y_mm : float, optional
        Acceptable waist displacement from z_mm (mm) for each axis independently.
        Converted to curvature_weight so that a waist displaced by waist_tol_*_mm
        contributes a residual equal to the corresponding size operand at its target.
        Ignored if the corresponding curvature_weight_* is given explicitly.
    curvature_weight_x, curvature_weight_y : float, optional
        Explicit weights for cvx/cvy. Override waist_tol_*_mm if given.
        Default (all None): same as size_weight.

    Notes
    -----
    The cvx/cvy residual is scaled internally as z_mm/R (dimensionless),
    which approximates Δz/z_mm for small waist displacements Δz.
    """
    def _cw(w_mm, tol_mm, explicit):
        if explicit is not None:
            return explicit
        if tol_mm is not None:
            curv_scale = tol_mm / z_mm
            return size_weight * (w_mm / curv_scale) ** 2 if curv_scale > 0 else 1.0
        return size_weight

    cwx = _cw(wx_mm, waist_tol_x_mm, curvature_weight_x)
    cwy = _cw(wy_mm, waist_tol_y_mm, curvature_weight_y)

    return [
        Operand('wx',  z_mm, wx_mm, size_weight),
        Operand('wy',  z_mm, wy_mm, size_weight),
        Operand('cvx', z_mm, 0.0,   cwx),
        Operand('cvy', z_mm, 0.0,   cwy),
    ]


# ── Minimum-configuration search ─────────────────────────────────────────────

def find_minimum_system(
    beam: GaussBeam,
    operands: list[Operand],
    lens_type: str | list[str] = 'spherical',
    max_lenses: int = 5,
    merit_threshold: float = 1e-3,
    z_bounds: tuple[float, float] | None = None,
    f_bounds: tuple[float, float] | None = None,
    f_abs_bounds: tuple[float, float] | None = None,
    min_lens_sep_mm: float = 0.0,
    algorithm: str = 'de',
    seed: int = 42,
    maxiter: int = 1000,
    popsize: int = 15,
    tol: float = 1e-12,
    verbose: bool = True,
) -> tuple[OptResult, int]:
    """Find the minimum number of lenses that satisfies the merit function.

    Tries n = 1, 2, ... up to max_lenses. Returns the first solution whose
    merit < merit_threshold.

    Parameters
    ----------
    lens_type : str or list of str
        Lens type(s) to repeat for each n. If a str, uses [lens_type] * n.
        If a list, its length must equal max_lenses; the first n elements
        are used for each trial (so ordering matters).
    merit_threshold : float
        Target merit value to consider the problem solved.
    verbose : bool
        Print progress per trial.

    Returns
    -------
    (OptResult, n_lenses)
        Best result found at the minimum n satisfying the threshold,
        or the best result overall if threshold was never reached.
    """
    best_result: OptResult | None = None
    best_n = 0

    for n in range(1, max_lenses + 1):
        if isinstance(lens_type, str):
            types = [lens_type] * n
        else:
            if n > len(lens_type):
                break
            types = list(lens_type[:n])

        result = optimize_astigmatic(
            beam, types, operands,
            z_bounds=z_bounds, f_bounds=f_bounds,
            f_abs_bounds=f_abs_bounds,
            min_lens_sep_mm=min_lens_sep_mm,
            algorithm=algorithm, seed=seed,
            maxiter=maxiter, popsize=popsize, tol=tol,
        )

        if verbose:
            vals = ', '.join(f'{v:.4g}' for v in result.operand_values)
            print(f'n={n}  merit={result.merit:.3e}  operands=[{vals}]')

        if best_result is None or result.merit < best_result.merit:
            best_result = result
            best_n = n

        if result.merit < merit_threshold:
            break

    return best_result, best_n


# ── Public helpers ────────────────────────────────────────────────────────────

def build_xy_systems(
    beam: GaussBeam,
    lens_specs: list[LensSpec],
    z_target_mm: float,
) -> tuple[OpticalSystem, OpticalSystem]:
    """Return (sys_x, sys_y) OpticalSystems for the given lens specs."""
    return _build_xy_systems(beam, lens_specs, z_target_mm)


# ── Simple 1-axis optimizer (unchanged) ──────────────────────────────────────

def find_lens_system(
    beam: GaussBeam,
    z_target_mm: float,
    w_target_mm: float,
    n_lenses: int = 1,
    z_bounds: tuple[float, float] | None = None,
    f_bounds: tuple[float, float] = (1.0, 1000.0),
    tol: float = 1e-10,
    seed: int = 42,
) -> OpticalSystem:
    """Find a rotationally symmetric lens system achieving w(z_target) = w_target."""
    try:
        from scipy.optimize import differential_evolution
    except ImportError:
        raise ImportError("scipy が必要です: pip install gbeampro[optimize]")

    if z_target_mm <= beam.z_mm:
        raise ValueError(
            f"z_target_mm ({z_target_mm}) は beam.z_mm ({beam.z_mm}) より大きくなければなりません"
        )

    z_lo, z_hi = z_bounds if z_bounds is not None else (beam.z_mm, z_target_mm)
    f_lo, f_hi = f_bounds
    trace_dz = max(0.5, (z_target_mm - beam.z_mm) / 200)

    def _build(params: np.ndarray) -> OpticalSystem:
        zs = np.sort(params[0::2])
        fs = params[1::2]
        sys = OpticalSystem()
        z_prev = beam.z_mm
        for z, f in zip(zs, fs):
            d = float(z - z_prev)
            if d > 0:
                sys.add(Propagation(d))
            sys.add(ThinLens(float(f)))
            z_prev = float(z)
        d_final = z_target_mm - z_prev
        if d_final > 0:
            sys.add(Propagation(d_final))
        return sys

    def _objective(params: np.ndarray) -> float:
        sys = _build(params)
        final = sys.trace(beam, dz=trace_dz)[-1]
        return (final.w_mm - w_target_mm) ** 2

    bounds = [(z_lo, z_hi), (f_lo, f_hi)] * n_lenses
    result = differential_evolution(
        _objective, bounds,
        tol=tol, seed=seed, maxiter=2000, popsize=15, polish=True, workers=1,
    )
    return _build(result.x)
