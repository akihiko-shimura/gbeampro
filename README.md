# gbeampro

*gbeampro* is a Python package for Gaussian (TEM₀₀) laser beam propagation and optimization using the ABCD matrix method (q-parameter formalism). It supports beam tracing through arbitrary optical systems and Zemax-inspired merit-function optimization for astigmatic beam shaping.

![demo](assets/demo.png)

## Installation

```bash
pip install gbeampro           # core
pip install gbeampro[optimize] # + scipy for optimization
```

## Quick Start

```python
from gbeampro import GaussBeam, Propagation, ThinLens, Interface, OpticalSystem

# Define a beam at its waist: 1064 nm, w₀=1 mm
beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=1.0)

# Build an optical system
sys = (OpticalSystem()
       .add(Propagation(100))
       .add(ThinLens(f_mm=50))
       .add(Interface(n1=1.0, n2=1.5))
       .add(Propagation(30)))

# Print system layout and beam state at each element
print(sys)
print(sys.summary(beam))

# Trace the full caustic
traj = sys.trace(beam, dz=0.5)
```

## API Reference

### `GaussBeam`

Immutable Gaussian beam value object (`frozen dataclass`).

| Parameter | Description | Unit |
|-----------|-------------|------|
| `wl_um` | Wavelength | µm |
| `n` | Refractive index | — |
| `z_mm` | z-coordinate of wavefront | mm |
| `R_mm` | Wavefront curvature radius (`inf` at waist) | mm |
| `w_mm` | Beam radius (1/e² intensity half-width) | mm |

Key properties: `.q` (complex q-parameter), `.theta` (divergence half-angle in rad).

**Constructors**

```python
GaussBeam.from_waist(wl_um, w0_mm, z_mm=0.0, n=1.0)  # from beam waist
GaussBeam.from_q(wl_um, n, q, z_mm=0.0)               # from complex q-parameter
```

### Optical Elements

Each element implements `apply(beam) -> GaussBeam` based on its ABCD matrix.

| Class | Parameters | Description |
|-------|-----------|-------------|
| `Propagation(d_mm)` | `d` — distance (mm) | Free-space propagation |
| `ThinLens(f_mm)` | `f` — focal length (mm) | Thin lens |
| `Interface(n1, n2)` | `n1`, `n2` — refractive indices | Flat dielectric interface |
| `InterfaceCurved(n1, n2, r_mm)` | `r > 0` convex, `r < 0` concave (mm) | Curved dielectric interface |
| `CurvedMirrorTan(r_mm, theta_deg)` | `r` — radius (mm), `θ` — angle of incidence (deg) | Curved mirror, tangential |
| `CurvedMirrorSag(r_mm, theta_deg)` | `r` — radius (mm), `θ` — angle of incidence (deg) | Curved mirror, sagittal |

Custom elements can be added by subclassing `Element` and implementing the `matrix` property.

### `OpticalSystem`

```python
sys = OpticalSystem().add(element1).add(element2)  # fluent API

sys.trace(beam, dz=0.5)    # -> list[GaussBeam], full caustic trajectory
str(sys)                   # element layout table
sys.summary(beam)          # beam state at each element + waist report
```

### Analysis (`gbeampro.analysis`)

```python
from gbeampro.analysis import find_waists, rayleigh_range, confocal_parameter

find_waists(trajectory)      # -> list[GaussBeam] at waist locations
rayleigh_range(beam)         # -> float, z_R (mm)
confocal_parameter(beam)     # -> float, 2*z_R (mm)
```

### Plot (`gbeampro.plot`)

```python
import gbeampro.plot as gplot

gplot.plot_system(sys, trajectory, ax, label="beam")
```

Multiple beams can be overlaid by calling `plot_system` on the same `ax`; each label gets a distinct color from the matplotlib color cycle.

### Optimization (`gbeampro.optimize`)

Requires `pip install gbeampro[optimize]`.

> Full API reference: [docs/optimize_api.md](docs/optimize_api.md)

Lens system optimization inspired by Zemax OpticStudio's merit function approach.
The merit function is defined as a list of **operands** — each specifying a beam
property, an evaluation position, a target value, and a weight.

**Operand types**

| Type | Quantity | Unit |
|------|----------|------|
| `wx` / `wy` | Beam radius at `z_mm` for x / y axis | mm |
| `cvx` / `cvy` | Wavefront curvature 1/R at `z_mm`; target=0 → waist | mm⁻¹ |
| `thx` / `thy` | Half-divergence angle at `z_mm` | mrad |

**Algorithms**

| `algorithm=` | Description |
|---|---|
| `'de'` | Differential evolution — global, robust |
| `'lm'` | Levenberg-Marquardt / TRF — local, fast |
| `'hammer'` | DE global search → TRF polish (like Zemax Hammer) |

**Example 1: astigmatic beam shaping**

Three-lens configuration (spherical + cyl_x + cyl_y) focusing a circular beam to an astigmatic waist.
`f_abs_bounds` restricts focal lengths to catalog range; `f_step_mm` enforces discrete steps.

```python
from gbeampro import GaussBeam
from gbeampro.optimize import waist_operands, optimize_astigmatic, build_xy_systems

# input beam: 800 nm wavelength, 1.0 mm waist radius
beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=1.0)

# merit function: target wx=120 µm, wy=400 µm at z=500 mm with waist condition (cvx=cvy=0)
operands = waist_operands(
    z_mm=500, wx_mm=0.12, wy_mm=0.40,
    size_weight=1.0,
    curvature_weight_x=1.0,   # weight for x-axis waist condition
    curvature_weight_y=1.0,   # weight for y-axis waist condition
)

result = optimize_astigmatic(
    beam,
    lens_types=['spherical', 'cyl_x', 'cyl_y'],  # spherical + two cylindrical lenses
    operands=operands,
    f_abs_bounds=(30, 1000),   # restrict |f| to [30, 1000] mm (catalog range)
    f_step_mm=5.0,             # round focal lengths to 5 mm steps
    z_max_mm=480.0,            # place all lenses before z=480 mm
    min_lens_sep_mm=20.0,      # minimum separation between adjacent lenses
    algorithm='de',            # differential evolution — global search
    maxiter=2000, popsize=20, seed=42,
)

print(result.specs)    # [{'type': ..., 'z_mm': ..., 'f_mm': ...}, ...]
print(f'merit={result.merit:.2e}  elapsed={result.elapsed_s:.1f}s')
```

**Example 2: find minimum lens count automatically**

Tries n=1, 2, ... lenses in sequence and stops as soon as `merit < merit_threshold`.

```python
from gbeampro.optimize import find_minimum_system

result, n = find_minimum_system(
    beam, operands,
    lens_type='spherical',     # repeat this lens type n times
    max_lenses=5,              # upper limit on lens count to try
    merit_threshold=1e-3,      # stop when merit falls below this value
    f_abs_bounds=(30, 1000),   # restrict |f| to [30, 1000] mm
    f_step_mm=5.0,             # round focal lengths to 5 mm steps
    z_max_mm=140.0,            # place all lenses before z=140 mm
    min_lens_sep_mm=20.0,      # minimum separation between adjacent lenses
    algorithm='de',
)
print(f'Minimum: {n} lens(es),  merit={result.merit:.2e}')
```

## Examples

- [All elements test](examples/test_all_elements.ipynb)
- [plot_system test](examples/test_plot_system.ipynb)
- [Beam focusing into a crystal](examples/beam_focusing_into_crystal.ipynb)
- [Minimum lens system search](examples/minimum_lens_system.ipynb)
- [Astigmatic beam shaping](examples/astigmatic_beam_shaping.ipynb)

## License

See [LICENSE](LICENSE).
