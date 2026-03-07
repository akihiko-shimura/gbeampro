# gbeampro

*gbeampro* is a small Python package for simulating Gaussian (TEM₀₀) laser beam propagation and transformations using the ABCD matrix method (q-parameter formalism).

## Installation

```bash
pip install gbeampro
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

sys.trace(beam, dz=0.01)   # -> list[GaussBeam], full caustic trajectory
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

gplot.plot_caustic(trajectory, ax)              # w vs z
gplot.plot_system(sys, trajectory, ax,
                  label="beam", beam_kw={})     # caustic + element symbols
```

Multiple beams can be overlaid by calling `plot_system` on the same `ax`; each label gets a distinct color from the matplotlib color cycle.

## Display Example

```
OpticalSystem
=========================================================
   #  Type                Parameters                z (mm)
---------------------------------------------------------
   0  --- input ---                                  0.000
   1  Propagation         d =  100.000 mm          100.000
   2  ThinLens            f =   50.000 mm          100.000
   3  Propagation         d =  100.000 mm          200.000
=========================================================
Total length: 200.000 mm  |  3 elements

OpticalSystem trace  [wl=1.064 um]
========================================================================
   #  Type                   z (mm)     w (um)      R (mm)   th (urad)
------------------------------------------------------------------------
   0  --- input ---           0.000    1000.00         inf      338.68
   1  Propagation           100.000    1000.57   8.728e+04      338.49
   2  ThinLens              100.000    1000.57  -5.003e+01      338.49
   3  Propagation           200.000    1000.00   5.000e+01      338.68
========================================================================
Beam waists:  z=150.500 mm (2w0=39.0 um)
```

## Examples

- [All elements test](gbeampro/examples/test_all_elements.ipynb)
- [plot_system test](gbeampro/examples/test_plot_system.ipynb)
- [Beam focusing into a crystal (v1 API)](gbeampro/examples/beam_focusing_into_crystal.ipynb)

## License

See [LICENSE](LICENSE).
