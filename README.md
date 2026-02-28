# gbeampro

*gbeampro* is a small Python package for simulating Gaussian (TEM₀₀) laser beam propagation and transformations using the ABCD matrix method (q-parameter formalism).

## Installation

```bash
pip install gbeampro
```

## Quick Start

```python
from gbeampro import beambase

# Define a beam: 1064 nm, w=1 mm, flat wavefront
beam = beambase.GaussBeam(wl_um=1.064, n=1.0, w_mm=1.0)

# Propagate 100 mm, focus with f=50 mm lens, propagate into crystal (n=1.5)
beam.propagate(50).thinlens(50).interface(1.5).propagate(10)

# Find beam waists
beam.search_BeamWaists()
```

## API Reference

### `GaussBeam(wl_um, n, z_mm, R_mm, w_mm, label)`

Fundamental Gaussian beam object.

| Parameter | Description | Unit |
|-----------|-------------|------|
| `wl_um` | Wavelength | µm |
| `n` | Refractive index of medium | — |
| `z_mm` | z-coordinate of wavefront | mm |
| `R_mm` | Wavefront curvature radius | mm |
| `w_mm` | Beam radius (1/e² intensity half-width) | mm |
| `label` | Label string | — |

Key attributes: `.q` (complex q-parameter), `.theta` (divergence half-angle in rad).

### Transformation Methods

All methods mutate the beam in-place and return `self`, enabling method chaining.

| Method | Description |
|--------|-------------|
| `propagate(d_total, dz=0.01)` | Free-space propagation over distance `d_total` (mm), recorded at steps of `dz` |
| `thinlens(f)` | Thin lens with effective focal length `f` (mm) |
| `interface(n2)` | Refraction at a flat dielectric interface into medium with index `n2` |
| `interface_curved(n2, r)` | Refraction at a curved dielectric interface; `r > 0` convex, `r < 0` concave (mm) |
| `curved_mirror_tan(r, theta_deg)` | Reflection from a curved mirror, tangential (in-plane) component |
| `curved_mirror_sag(r, theta_deg)` | Reflection from a curved mirror, sagittal (out-of-plane) component |

### Analysis Methods

| Method | Description |
|--------|-------------|
| `search_BeamWaists(out=False)` | Locate beam waists (sign change of R). Prints results; returns `(n, z, w)` arrays if `out=True` |

### Plot Methods

| Method | Description |
|--------|-------------|
| `plot_w(ax)` | Beam radius w vs z |
| `plot_R(ax)` | Wavefront curvature radius R vs z |
| `plot_n(ax)` | Refractive index n vs z |
| `plot_theta(ax)` | Divergence half-angle θ vs z |

## Trajectory

Every transformation appends to `beam.traj`, a dict of lists:

| Key | Description |
|-----|-------------|
| `z_mm` | z-coordinate (mm) |
| `n` | Refractive index |
| `R_mm` | Wavefront curvature radius (mm) |
| `w_mm` | Beam radius (mm) |
| `q` | Complex q-parameter |
| `theta_rad` | Divergence half-angle (rad) |

## Examples

- [Beam focusing into a slab of crystal](gbeampro/examples/beam_focusing_into_crystal.ipynb)

## License

See [LICENSE](LICENSE).
