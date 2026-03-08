import marimo

__generated_with = "0.20.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Minimum Lens System Search — Progressive Constraint Tightening

    **Problem**: Find the minimum number of lenses to focus a circular Gaussian beam to a target waist.

    | Parameter | Value |
    |-----------|-------|
    | Wavelength λ | 1.064 µm |
    | Input waist w₀ | 2.0 mm |
    | Target position z | 200 mm |
    | Target waist w | 0.1 mm |
    | Waist tolerance Δz | ±10 mm |

    **Strategy**: Start with minimal constraints to find the theoretical minimum,
    then progressively tighten towards a realistic, manufacturable lens system.

    | Step | Added constraint | Physical meaning |
    |------|-----------------|------------------|
    | 1 | f ∈ (−10000, 10000) mm, no sep. | Theoretical baseline |
    | 2 | sep ≥ 20 mm, z ≤ 140 mm, Δf = 5 mm | Mechanical + catalog discretization |
    | 3 | + \|f\| ∈ [30, 1000] mm | Standard catalog focal length range |
    """)
    return


@app.cell
def _():
    # '%matplotlib inline' command supported automatically in marimo
    import numpy as np
    import matplotlib.pyplot as plt

    from gbeampro import GaussBeam, OpticalSystem, Propagation
    from gbeampro.optimize import waist_operands, find_minimum_system, build_xy_systems
    import gbeampro.plot as gplot
    import gbeampro.analysis as ga
    import gbeampro
    print('gbeampro version:', gbeampro.__version__)
    return (
        GaussBeam,
        build_xy_systems,
        find_minimum_system,
        ga,
        gplot,
        np,
        plt,
        waist_operands,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Problem Setup
    """)
    return


@app.cell
def _(GaussBeam, np, waist_operands):
    WL        = 1.064   # µm
    W0        = 2.0     # mm
    Z_TARGET  = 200.0   # mm
    W_TARGET  = 0.1     # mm
    WAIST_TOL = 10.0    # mm

    beam = GaussBeam.from_waist(wl_um=WL, w0_mm=W0)
    print(beam)

    z_R = np.pi * W0**2 / (WL * 1e-3)
    print(f'\nRayleigh range:       z_R = {z_R:.0f} mm')
    print(f'Beam size at z=200mm: w   = {beam.w_mm * np.sqrt(1+(Z_TARGET/z_R)**2)*1e3:.2f} µm  (no optics)')
    print(f'Target:               w   = {W_TARGET*1e3:.0f} µm')

    operands = waist_operands(
        z_mm=Z_TARGET, wx_mm=W_TARGET, wy_mm=W_TARGET,
        size_weight=1.0,
        waist_tol_x_mm=WAIST_TOL,
        waist_tol_y_mm=WAIST_TOL,
    )
    print(f'\nMerit function operands:')
    for op in operands:
        print(f'  {op.type:<4}  target={op.target:.4g}  weight={op.weight:.3g}')
    return W_TARGET, Z_TARGET, beam, operands


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Helper functions
    """)
    return


@app.cell
def _(W_TARGET, Z_TARGET, beam, build_xy_systems, ga, gplot, plt):
    def print_specs(result):
        print(f"{'#':>2}  {'z (mm)':>10}  {'f (mm)':>10}  {'|f| (mm)':>10}")
        print('   ' + '-' * 36)
        for i, s in enumerate(result.specs, 1):
            ftype = 'conv.' if s['f_mm'] > 0 else 'div.'
            print(f"{i:>2}  {s['z_mm']:>10.2f}  {s['f_mm']:>10.2f}  {abs(s['f_mm']):>10.2f}  ({ftype})")
        if len(result.specs) > 1:
            seps = [result.specs[i+1]['z_mm'] - result.specs[i]['z_mm']
                    for i in range(len(result.specs)-1)]
            print(f'  separations: {[f"{d:.1f} mm" for d in seps]}')

    def plot_result(result, n, title, z_end=Z_TARGET+50):
        sx, _ = build_xy_systems(beam, result.specs, z_end)
        traj = sx.trace(beam, dz=0.5)
        waists = ga.find_waists(traj)
        fig, ax = plt.subplots(figsize=(11, 3.5))
        gplot.plot_system(sx, traj, ax, label='beam')
        ax.axvline(Z_TARGET, color='k', ls=':', lw=1.2, label=f'z={Z_TARGET:.0f} mm')
        ax.scatter([Z_TARGET], [W_TARGET*1e3], color='r', zorder=6, s=80, label='target')
        ax.set_title(title)
        plt.tight_layout()
        print('Waists found:')
        for w in waists:
            dz = w.z_mm - Z_TARGET
            dw = (w.w_mm - W_TARGET) / W_TARGET * 100
            print(f'  z={w.z_mm:.1f} mm (Δz={dz:+.1f})  '
                  f'w₀={w.w_mm*1e3:.2f} µm ({dw:+.1f}%)  '
                  f'2z_R={2*ga.rayleigh_range(w):.1f} mm')
        return sx, traj

    return plot_result, print_specs


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Step 1: Theoretical Baseline

    No mechanical constraints — lenses can be placed anywhere in [0, 200 mm],
    any focal length is allowed, continuous f.
    This gives the **theoretical minimum** number of lenses.
    """)
    return


@app.cell
def _(beam, find_minimum_system, operands, print_specs):
    r1, n1 = find_minimum_system(
        beam, operands,
        lens_type='spherical',
        max_lenses=5,
        merit_threshold=1e-3,
        f_bounds=(-10000, 10000),
        min_lens_sep_mm=0.0,
        algorithm='de',
        verbose=True,
    )
    print(f'\n→ Minimum: {n1} lens(es),  merit = {r1.merit:.2e}')
    print_specs(r1)
    return n1, r1


@app.cell
def _(n1, plot_result, r1):
    sx1, traj1 = plot_result(r1, n1,
        f'Step 1: {n1} lens(es)  |  f ∈ (−10000, 10000) mm, no sep. constraint  |  merit={r1.merit:.2e}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Step 2: Mechanical + Catalog Discretization Constraints

    Added simultaneously:
    - **Lens separation ≥ 20 mm** — mechanical clearance
    - **Lens z ≤ 140 mm** — upstream optical path budget
    - **Δf = 5 mm** — focal length available only in 5 mm steps (catalog discretization)
    """)
    return


@app.cell
def _(beam, find_minimum_system, operands, print_specs):
    r2, n2 = find_minimum_system(
        beam, operands,
        lens_type='spherical',
        max_lenses=5,
        merit_threshold=1e-3,
        f_bounds=(-10000, 10000),
        f_step_mm=5.0,
        z_max_mm=140.0,
        min_lens_sep_mm=20.0,
        algorithm='de',
        verbose=True,
    )
    print(f'\n→ Minimum: {n2} lens(es),  merit = {r2.merit:.2e}')
    print_specs(r2)
    return n2, r2


@app.cell
def _(n2, plot_result, r2):
    sx2, traj2 = plot_result(r2, n2,
        f'Step 2: {n2} lens(es)  |  sep ≥ 20 mm, z ≤ 140 mm, Δf = 5 mm  |  merit={r2.merit:.2e}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Step 3: Add Focal Length Range Constraint

    Further restrict to **|f| ∈ [30, 1000] mm** — a range representative of
    standard catalog lenses (both converging and diverging allowed).

    All Step 2 constraints remain active.
    """)
    return


@app.cell
def _(beam, find_minimum_system, operands, print_specs):
    r3, n3 = find_minimum_system(
        beam, operands,
        lens_type='spherical',
        max_lenses=5,
        merit_threshold=1e-3,
        f_abs_bounds=(30, 1000),
        f_step_mm=5.0,
        z_max_mm=140.0,
        min_lens_sep_mm=20.0,
        algorithm='de',
        verbose=True,
    )
    print(f'\n→ Minimum: {n3} lens(es),  merit = {r3.merit:.2e}')
    print_specs(r3)
    return n3, r3


@app.cell
def _(n3, plot_result, r3):
    sx3, traj3 = plot_result(r3, n3,
        f'Step 3: {n3} lens(es)  |  + |f| ∈ [30,1000] mm  |  merit={r3.merit:.2e}')
    return (sx3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary
    """)
    return


@app.cell
def _(beam, n1, n2, n3, r1, r2, r3, sx3):
    print(f"{'Step':<6}  {'Constraints':<52}  {'N':>4}  {'merit':>10}")
    print('-' * 78)
    for step, desc, _n, _res in [('1', 'f free, no sep.', n1, r1), ('2', 'sep≥20mm, z≤140mm, Δf=5mm', n2, r2), ('3', 'sep≥20mm, z≤140mm, Δf=5mm, |f|∈[30,1000]mm', n3, r3)]:
        print(f'{step:<6}  {desc:<52}  {_n:>4}  {_res.merit:>10.2e}')
    print()
    print('Step 3 (most realistic) system:')
    print(sx3.summary(beam))
    return


@app.cell
def _(
    W_TARGET,
    Z_TARGET,
    beam,
    build_xy_systems,
    gplot,
    n1,
    n2,
    n3,
    plt,
    r1,
    r2,
    r3,
):
    Z_END = Z_TARGET + 50
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, (_res, _n, label) in zip(axes, [(r1, n1, 'Step 1: f free, no sep.'), (r2, n2, 'Step 2: sep≥20mm, z≤140mm, Δf=5mm'), (r3, n3, 'Step 3: + |f|∈[30,1000]mm')]):
        sx, _ = build_xy_systems(beam, _res.specs, Z_END)
        traj = sx.trace(beam, dz=0.5)
        gplot.plot_system(sx, traj, ax, label='beam')
        ax.axvline(Z_TARGET, color='k', ls=':', lw=1.2)
        ax.scatter([Z_TARGET], [W_TARGET * 1000.0], color='r', zorder=6, s=60, label='target')
        ax.set_title(f'{label}  →  {_n} lens(es),  merit={_res.merit:.1e}')
    plt.tight_layout()
    return


if __name__ == "__main__":
    app.run()
