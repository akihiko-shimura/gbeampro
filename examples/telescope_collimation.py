"""
Beam shrinking and collimation with Galilean telescope.
"""
import matplotlib.pyplot as plt
from gbeampro.beam_base import GaussBeam

wl = 1.064 #µm  laser wavelength
n1 = 1.0 #      refractive index of medium (air)
f1 = 100 #mm    convex lens
f2 = -50 #mm    concave lens

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(8, 14))
beam1 = GaussBeam(wl_um=wl, n=n1, z_mm=0, R_mm=1.0e7, w_mm=0.5*4)
beam1.print_params(label="Initial, collimated")
beam1.thinlens(f1)
beam1.propagate(f1 + f2)
beam1.thinlens(f2)
beam1.propagate(400)
beam1.print_params(label='after collimation by lens2')
beam1.plot_n(ax1)
beam1.plot_w(ax2)
beam1.plot_R(ax3)
beam1.plot_theta(ax4)
beam1.search_BeamWaists()
# %%
