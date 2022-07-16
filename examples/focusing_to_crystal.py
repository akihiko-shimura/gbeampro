"""
Beam focusing into a crystal block.
Compute the waist diameter inside the crystal and find the confocal parameter (collimated range).
"""
import ndispers as nd
import matplotlib.pyplot as plt
from gbeampro.beam_base import GaussBeam

wl = 1.064 #µm  lase wavelength
n1 = 1.0 #      refractive index (air)
R = 1.0e7 #mm    wavefront curvature radius (very large if collimated)
w = 2.0 #mm     beam radius (1/e**2 intensity half width)
f1 = 150 #mm        Focusing lens
f2 = 100 #mm        Collimating lens
Xtal = nd.media.crystals.LBO_Newlight_xy() # crystal object
n2 = Xtal.n(wl, 0, 149, pol='o') # refractive index of the crystal
L = 20.0 #mm    crystal length

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(8, 16))
beam1 = GaussBeam(wl_um=wl, n=n1, z_mm=0, R_mm=R, w_mm=w) # setup a GaussBeam instance
beam1.thinlens(f1) # focusing by lens 1
beam1.propagate(f1 - L*0.5 + 3) # propagation upto the entrance face of crystal (+3 is a small adjuctment trying to locate the beam waist at the center of the crystal.)
beam1.interface(n2) # Air-Crystal interface
beam1.propagate(L) # propagation inside crystal
beam1.interface(n1) # Crystal-Air interface
beam1.propagate(f2 - L*0.5 + 3) # # propagation in air
beam1.thinlens(f2) # collimation by lens 2
beam1.propagate(100) # propagation in air
beam1.plot_n(ax1)
beam1.plot_w(ax2)
beam1.plot_R(ax3)
beam1.plot_theta(ax4)
beam1.search_BeamWaists()

