# gbeampro

*gbeampro* is a small Python package for designing gaussian laser beam propagation and transformations (focusing and collimation, refraction at crystal interface, etc., for instance)  according to the ABCD law[^1].

[^1]: Kogelnik, Herwig. "Imaging of optical modes—resonators with internal lenses." Bell System Technical Journal 44.3 (1965): 455-494.

## Install
```shell
pip install gbeampro
```

## Example: Beam focusing into a slab of crystal

This is an example of how to compute the waist diameter inside the crystal and to find the **confocal parameter**[^2] (twice the Rayleigh range):

$$ 2 z_0 = \frac{2\pi w_0^2 n}{\lambda} $$

[^2]: Yariv, Amnon. *Quantum electronics, Third edition*. John Wiley & Sons, 1989.)

At first, import a fundamental (TEM00) Gaussian beam class from `gbeampro.beambase`.
```python
import ndispers as nd
import numpy as np
import matplotlib.pyplot as plt
from gbeampro.beambase import GaussBeam
```

Secondly, make an instance of `GaussBeam` object by giving the five parameters:
1. `wl_um` : Wavelength of the laser beam (unit: µm)
2. `n` : Refractive index of the medium in which the beam starts to propagate
3. `z_mm` : z coordinate of the wavefront. (unit: mm)
            Wavevector of the gaussian beam is assumed to be along to z axis.
4. `R_mm` :  Wavefront curvature radius. (unit: mm)
5. `w_mm` : Beam radius, defined as the 1/e**2 intensity half width. (unit: mm)

We will set a beam, named `b1`,  whose waist is located at z=0 and beam radius is 1.5 mm at z=0, as follows.
```python
# setup a GaussBeam instance
b1 = GaussBeam(wl_um=1.064, n=1.0, z_mm=0, R_mm=np.infty, w_mm=1.5)
```

By entering the beam instance `b1`, we can check the beam parameters at the present state.
```python
b1
```


Output:
```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=0.00000, R_mm=inf, w_mm=1.50000)
      q : 0.00000e+00 -6.64341e+03i
      theta : 2.257878e-01 mrad
```

Here a complex q-parameter `q` and beam divergence (half) angle `theta` are also printed.


Next, let's put a thin lens (focal length, f=150 mm) at z=0 by using a method,
```python
b1.thinlens(150)
```

Then the state of `b1` has been changed and the transformed beam parameters are shown:

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=0.00000, R_mm=-1.50000e+02, w_mm=1.50000)
      q : -1.49924e+02 -3.38509e+00i
      theta : 2.257878e-01 mrad
```

Let the beam propagate in the air (n=1.0) upto the entrance face of a crystal,


```python
b1.propagate(150 - 20*0.5 + 5)
```
(Here `+ 5` distance is trying to bring the wasit location at the center of the crystal.)

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=145.00000, R_mm=-7.25091e+00, w_mm=0.05977)
      q : -4.92357e+00 -3.38509e+00i
      theta : 5.666828e+00 mrad
```

The beam enters a slab of LBO crystal (a biaxial crystal). Here we need to compute the refractive index `n2` of the crystal. 


```python
Xtal = nd.media.crystals.LBO_Newlight_xy() # principal dielectric plane: xy
n2 = Xtal.n(1.064, 0, 149, pol='o') # for ordinary ray of wavelength 1.064 µm and temperature 149 degC.
b1.interface(n2)
```
(Note that this propagation method can not be applied for extraordinary rays in ansiotropic crystals.)

```python
    GaussBeam(wl_um=1.06400, n=1.604333, z_mm=145.00000, R_mm=-1.16329e+01, w_mm=0.05977)
      q : -7.89905e+00 -5.43082e+00i
      theta : 3.532224e+00 mrad
```

Let the beam propagate inside the crystal upto its end faceof the crystal, whose length is 20 mm.

```python
b1.propagate(20)
```

```python
    GaussBeam(wl_um=1.06400, n=1.604333, z_mm=165.00000, R_mm=1.45383e+01, w_mm=0.08270)
      q : 1.21010e+01 -5.43082e+00i
      theta : 2.552784e+00 mrad
```

At the end face of the crystal, the beam goes out into the air (with n=1.0).

```python
b1.interface(1.0)
```

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=165.00000, R_mm=9.06187e+00, w_mm=0.08270)
      q : 7.54267e+00 -3.38509e+00i
      theta : 4.095503e+00 mrad
```

We want the diverging beam to be collimated. So let the beam propagate in the air upto the second lens (f=200 mm).

```python
b1.propagate(200 - 20*0.5 + 5)
```
(Here `+ 5` distance is trying to collimate the beam as much as possible after the second lens.)

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=360.00000, R_mm=2.02599e+02, w_mm=2.02623)
      q : 2.02543e+02 -3.38509e+00i
      theta : 1.671490e-01 mrad
```

The second thin lens collimates the beam.

```python
b1.thinlens(200)
```

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=360.00000, R_mm=-1.55891e+04, w_mm=2.02623)
      q : -5.87433e+03 -7.55432e+03i
      theta : 1.671490e-01 mrad
```

Finally, let the beam propagate in air for an illustration purpose.

```python
b1.propagate(100)
```

```python
    GaussBeam(wl_um=1.06400, n=1.000000, z_mm=460.00000, R_mm=-1.56573e+04, w_mm=2.01330)
      q : -5.77433e+03 -7.55432e+03i
      theta : 1.682224e-01 mrad
```

Let's plot the beam trajectory in terms of its four parameters:

```python
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(8, 16), facecolor="white")
b1.plot_n(ax1)
b1.plot_w(ax2)
b1.plot_R(ax3)
b1.plot_theta(ax4)
```


    
![output_12_0](https://user-images.githubusercontent.com/88579896/179395036-151ba08b-37c3-4ecd-b0d0-ac2598329c24.png)

    
We want to find the beam waist inside the crystal. So use `search_BeamWaists` method. It searches all the sign-flip points of wavefront curvature (*R*), at which the beam becomes its waist.


```python
b1.search_BeamWaists()
```

```
    Beam waists in z range [0.000, 460.000000] mm
    --------------------------------------------------
    No.0:
      Waist location, z  : 0.0000 mm
      Refractive index, n  : 1.0000
      Waist spot diamter, 2*w0 : 3000.0 µm
      Confocal parameter (2*Rayleigh range) : 13286.81 mm
      2.84 * Confocal parameter : 37734.54 mm
    No.1:
      Waist location, z  : 152.9000 mm
      Refractive index, n  : 1.6043
      Waist spot diamter, 2*w0 : 67.7 µm
      Confocal parameter (2*Rayleigh range) : 10.86 mm
      2.84 * Confocal parameter : 30.85 mm
    No.2:
      Waist location, z  : 360.0000 mm
      Refractive index, n  : 1.0000
      Waist spot diamter, 2*w0 : 4052.5 µm
      Confocal parameter (2*Rayleigh range) : 24244.54 mm
      2.84 * Confocal parameter : 68854.49 mm
```

We can see the focusing condition inside the crystal at the above result No.1.

