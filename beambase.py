""" 
Beam class
"""
import numpy as np
from gbeampro.helper import arg_signchange

pi = np.pi

class GaussBeam(object):
    __slots__ = ["WL", "n", "z", "R", "w", "traj"]

    def __init__(self, wl_um=1.064, n=1.0, z_mm=0, R_mm=np.infty, w_mm=1.0):
        """Fundamental (TEM00) gaussian beam class.

        Parameters
        ----------
        wl_um : float
            Wavelength of the laser beam. (unit: µm)
        n : float
            Refractive index of medium.
        z_mm : float
            z coordinate of the wavefront. (unit: mm)
        R_mm : float
            Wavefront curvature radius. (unit: mm)
        w : flaot
            Beam radius. Defined as the 1/e**2 intensity half width. (unit: mm)
        
        Attributes
        ----------
        WL : float
            Wavelength of the laser beam. (unit: µm)
        n : float
            Refractive index of medium.
        z : float
            z coordinate of the wavefront. (unit: mm)
            Wavevector of the gaussian beam is assumed to be along to z axis.
        R : float
            Wavefront curvature radius. (unit: mm)
        w : flaot
            Beam radius. Defined as the 1/e**2 intensity half width. (unit: mm)
        traj : dict of list
            Beam trajectory containing all the tracks of the attributes' values.
            traj['z_mm'] : track of z coordinate.
            traj['n'] : track of refractive index.
            traj['R_mm'] : track of Wavefront curvature radius.
            traj['w_mm'] : track of beam radius.
            traj['q'] : track of q-parameter.
            traj['thtea'] : track of beam divergence half angle.

        """
        self.WL = float(wl_um)
        self.n = float(n)
        self.z = float(z_mm)
        self.R = float(R_mm)
        self.w = float(w_mm)

        self.traj = {'n':[float(n)], 
                     'z_mm': [float(z_mm)], 
                     'R_mm': [float(R_mm)], 
                     'w_mm': [float(w_mm)],
                     'q': [self.q],
                     'theta_rad': [self.theta]}
    
    @property
    def q(self):
        """q parameter of the gaussian beam.

        q parameter is defined with three attributes, WL, R and w.
        
        Retruns
        -------
        complex
            q parameter of the gaussian beam.

        """
        return 1./(1j * self.WL*1e-3 / (pi * self.n * self.w**2) + 1./self.R)
    
    @property
    def theta(self):
        """Beam divergence half angle.

        Returns
        -------
        float
            Beam divergence half angle (unit: radians)

        """
        return np.arctan(self.WL*1e-3 / (pi * self.n * self.w))

    def make_w_from_q(self, q):
        """Method to deduce `w` from a given q-parameter.

        Parameters
        ----------
        q : complex
            q parameter of a gaussian beam.
        n : float
            refactive index of a medium in which the beam propagates.

        Returns
        -------
        float
            w, beam radius.

        Note
        ----
        Here the imaginary part of q inverse (1/q) is forced to be positive so that `w_new` is always real number.

        """
        qinv_imag = np.imag(1./q)
        w_new = np.sqrt(self.WL*1e-3 / (pi * self.n * np.abs(qinv_imag)))
        return w_new
    
    def make_R_from_q(self, q):
        """Method to deduce `R` from a given q-parameter.

        Parameters
        ----------
        q : complex
            q parameter of a gaussian beam.
        n : float
            refactive index of a medium in which the beam propagates.

        Returns
        -------
        float
            R, wavefront curvature radius.

        """
        qinv_real = np.real(1./q)
        try:
            R_new = 1./qinv_real
        except ZeroDivisionError:
            R_new = np.inf
        
        return R_new

    def _add_to_traj(self, **kwargs):
        # update key arguments
        for key, val in kwargs.items():
            self.traj[key].append(val)
        # update properties
        self.traj['q'].append(self.q)
        self.traj['theta_rad'].append(self.theta)
    
    def _set(self, n_new, z_new, R_new, w_new):
        """Set new parameters after beam transformation.

        Note
        ----
        q and theta will be deduced from newly set R and w. So there is no need to input them.

        """
        self.n = n_new
        self.z = z_new
        self.R = R_new
        self.w = w_new
        self._add_to_traj(n=n_new, z_mm=z_new, R_mm=R_new, w_mm=w_new)
    
    def __repr__(self):
        p1 = "{}(wl_um={:.5f}, n={:.6f}, z_mm={:.5f}, R_mm={:.5e}, w_mm={:.5f})".format(
            self.__class__.__name__,
            self.WL, self.n, self.z, self.R, self.w
        )
        p2 = "  q : {c.real:.5e} {c.imag:+.5e}i".format(c=self.q)
        p3 = "  theta : {:.6e} mrad".format(self.theta*1e3)
        return '\n'.join([p1, p2, p3])
    
    def _propagate1(self, d):
        """Beam transformation method for propagation in homogeneous medium.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = d, C = 0, D = 1

        Parameters
        ----------
        d : float
            Length of the medium in which the beam propagates.

        Returns
        -------
        object, transformed GaussBeam object.

        """
        n_new = self.n
        z_new = self.z + d
        q_new = self.q + d
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self._set(n_new, z_new, R_new, w_new)
        return self

    def propagate(self, d_total, dz=0.01):
        """Iterative beam propagation in homegeneous medium with a given z increment.

        This method iteratively apply `homo_medium` method by total distance `D`.
        This is usefull to simulate the caustic curve along z.

        Parameters
        ----------
        d_total : float
            Total distance of propagation in a homogeneous medium.
        dz : float
            Finite step size of the iterative propagation along z.
            Default: 0.01
        
        Returns
        -------
        object, transformed GaussBeam object.

        """
        for _ in range(int(d_total/dz)):
            self._propagate1(dz)
        return self
    
    def thinlens(self, f):
        """Beam transformation method of thin lens.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = 0, C = -1/f, D = 1

        Parameters
        ----------
        f : float
            Effective focal length of the thin lens.

        Returns
        -------
        object, transformed GaussBeam object.

        """
        n_new = self.n
        z_new = self.z
        q_new = self.q / (- self.q / f + 1.0)
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self._set(n_new, z_new, R_new, w_new)
        return self

    def interface(self, n2):
        """Beam transformation method of refraction at a dielectric interface.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = 0, C = 0, D = n1/n2

        Parameters
        ----------
        n2 : float
            Refractive index of the next medium.
        
        Note
        ----
        `self.n` need to be changed from n1 to n2 before make R and w in order to make w continuous 
        at the interface.

        Returns
        -------
        object, transformed GaussBeam object.

        """
        n1 = self.n
        z_new = self.z
        q_new = self.q * (n2/n1)
        self.n = n2
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self._set(n2, z_new, R_new, w_new)
        return self
    
    def interface_curved(self, n2, r):
        """Beam transformation method of refraction at a curved dielectric interface.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = 0, C = (n2 - n1)/(n2*r), D = n1/n2

        Parameters
        ----------
        n2 : float
            Refractive index of the next medium.
        r  : float
            Curvature radius of the interface (unit: mm). r > 0, convex; r < 0, concave.
            Focal length,  f = r/2.
        
        Returns
        -------
        object, transformed GaussBeam object.

        """
        n1 = self.n
        z_new = self.z
        q_new = self.q / (self.q * (n2 - n2)/(n2*r) + n1/n2)
        self.n = n2
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self._set(n2, z_new, R_new, w_new)
        return self
    
    def Amplitude(self, r):
        """Field amplitude of a gaussian beam.

        This method calcurates the E-field amplitude as a function of transverse radius from the peak.
        At beam waist, q0 = -1j*z0.
        
        """
        k = 2*pi/(self.WL*1e-3)
        q = self.q
        z = self.z
        return np.exp( - np.log(1.0 + z/q) + 1j*k*r**2/(2*(z + q)) )
    
    def plot_w(self, ax, marker=".", size=3):
        """Plot of (z, w) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        w_track = np.asarray(self.traj['w_mm'])
        z_track = np.asarray(self.traj['z_mm'])
        ax.scatter(z_track, w_track * 1e3, s=size, marker=marker)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$w$ (µm)")
        ax.set_ylim(0, np.max(w_track) * 1e3 * 1.05)
    
    def plot_R(self, ax):
        """Plot of (z, R) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        z_track = np.asarray(self.traj['z_mm'])
        R_track = np.asarray(self.traj['R_mm'])
        ax.set_yscale("symlog")
        ax.plot(z_track, R_track)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$R$ (mm)")
        ax.axhline(y=0, lw=0.5, color='k')
    
    def plot_n(self, ax):
        """Plot of (z, n) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        z_track = np.asarray(self.traj['z_mm'])
        n_track = np.asarray(self.traj['n'])
        ax.plot(z_track, n_track)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$n$")
    
    def plot_theta(self, ax):
        """Plot of (z, theta) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        z_track = np.asarray(self.traj['z_mm'])
        theta_track = np.asarray(self.traj['theta_rad'])
        ax.semilogy(z_track, theta_track*1e3)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel(r"$\theta$ (mrad)")
    
    def search_BeamWaists(self):
        """Search for and report the beam waists.

        Here beam waists are defined as the points where the sign of wavefront curvature flips.
        """
        idx_waists = arg_signchange(np.asarray(self.traj['R_mm']))
        n_waists = np.asarray(self.traj['n'])[idx_waists]
        z_waists = np.asarray(self.traj['z_mm'])[idx_waists]
        w_waists = np.asarray(self.traj['w_mm'])[idx_waists]

        print("\nBeam waists in z range [{:.3f}, {:3f}] mm".format(self.traj['z_mm'][0], self.traj['z_mm'][-1]))
        print("--------------------------------------------------")
        for i in range(idx_waists[0].size):
            print("No.%d:" % i)
            print("  Waist location, z  : {:.4f} mm".format(z_waists[i]))
            print("  Refractive index, n  : {:.4f}".format(n_waists[i]))
            print("  Waist spot diamter, 2*w0 : {:.1f} µm".format((2 * w_waists[i]*1e3)))
            print("  Confocal parameter (2*Rayleigh range) : {:.2f} mm".format(2*pi * w_waists[i]**2 * n_waists[i] / (self.WL*1e-3)))
            print("  2.84 * Confocal parameter : {:.2f} mm".format( 2.84 * 2*pi* (w_waists[i])**2 * n_waists[i] / (self.WL*1e-3) ))

