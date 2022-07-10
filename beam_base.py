""" 
Beam class
"""
import numpy as np
from gbeampro.helper import arg_signchange

pi = np.pi

class GaussBeam(object):
    __slots__ = ["WL", "n", "z", "R", "w", "n_hist", "z_hist", "R_hist", "w_hist", "q_hist", "theta_hist"]

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
        w_hist : list of float
            List of w values as the gaussian beam propagates.
        R_hist : list of float
            List of R values as the gaussian beam propagates.
        z_hist : list of float
            List of z values as the gaussian beam propagates.
        q_hist : list of complex
            List of q values as the gaussian beam propagates.
        theta_hist : list of float
            List of theta values as the gaussian beam propagates.

        """
        self.WL = float(wl_um)
        self.n = float(n)
        self.z = float(z_mm)
        self.R = float(R_mm)
        self.w = float(w_mm)

        self.n_hist = [float(n)]
        self.z_hist = [float(z_mm)]
        self.R_hist = [float(R_mm)]
        self.w_hist = [float(w_mm)]

        self.q_hist = [self.q()]
        self.theta_hist = [self.theta()]
    

    def q(self):
        """q parameter of the gaussian beam.

        q parameter is defined with three attributes, WL, R and w.
        
        Retruns
        -------
        complex
            q parameter of the gaussian beam.

        """
        return 1./(1j * self.WL*1e-3 / (pi * self.n * self.w**2) + 1./self.R)
    
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
    
    def _update_hist(self, n_new, z_new, R_new, w_new):
        self.n_hist.append(n_new)
        self.z_hist.append(z_new)
        self.R_hist.append(R_new)
        self.w_hist.append(w_new)
        self.q_hist.append(self.q())
        self.theta_hist.append(self.theta())
    
    def set(self, n_new, z_new, R_new, w_new):
        """Set new parameters after beam transformation.

        Note
        ----
        q will be deduced from newly set R and w. So there is no need to input q here.

        """
        self.n = n_new
        self.z = z_new
        self.R = R_new
        self.w = w_new
        self._update_hist(n_new, z_new, R_new, w_new)
    
    # Ray propagation methods
    def _propagate1(self, d):
        """Beam transformation method for propagation in homogeneous medium.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = d, C = 0, D = 1

        Parameters
        ----------
        d : float
            Length of the medium in which the beam propagates.

        """
        n_new = self.n
        z_new = self.z + d
        q_new = self.q() + d
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self.set(n_new, z_new, R_new, w_new)
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

        """
        for _ in range(int(d_total/dz)):
            self._propagate1(dz)
    
    def thinlens(self, f):
        """Beam transformation method of thin lens.

        This method converts the q-parameter of the gaussian beam according to the ABCD law:
        A = 1, B = 0, C = -1/f, D = 1

        Parameters
        ----------
        f : float
            Effective focal length of the thin lens.

        """
        n_new = self.n
        z_new = self.z
        q_new = self.q() / (- self.q() / f + 1.0)
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self.set(n_new, z_new, R_new, w_new)
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

        """
        n1 = self.n
        z_new = self.z
        q_new = self.q() * (n2/n1)
        self.n = n2
        R_new = self.make_R_from_q(q_new)
        w_new = self.make_w_from_q(q_new)
        self.set(n2, z_new, R_new, w_new)
        return self
    
    def Amplitude(self, r):
        """Field amplitude of a gaussian beam.

        This method calcurates the E-field amplitude as a function of transverse radius from the peak.
        At beam waist, q0 = -1j*z0.
        
        """
        k = 2*pi/(self.WL*1e-3)
        q = self.q()
        z = self.z
        return np.exp( - np.log(1.0 + z/q) + 1j*k*r**2/(2*(z + q)) )
    
    def plot_w(self, ax, marker=".", size=3):
        """Plot of (z, w) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        w_hist = np.asarray(self.w_hist)
        z_hist = np.asarray(self.z_hist)
        ax.scatter(z_hist, w_hist * 1e3, s=size, marker=marker)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$w$ (µm)")
        ax.set_ylim(0, np.max(w_hist) * 1e3 * 1.05)
        # ax.axhline(y=0, lw=0.5, color='k')
        # ax.legend(bbox_to_anchor=(1.02, 1))
    
    def plot_R(self, ax):
        """Plot of (z, R) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        R_hist = np.asarray(self.R_hist)
        z_hist = np.asarray(self.z_hist)
        ax.set_yscale("symlog")
        # ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
        ax.plot(z_hist, R_hist)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$R$ (mm)")
        ax.axhline(y=0, lw=0.5, color='k')
    
    def plot_n(self, ax):
        """Plot of (z, n) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        n_hist = np.asarray(self.n_hist)
        z_hist = np.asarray(self.z_hist)
        ax.plot(z_hist, n_hist)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel("$n$")
    
    def plot_theta(self, ax):
        """Plot of (z, theta) points to a given figure axis.

        Parameters
        ----------
        ax : matplotlib Axes object
   
        """
        z_hist = np.asarray(self.z_hist)
        theta_hist = np.asarray(self.theta_hist)
        ax.semilogy(z_hist, theta_hist*1e3)
        ax.set_xlabel("$z$ (mm)")
        ax.set_ylabel(r"$\theta$ (mrad)")
    
    def print_params(self, label=""):
        """Report the present beam parameters.
        """
        print("\nGaussBeam (%s)" % label)
        print("----------------------------------")
        print("  z : %.3f mm" % self.z)
        print("  n : %.4f" % self.n)
        print("  q : {c.real:.3f} {c.imag:+.3f}i".format(c=self.q()))
        print("  w : %.3f mm" % self.w)
        print("  R : %.3f m" %(self.R*1e-3))
        print("  theta : %.3f mrad" %(self.theta()*1e3))
        print("----------------------------------")
    
    def search_BeamWaists(self):
        """Search for and report the beam waists.

        Here beam waists are defined as the points where the sign of wavefront curvature flips.
        """
        idx_waists = arg_signchange(np.asarray(self.R_hist))
        n_waists = np.asarray(self.n_hist)[idx_waists]
        z_waists = np.asarray(self.z_hist)[idx_waists]
        w_waists = np.asarray(self.w_hist)[idx_waists]

        print("\nBeam waists in z-range [%.2f, %.2f] mm" %(self.z_hist[0], self.z_hist[-1]))
        print("--------------------------------------------------")
        for i in range(idx_waists[0].size):
            print("No.%d:" % i)
            print("  Waist location, z  : %.4f mm" % z_waists[i])
            print("  Refractive index, n  : %.4f" % n_waists[i])
            print("  Waist spot diamter, 2*w0 : %.1f µm" %(2 * w_waists[i]*1e3))
            print("  Confocal parameter (2*Rayleigh range) : %.2f mm" %(2*pi * w_waists[i]**2 * n_waists[i] / (self.WL*1e-3)))
            print("  2.84 * Confocal parameter : %.2f mm" %(2.84 * 2*pi* (w_waists[i])**2 * n_waists[i] / (self.WL*1e-3)))

