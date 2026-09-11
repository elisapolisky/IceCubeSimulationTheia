from theia.material import HenyeyGreensteinPhaseFunction
from theia.material import DispersionFreeMedium

class LayerMediumModel_old(
    HenyeyGreensteinPhaseFunction,
    DispersionFreeMedium
):
    def __init__(
        self,
        mu_a,
        mu_s,
        *,
        g=0.9,
        n=1.15,
        ng= 1.35,
        name= "layer"

    ):
        super().__init__(
            n=n,
            ng=ng,
            mu_a=mu_a,
            mu_s = mu_s,
            g=g,
            name=name,
        )

from theia.material import HenyeyGreensteinPhaseFunction, MediumModel
from theia.material import medium_property
import numpy as np
import theia.units as u

import numpy as np
from numpy.typing import ArrayLike, NDArray

import theia.units as u
from theia.material import MediumModel
from theia.material import DispersionFreeMedium

import numpy as np
import theia.units as u
from theia.material import MediumModel
from theia.material import medium_property


from theia.material import NumericalPhaseSamplingMixin
from theia.material import medium_property

class HGSAMPhaseFunction(NumericalPhaseSamplingMixin):

    def __init__(
        self,
        g=0.9,
        fSL=0.35,
        *,
        name="HGSAM",
        **kwargs,
    ):
        super().__init__(name=name, **kwargs)

        self.g = g
        self.fSL = fSL

    def phase_function(self, cos_theta):

        x = np.asarray(cos_theta, dtype=float)

        g = self.g
        f = self.fSL

        hg = (
            (1.0 - g**2)
            / (
                4.0 * np.pi
                * (1.0 + g**2 - 2.0 * g * x) ** 1.5
            )
        )

        sam_alpha = 2.0 * g / (1.0 - g)

        sam = (
            (1.0 + sam_alpha)
            / (4.0 * np.pi)
            * ((1.0 + x) / 2.0) ** sam_alpha
        )

        return (1.0 - f) * hg + f * sam

    @medium_property(range=(-1.0, 1.0))
    def log_phase_function(self, cos_theta):

        return np.log(
            self.phase_function(cos_theta)
        )
class LayerMediumModel(HGSAMPhaseFunction):

    def __init__(
        self,
        be400,
        adust400,
        delta_tau,
        *,
        alpha=0.898608505726,
        kappa=1.084106802940,
        A=6954.090332031250,
        B=6617.754394531250,
        g=0.9,
        fSL=0.35,
        n=1.15,
        ng=1.35,
        name="layer",
    ):
        super().__init__(
            g=g,
            fSL = fSL,
            name=name,
        )

        self.be400 = be400
        self.adust400 = adust400
        self.delta_tau = delta_tau

        self.alpha = alpha
        self.kappa = kappa
        self.A = A
        self.B = B

        self.n = n
        self.ng = ng

    @medium_property
    def refractive_index(self, wavelength):
        return np.ones_like(wavelength) * self.n

    @medium_property
    def group_velocity(self, wavelength):
        return np.ones_like(wavelength) / self.ng * u.c

    @medium_property
    def scattering_coef(self, wavelength):
        wavelength = np.asarray(wavelength)
    
        be = (
            self.be400
            * (wavelength / (400 * u.nm))**(-self.alpha)
        )
    
        return be / (1.0 - self.g)
    
    
    @medium_property
    def absorption_coef(self, wavelength):
        wavelength = np.asarray(wavelength)
    
        a_dust = (
            self.adust400
            * (wavelength / (400 * u.nm))**(-self.kappa)
        )
    
        a_ice = (
            self.A
            * np.exp(-self.B / (wavelength / u.nm))
            * (1.0 + 0.01 * self.delta_tau)
        )
    
        return a_dust + a_ice