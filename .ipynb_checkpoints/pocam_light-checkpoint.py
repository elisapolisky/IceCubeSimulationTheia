from __future__ import annotations
from theia.lookup import Table
import numpy as np
import scipy.constants as consts
from scipy.integrate import quad
from scipy.stats.sampling import NumericalInversePolynomial
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats.sampling import NumericalInversePolynomial
from theia.light import LightSource, WavelengthSource
import hephaistos as hp
from hephaistos.pipeline import PipelineStage, SourceCodeMixin
from hephaistos.queue import IOQueue

from ctypes import Structure, c_float, c_int64, c_uint32
from hephaistos.glsl import buffer_reference, vec2, vec3
from theia.util import createCType

from theia.compiler import compileShader, createPreamble, loadShader
from theia.material import MaterialStore
from theia.random import RNG
from theia.ray import RayModel
import theia.units as u

from warnings import warn

from collections.abc import Callable
from numpy.typing import ArrayLike, NDArray
from theia.light import LightSource, WavelengthSource, FunctionWavelengthSource

_pulse_FB_WIDTH15 = np.array(
          [[  0.00000000e+00,   1.00000000e+00,   2.00000000e+00,
              3.00000000e+00,   4.00000000e+00,   5.00000000e+00,
              6.00000000e+00,   7.00000000e+00,   8.00000000e+00,
              9.00000000e+00,   1.00000000e+01,   1.10000000e+01,
              1.20000000e+01,   1.30000000e+01,   1.40000000e+01,
              1.50000000e+01,   1.60000000e+01,   1.70000000e+01,
              1.80000000e+01,   1.90000000e+01,   2.00000000e+01,
              2.10000000e+01,   2.20000000e+01,   2.30000000e+01,
              2.40000000e+01,   2.50000000e+01,   2.60000000e+01,
              2.70000000e+01,   2.80000000e+01,   2.90000000e+01,
              3.00000000e+01,   3.10000000e+01,   3.20000000e+01,
              3.30000000e+01,   3.40000000e+01,   3.50000000e+01,
              3.60000000e+01,   3.70000000e+01,   3.80000000e+01,
              3.90000000e+01,   4.00000000e+01,   4.10000000e+01,
              4.20000000e+01,   4.30000000e+01,   4.40000000e+01,
              4.50000000e+01,   4.60000000e+01,   4.70000000e+01,
              4.80000000e+01,   4.90000000e+01,   5.00000000e+01],
           [  1.18000000e-03,   2.76900000e-02,   1.25170000e-01,
              2.14840000e-01,   3.20890000e-01,   4.32390000e-01,
              4.64370000e-01,   5.00230000e-01,   4.31610000e-01,
              3.16210000e-01,   2.29650000e-01,   1.37640000e-01,
              8.77400000e-02,   7.21400000e-02,   5.96600000e-02,
              4.79700000e-02,   4.09500000e-02,   2.92500000e-02,
              3.08100000e-02,   2.84700000e-02,   2.61300000e-02,
              1.83400000e-02,   1.83400000e-02,   1.99000000e-02,
              1.28800000e-02,   1.28800000e-02,   1.28800000e-02,
              1.60000000e-02,   1.44400000e-02,   1.67800000e-02,
              7.42000000e-03,   6.64000000e-03,   9.76000000e-03,
              1.13200000e-02,   7.42000000e-03,   9.76000000e-03,
              4.30000000e-03,   5.86000000e-03,   7.42000000e-03,
              4.30000000e-03,   8.20000000e-03,   5.86000000e-03,
              3.52000000e-03,   1.96000000e-03,   2.74000000e-03,
              4.30000000e-03,   5.08000000e-03,   2.74000000e-03,
              3.52000000e-03,   4.30000000e-03,   2.74000000e-03]])
    #adjust zero offset and re-scale to one
_pulse_FB_WIDTH15[1]=(_pulse_FB_WIDTH15[1]-0.00118)/0.49905
_pulse_narrow = interp1d(_pulse_FB_WIDTH15[0], _pulse_FB_WIDTH15[1], kind='linear', bounds_error=False,
                         fill_value=0.)

def _pulse_rising_edge(x, width):
    templateWidth=7.
    scaledX = templateWidth*x/width

    scaledX = np.where(scaledX>templateWidth,templateWidth,scaledX)
    scaledX = np.where(scaledX<0.,0.,scaledX)

    return _pulse_narrow(scaledX)

def _pulse_falling_edge(x):
    templateStart=7.

    scaledX = x+templateStart
    scaledX = np.where(scaledX<templateStart,templateStart,scaledX)

    return _pulse_narrow(scaledX)

def _pulse_plateau_width(FB_WIDTH):
    return (FB_WIDTH-15.)*59.5/(124.-15.)

def _pulse_rising_edge_width(FB_WIDTH):
    return np.log(FB_WIDTH-12.)*1.91+5.

def _the_pulse(x, FB_WIDTH):
    scalar_input = np.isscalar(x)
    useX = np.atleast_1d(x)

    if FB_WIDTH <= 15:
        result = _pulse_narrow(useX * (15.0 / FB_WIDTH))
    else:
        plateau_width = _pulse_plateau_width(FB_WIDTH)
        rising_edge_width = _pulse_rising_edge_width(FB_WIDTH)

        result = np.where(
            useX <= rising_edge_width,
            _pulse_rising_edge(useX, rising_edge_width),
            np.where(
                useX <= rising_edge_width + plateau_width,
                np.ones_like(useX),
                _pulse_falling_edge(
                    useX - rising_edge_width - plateau_width
                ),
            ),
        )

    if scalar_input:
        return float(result[0])

    return result

class PocamKapu405WavelengthSource(WavelengthSource):
    """
    Wavelength sampler for IceTray/CLSim
    POCAMKapu405nmIsotropic.

    Uses the tabulated spectrum
    POCAM_XRL-400-5E_datasheet.txt.
    """

    class WavelengthParams(Structure):
        _fields_ = [
            ("_table", c_int64),
            ("_contrib", c_float),
        ]

    def __init__(self, numSamples: int = 4096) -> None:
        super().__init__(
            nRNGSamples=1,
            params={
                "WavelengthParams":
                PocamKapu405WavelengthSource.WavelengthParams
            },
        )

        wavelengths_nm = np.array([
            380.00, 381.25, 382.50, 383.75, 385.00,
            386.25, 387.50, 388.75, 390.00, 391.25,
            392.50, 393.75, 395.00, 396.25, 397.50,
            398.75, 400.00, 401.25, 402.50, 403.75,
            405.00, 406.25, 407.50, 408.75, 410.00,
            411.25, 412.50, 413.75, 415.00, 416.25,
            417.50, 418.75, 420.00, 421.25, 422.50,
            423.75, 425.00, 426.25, 427.50, 428.75,
            430.00,
        ])

        intensity = np.array([
            0.0147, 0.0153, 0.0184, 0.0249, 0.0364,
            0.0514, 0.0669, 0.0953, 0.1234, 0.1622,
            0.2335, 0.3073, 0.4108, 0.5157, 0.6145,
            0.7925, 0.8858, 0.9912, 0.9971, 0.9724,
            0.8651, 0.7674, 0.6179, 0.5160, 0.4123,
            0.3119, 0.2687, 0.2134, 0.1881, 0.1496,
            0.1303, 0.1035, 0.0830, 0.0705, 0.0575,
            0.0474, 0.0407, 0.0348, 0.0296, 0.0249,
            0.0208,
        ])

        # Same normalization used by CLSim
        intensity = intensity / 15.33201835267

        # Convert wavelength grid into Theia internal units
        wavelength = wavelengths_nm * u.nm

        # Integral of each linear segment
        dx = np.diff(wavelength)
        segment_area = (
            0.5
            * (intensity[:-1] + intensity[1:])
            * dx
        )

        cumulative = np.concatenate(
            ([0.0], np.cumsum(segment_area))
        )

        total_area = cumulative[-1]

        # Uniform CDF values for inverse-CDF lookup table
        probabilities = np.linspace(
            0.0,
            1.0,
            numSamples,
        )

        target_area = probabilities * total_area

        # Determine in which tabulated segment each CDF value lies
        segment_idx = np.searchsorted(
            cumulative,
            target_area,
            side="right",
        ) - 1

        segment_idx = np.clip(
            segment_idx,
            0,
            len(wavelength) - 2,
        )

        x0 = wavelength[segment_idx]
        y0 = intensity[segment_idx]

        x1 = wavelength[segment_idx + 1]
        y1 = intensity[segment_idx + 1]

        local_area = (
            target_area
            - cumulative[segment_idx]
        )

        local_dx = x1 - x0
        slope = (y1 - y0) / local_dx

        # Invert the integral of the linear PDF:
        #
        # local_area = y0*t + 0.5*slope*t^2
        #
        # The rationalized form below is numerically stable.
        disc = np.sqrt(
            y0**2
            + 2.0 * slope * local_area
        )

        t = np.where(
            np.abs(slope) < 1e-15,
            local_area / y0,
            2.0 * local_area / (y0 + disc),
        )

        inverse_cdf = x0 + t

        # Ensure exact endpoints
        inverse_cdf[0] = wavelength[0]
        inverse_cdf[-1] = wavelength[-1]

        # Upload inverse CDF as lookup table
        table = Table(inverse_cdf)
        self._table_gpu = table.upload()

        self.setParams(
            _table=self._table_gpu.address,
            _contrib=total_area,
        )

    @property
    def sourceCode(self) -> str:
        return loadShader(
            "wavelengthsource/function.glsl"
        )


class PocamLightSource(LightSource):
    """
    Isotropic point light source radiating light in all directions uniformly.

    Parameters
    ----------
    wavelengthSource: WavelengthSource
        Source to sample wavelengths from. Required for forward mode.
    mediumIdx: int
        Index of the medium the light source is located in
    position: (float, float, float), default=(0.0, 0.0, 0.0)
        Position of the light source.
    timeRange: (float, float), default=(0.0, 100.0)
        start and stop time of the light source
    budget: float, default=1.0
        Total amount of energy or photons the light source distributes among the
        sampled photons per wavelength. Ignored if source emits particles.
    emitParticles: bool, default=False
        Whether to emit particles instead of rays. If so, `budget` will be
        ignored.
    """

    class LightParams(Structure):
        _fields_ = [
            ("position", vec3),
            ("mediumIdx", c_uint32),
            ("_timeProfileTable", buffer_reference), # das noch ersetzen!!
            ("_contribFwd", c_float),
            ("_contribBwd", c_float),
        ]

    def __init__(
        self,
        wavelengthSource: WavelengthSource | None = None,
        *,
        mediumIdx: int,
        position: tuple[float, float, float] = (0.0, 0.0, 0.0),
        budget: float = 1.0,
        emitParticles: bool = False,
    ) -> None:
        lam = wavelengthSource
        super().__init__(
            nRNGForward=0 if lam is None else lam.nRNGSamples + 4,
            nRNGBackward=2,
            wavelengthSource=wavelengthSource,
            params={"LightParams": PocamLightSource.LightParams},
            extra={"budget"}, 
        )
        if emitParticles and budget != 1.0:
            warn("Cannot set budget: Light source emits particles")
        # save params
        self._emitParticles = emitParticles
        self._budget = budget
        self._timeProfileTable = self._createTimeProfile()
        self.setParams(
            position=position,
            mediumIdx=mediumIdx,
            _timeProfileTable=self._timeProfileTable,
        )
    def _createTimeProfile(self):

        width = 2.0 # FB_WIDTH = 2 corresponds to IceTray PulseWidth = 1 ns
    
        class Dist:                                       #das ist die Wahrscheinlichkeitsdichte pdf
            def pdf(self, x):
                return _the_pulse(x, width)
    
        inv_cdf = NumericalInversePolynomial(             #Umwandlung in eine inverse cdf 
            Dist(),
            domain=(0.0, 120.0),                          #Tobias fragen: was macht diese Zeile?
        )
    
        u = np.linspace(0.0, 1.0, 1024, endpoint=False)   # Array mit 1024 äquidistanten Werten zwischen 0 und 1
        ppf = inv_cdf.ppf(u)   
        #Funktion ppf(u): (u zw. 0 und 1) -> Zeit
        
        ppf_table = Table(ppf)                             # als Lookup-Tabelle auf GPU kopieren(vereinfacht Arbeit)
        self._ppf_gpu = ppf_table.upload()
        return self._ppf_gpu.address
        
    @property
    def budget(self) -> float:
        """
        Total energy or amount of photons the light source distributes among all
        sampled rays per wavelength. Ignored if source emits particles.
        """
        return self._budget

    @budget.setter
    def budget(self, value: float) -> None:
        if self.emitsParticles:
            warn("Cannot set budget: Light source emits particles")
        self._budget = value

    @property
    def emitsParticles(self) -> bool:
        """Whether this light source emits particles instead of light rays"""
        return self._emitParticles

    @property
    def forwardSourceCode(self) -> str | None:
        if self.wavelengthSource is None:
            return None
        else:
            preamble = createPreamble(LIGHT_SOURCE_EMIT_PARTICLE=self.emitsParticles)
            code = loadShader("lightsource/pocam/forward.glsl") # hier habe ich was geändert
            return preamble + code

    @property
    def backwardSourceCode(self) -> str | None:
        if self.emitsParticles:
            return None
        else:
            return loadShader("lightsource/pocam/backward.glsl")

    def _finishParams(self, i: int) -> None:
        super()._finishParams(i)
        c = self.budget
        self.setParam("_contribFwd", c)
        # in forward parameter volume cancels with the probability
        # This is not the case in backward mode
        # -> divide contrib by parameter space volume
        c /= 4.0 * np.pi
        self.setParam("_contribBwd", c)
