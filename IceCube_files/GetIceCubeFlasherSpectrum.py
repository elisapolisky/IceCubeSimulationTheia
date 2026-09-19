#
# Copyright (c) 2011, 2012 Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# Copyright (c) 2011, 2012 the IceCube Collaboration <http://www.icecube.wisc.edu>
# SPDX-License-Identifier: ISC
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#
# $Id: GetIceCubeFlasherSpectrum.py 108199 2013-07-12 21:33:08Z nwhitehorn $
#
# @file GetIceCubeFlasherSpectrum.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from icecube.simclasses import I3CLSimFunctionFromTable
from icecube.simclasses import I3CLSimFunctionDeltaPeak
from icecube.simclasses import I3FlasherPulse

from icecube.icetray import I3Units

import numpy, math
from os.path import expandvars

def GetIceCubeFlasherSpectrumData(spectrumType):
    if spectrumType == I3FlasherPulse.FlasherPulseType.LED340nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_20mA_pulseCurrent.txt"), unpack=True)
        #data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_200mA_pulseCurrent.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 24.306508         # (20mA pulse)  pre-calculated normalization constant (not really necessary for the generation algorithm)
        #data[1] /= 22.323254        # (200ma pulse) pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.LED370nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_370nm_emission_spectrum_cw_measured.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 15.7001863        # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.LED405nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_405nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 8541585.10324     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.LED450nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_450nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 21.9792812618     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.LED505nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_505nm_emission_spectrum_cw_measured.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 38.1881           # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMKapu405nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_XRL-400-5E_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 15.33201835267    # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMKapu465nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_NSPB300B_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 24.5122477        # pre-calculated normalization constant (not really necessary for the generation algorithm)  
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMLMG365nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_XSL-365-5E_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 15.642481874785    # pre-calculated normalization constant (not really necessary for the generation algorithm)   
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMLMG405nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_RLT405500MG_analytic.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 1.63987039        # pre-calculated normalization constant (not really necessary for the generation algorithm)  
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMLMG450nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_PL-TB450B_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 1.547889456       # pre-calculated normalization constant (not really necessary for the generation algorithm)  
    elif spectrumType == I3FlasherPulse.FlasherPulseType.POCAMLMG520nmIsotropic:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/POCAM_PLT5_520_B1-3_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 1.0670525         # pre-calculated normalization constant (not really necessary for the generation algorithm)  
    else:
        raise RuntimeError("invalid spectrumType")

    return data

def GetIceCubeFlasherSpectrum(spectrumType = I3FlasherPulse.FlasherPulseType.LED405nm):
    if spectrumType in [I3FlasherPulse.FlasherPulseType.SC1, I3FlasherPulse.FlasherPulseType.SC2]:
        # special handling for Standard Candles:
        # these currently have single-wavelength spectra
        return I3CLSimFunctionDeltaPeak(337.*I3Units.nanometer)
    else:
        data = GetIceCubeFlasherSpectrumData(spectrumType)
        return I3CLSimFunctionFromTable(data[0], data[1])
