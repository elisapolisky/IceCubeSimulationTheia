#!/usr/bin/env python3
from optparse import OptionParser
from os.path import expandvars
#
usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_flashes.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12344,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",
                  default=expandvars("$I3_DATA/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-p", "--position", dest="POSITION",
                  help="Position of the flasher in the detector 'x,y,z'")
parser.add_option("-f", "--flashertype", type="string", default="POCAMKapu405nmIsotropic", dest="FLASHER_TYPE",
                  help = "The type of flasher to use. See /src/simclasses/public/simclasses/I3FlasherPulse.h for a list of options.")
parser.add_option("","--nphotons", type="int", default=1_000_000_000, dest="NPHOTONS",
                  help = "The number of photons to generate per event.")
parser.add_option("","--pulsewidth", type="int", default=1, dest="PULSEWIDTH",
                  help = "The pulse width of a single flash.")
parser.add_option("", "--weighted", action="store_true", default=False, dest="WEIGHTED",
                  help="Use weighted photons.")
parser.add_option("", "--unweighted", action="store_false", dest="WEIGHTED",
                    help="Use unweighted photons.")
parser.add_option("", "--holeiceparam", type="string", 
                  default=expandvars("/data/user/leidensc/software/icetray/pocam/build/ice-models/resources/models/ANGSENS/angsens_unified_ftp/as.0"), 
                  dest="HOLEICEPARAM",
                  help="Full path to the Hole Ice Parameterization file (z.B. .../as.0 oder .../as.1)")

# ... (Rest des OptionParsers)
#parser.add_option("", "--no-gpu", action="store_false", dest="USEGPU",
#                  help="Do not use GPUs.")
# parse cmd line args, bail out if anything is not understood

(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

from icecube.icetray import I3Tray, I3Units
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services, simclasses

tray = I3Tray()

flasher_pos = dataclasses.I3Position(*(float(v) for v in options.POSITION.split(",")))
print("Source position: ", flasher_pos)
flasher_type = simclasses.I3FlasherPulse.FlasherPulseType.names[options.FLASHER_TYPE]
print("Flasher type: ", flasher_type)
print("Number of flash events: ", options.NUMEVENTS)
print("Number of photons per event: ", options.NPHOTONS)
print("Weighted photons: ", options.WEIGHTED)
# a random number generator
try:
    randomService = phys_services.I3SPRNGRandomService(
        seed = options.SEED,
        nstreams = 10000,
        streamnum = options.RUNNUMBER)
except AttributeError:
    randomService = phys_services.I3GSLRandomService(
        seed = options.SEED*10000 + options.RUNNUMBER,
    )

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2012,
               DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

from icecube.clsim import POCAMFlashInfoIsotropic
tray.AddModule(POCAMFlashInfoIsotropic, "POCAMFlasher", 
               FlasherPulseType=flasher_type,
               PhotonPulseSeriesName="I3FlasherPulseSeries",
               NumberOfPhotons=options.NPHOTONS,
               FlasherPosition=flasher_pos,
               PulseWidth= options.PULSEWIDTH*I3Units.ns, 
               NEvents=options.NUMEVENTS
               )
print("Uniform AS")
# start photon propagation with CLSim
tray.AddSegment(clsim.I3CLSimMakePhotons, "goCLSIM",
    UseGPUs=True,
    UseCPUs=False,
    UseGeant4=False,
    UseI3PropagatorService=False,
    RandomService=randomService,
    DoNotParallelize=False,
    UnweightedPhotons=(not options.WEIGHTED),
    StopDetectedPhotons=True,
    #OMKeyMaskName=I3Vector,
    PhotonSeriesName='I3PhotonSeriesMap',
    MCPESeriesName='',
    FlasherPulseSeriesName="I3FlasherPulseSeries",
    GCDFile=options.GCDFILE,
    IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_ftp-v3"),
    HoleIceParameterization=options.HOLEICEPARAM,
    DOMOversizeFactor=1.,
    DOMEfficiency=1.5,
    IgnoreSubdetectors=['IceTop', 'NotOpticalSensor'],
   )
def EmptyI3MCTree(frame):
    if not frame.Has('I3MCTree'):
        frame.Put('I3MCTree', dataclasses.I3MCTree())
    return
tray.AddModule(EmptyI3MCTree, "EmptyI3MCTree", streams=[icetray.I3Frame.DAQ])
tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo])

tray.Execute()
