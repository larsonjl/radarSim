#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 10:52:59 2017

@author: jakelarson
"""
import numpy as np
import scipy.io as sio
import math
import scipy.interpolate as interpolate

# Auxilary functions 
def degrees2radians(degrees):
    return (np.pi/180.) * degrees


def db2linear(dB):
    return (10 ** (dB/10))


def linear2db(linear):
    return 10 * np.log10(linear)


def freq2lambda(frequency):
    return (3e8) / frequency

def expNoise(samples, meanValue):
    return np.random.exponential(scale=meanValue, size=samples)

def chiNoise(samples, meanValue, freedomDegrees):
    xx = np.linspace(0, 4 * meanValue, 1000)
    pdf = ((freedomDegrees * xx) / (meanValue**2)) * np.exp((-2 * xx) / meanValue)
    pdf = pdf / np.sum(pdf)
    cumVals = np.cumsum(pdf)
    invCdf = interpolate.interp1d(cumVals, xx)
    r = np.random.rand(samples)
    return invCdf(r)

# Input gain patterns
gainDat = sio.loadmat('./InputData/antenna_pattern_N02_elements.mat')
gainPattern = gainDat['antenna_gain_pattern_linear']
radPattern = gainDat['normalized_radiation_pattern']
phirange = gainDat['phi_range']


class simulationAttributes:
    runtime = 1  # s

class radarVals:

    class fftPts:
        def __init__(self):
            self.v = 128
            self.nameString = "Points used for FFT"
            self.units = "nPts"
 
    class transPower:
        def __init__(self):
            self.v = 10e3
            self.nameString = "Transmitted Power"
            self.units = "W"

    class antGain:
        def __init__(self, targetLocation):
            angleToRadar = math.asin(targetLocation[1] / 
                                     (np.sqrt(targetLocation[1]**2 + 
                                              targetLocation[0]**2)))
            indx = np.argmin(np.abs(phirange - angleToRadar))
            self.v = linear2db(gainPattern[indx])
            self.nameString = "Antennae Gain Pattern"
            self.units = "dB"

    class frequency:
        def __init__(self):
            self.v = 2.835e9
            self.nameString = "Transmit Frequency"
            self.units = "Hz"
    
    class wavenumber:
        def __init__(self):
            self.v = (2 * np.pi) / freq2lambda(radarVals.frequency().v)
            self.nameString = "wavenumber"
            self.units = "1/m"
    
    class wavelength:
        def __init__(self):
            self.v = freq2lambda(radarVals.frequency().v)
            self.nameString = "wavelength"
            self.units = "m"


    class tIpp:
        def __init__(self):
            self.v = 100e-6
            self.nameString = "Inter-pulse Period"
            self.units = "sec"

    class pulseWidth:
        def __init__(self):
            self.v = 400e-9
            self.nameString = "Pulse  Width"
            self.units = "sec"

    class totGates:
        def __init__(self):
            self.v = 240
            self.nameString = "Number of range  gates"
            self.units = ""

    class tFirstRangeGate:
        def __init__(self):
            self.v = 800e-9
            self.nameString = "Time  of first  A/D sample"
            self.units = "sec"

    class tDiff:
        def __init__(self):
            self.v = 400e-9
            self.nameString = "Time  between A/D samples"
            self.units = "sec"

    class rxNoise:
        def __init__(self):
            self.v = 3
            self.nameString = "Receiver Noise  Figure"
            self.units = "dB"

    class snrDetect:
        def __init__(self):
            self.v = 1
            self.nameString = "SNR detectibility threshold"
            self.units = "dB"

    class rxBandwidth:
        def __init__(self):
            self.v = 1 / radarVals.pulseWidth().v
            self.nameString = "Receiver Bandwidth"
            self.units = "Hz"

    class unambRange:
        def __init__(self):
            self.v = 0.5 * (1e-3) * 3e8 * radarVals.tIpp().v
            self.nameString = "Unambiguous Range"
            self.units = "km"

    class rangeResolution:
        def __init__(self):
            self.v = 0.5 * 3e8 * radarVals.tDiff().v
            self.nameString = "Range  Resolution"
            self.units = "m"

    class timeLastSample:
        def __init__(self):
            self.v = radarVals.tIpp().v + radarVals.pulseWidth().v
            self.nameString = "Time  of Last Sample"
            self.units = "s"

    class range2FirstSample:
        def __init__(self):
            self.v = radarVals.tFirstRangeGate().v * 3e8 / 2
            self.nameString = "Range  to first  sample"
            self.units = "m"

    class range2LastSample:
        def __init__(self):
            self.v = (1e-3) * (radarVals.timeLastSample().v * 3e8 / 2)
            self.nameString = "Range  to last sample"
            self.units = "km"

    # Location of gate window  boundaries
    class gateTimeIntervalArray:
        def __init__(self):
            tBeginning = radarVals.tFirstRangeGate().v - radarVals.tDiff().v/2
            tEnd = tBeginning + radarVals.tDiff().v * radarVals.totGates().v +\
                   radarVals.tDiff().v/2
            self.v = np.arange(tBeginning, tEnd, radarVals.tDiff().v)

    class gateTimeRangeArray:
        def __init__(self):
            tBeginning = (radarVals.tFirstRangeGate().v)
            tEnd = tBeginning + radarVals.tDiff().v * radarVals.totGates().v
            self.v = 0.5 * 3e8 * np.arange(tBeginning, tEnd,
                                           radarVals.tDiff().v)

    # Coherent integration
    class numIntegrate:
        def __init__(self):
            self.v = 1
            self.nameString = "Coherent integration number"
            self.units = ""


class targetVals:
    class radCrossSection:
        def __init__(self):
            self.v = 2
            self.nameString = "Radar  cross section"
            self.units = "m^2"

    class initPosition:
        def __init__(self):
            self.v = [0, 6e3]
            self.nameString = "Initial target position (x, y)"
            self.units = "m"

    class targetSpeed:
        def __init__(self):
            self.v =  (freq2lambda(radarVals.frequency().v) / radarVals.tIpp().v ) \
                        * (2/radarVals.fftPts().v)
            self.nameString = "Target speed"
            self.units = "m/s"

    class approachAngle:
        def __init__(self):
            self.v = 90
            self.nameString = "Target angle of attack"
            self.units = "degrees"


class initialValues:
    class powerLinear:
        def __init__(self):
            targetRange = np.linalg.norm(targetVals.initPosition().v)
            numer = radarVals.transPower().v * \
                    (db2linear(radarVals.antGain(
                            targetVals.initPosition().v).v)**2) * \
                    (freq2lambda(radarVals.frequency().v)**2) * \
                    targetVals.radCrossSection().v
            denom = ((4 * np.pi)**3) * (targetRange**4)
            self.v = (numer/denom)
            self.nameString = 'Pr'
            self.units = "Watts"

    class powerDB:
        def __init__(self):
            self.v = 10 * np.log10(initialValues.powerLinear().v)
            self.nameString = 'Pr'
            self.units = 'dBW'

    class powerDBw:
        def __init__(self):
            self.v = initialValues.powerDB().v + 30
            self.nameString = 'Pr'
            self.units = 'dBm'

    class SNR:
        def __init__(self):
            denom = (1.38e-23) * 290 * db2linear(radarVals.rxNoise().v) * \
                    radarVals.rxBandwidth().v
            linear_snr = (initialValues.powerLinear().v) / denom
            self.v = 10 * np.log10(linear_snr)
            self.nameString = 'SNR'
            self.units = 'dB'

    class time2target:
        def __init__(self):
            targetRange = np.linalg.norm(targetVals.initPosition().v)
            self.v = 2 * targetRange / 3e8
            self.nameString = "Time  to target"
            self.units = 's'

    class range2target:
        def __init__(self):
            self.v  = 1e-3 * np.linalg.norm(targetVals.initPosition().v)
            self.nameString = "Range  to target"
            self.units = 'km'
    
    
class  calculateValues:
    class findPower:
        def __init__(self,  targetLocation):
            targetRange = np.linalg.norm(targetLocation)
            numer  = radarVals.transPower().v * \
                    (db2linear(radarVals.antGain(targetLocation).v)**2) * \
                    (freq2lambda(radarVals.frequency().v)**2) * \
                    targetVals.radCrossSection().v
            denom = ((4 * np.pi)**3) * (targetRange**4)
            self.v = numer/denom
            self.units = 'Watts'

    class expectedTime:
        def __init__(self,  targetLocation):
            self.v  = 2 * (np.linalg.norm(targetLocation)/3e8)
            self.units = 's'
    
    class findPosition:
        def __init__(self,  previousLocation, attackAngle):
            targetGradient = [np.cos(degrees2radians(
                                        attackAngle)),
                              np.sin(degrees2radians(
                                      attackAngle))]
        
            x_target = previousLocation[0] + \
                        targetGradient[0] * targetVals.targetSpeed().v * \
                        radarVals.tIpp().v
            y_target = previousLocation[1] + \
                       targetGradient[1] * targetVals.targetSpeed().v * \
                       radarVals.tIpp().v
   
            self.v  =  [x_target, y_target]
        
    
class  finalValues:
    class  timeofPulse:
        def __init__(self,  finalPulse): 
            self.v  = finalPulse
            self.units = 's'
            self.nameString = "Time  of transmitted pulse" 
    
    
    class finalPowerLin:
        def __init__(self,  finalPower): 
            self.v  = finalPower
            self.units = 'Watts'
            self.nameString = "Pr"
            

    class finalPowerDB:
        def __init__(self,  finalPower):
            self.v  = linear2db(finalPower)
            self.units = 'dB'
            self.nameString = "Pr"
    
    
    class finalPowerDBm:
        def __init__(self,  finalPower):
            self.v  = linear2db(finalPower) + 30
            self.units = 'dBm'
            self.nameString = "Pr"
    
    
    class  finalSNR:
        def __init__(self,  finalPower):
            denom  = (1.38e-23) * 290 * db2linear(radarVals.rxNoise().v) *  \
                     radarVals.rxBandwidth().v
            linear_snr = (finalPower) / denom
            self.v  = 10 * np.log10(linear_snr)
            self.nameString = 'SNR'
            self.units = 'dB'
    
    
    class finalTime2Target:
        def __init__(self,  time2target):
            self.v  = time2target
            self.units = 's'
            self.nameString = "Time  to Target"
    
    
    class  finalRange2Target:
        def __init__(self,  dist2target):
            self.v  = (1e-3)  * dist2target
            self.units = 'km'
            self.nameString = "Range  to Target"
