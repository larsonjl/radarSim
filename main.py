# -*- coding: utf-8 -*-
"""
This is the main method for an object oriented radar simulator
All variables and calculation methods are given as objects in the
radarAttributes.py file.
Auxiliary functions to create statistics text file and plots are given
in helperFunctions.py
@author: jakelarson
"""
import numpy as np
import helperFunctions
from modelAttributes import radarVals, targetVals, simulationAttributes, \
                            initialValues, calculateValues, linear2db, \
                            db2linear, finalValues

# Create text file with initial conditions
# helperFunctions.initattributes2text(radarVals, targetVals, initialValues)
# Initialize first target location
targetLocation = targetVals.initPosition().v
gateArr = np.zeros((int(simulationAttributes.runtime/radarVals.tIpp().v),
                        radarVals.totGates().v))

# Initialize time, index
time = 0
i = 0
nPts = 128
ptNum = 0

# Initialize arrays to store simulator data
# transLength = int(2 + simulationAttributes.runtime / radarVals.tIpp().v)
transLength = nPts * radarVals.numIntegrate().v
transmitArr = np.zeros((transLength))
sampleArr = np.zeros((radarVals.totGates().v, transLength))
snrArr = np.zeros((radarVals.totGates().v, transLength))
arrivalTimeArr = np.zeros((transLength))
targetRangeArr = np.zeros((transLength))
returnPowerArr = np.zeros((transLength))
trueRadialVelArr = np.zeros((transLength))
radialVelFromSNR = np.zeros((transLength))
iOut = np.zeros((nPts))
qOut = np.zeros((nPts))
integrateTime = np.zeros((nPts))

cohereTrack = 0
noiseOn = True

# Main method to perform simulation
# while time <= simulationAttributes.runtime + radarVals.tIpp().v:
while ptNum < nPts:
    cohereTrack = 0
    iIntegrate = np.zeros((radarVals.numIntegrate().v))
    qIntegrate = np.zeros((radarVals.numIntegrate().v))
    integrateTimeArr = np.zeros((radarVals.numIntegrate().v))
    while cohereTrack < radarVals.numIntegrate().v:
        # Update Target position, return time, and return power
        targetLocation = calculateValues.findPosition(targetLocation).v
        returnTime = calculateValues.expectedTime(targetLocation).v
        returnPower = calculateValues.findPower(targetLocation).v
        targetRange = np.linalg.norm(targetLocation)
        
        if noiseOn == True:
            returnPower += chiNoise(1, 2e-14, 4)

        # Store transmit, arrival times, and power
        transmitArr[i] = time
        arrivalTimeArr[i] = returnTime
        targetRangeArr[i] = targetRange
        returnPowerArr[i] = linear2db(returnPower) + 30  # dBm
    
        # Store coherent integration information
        iIntegrate[cohereTrack] = np.sqrt(returnPower) * np.cos(
                                    - 2 * targetRange * (radarVals.wavenumber().v))
        qIntegrate[cohereTrack] = np.sqrt(returnPower) * np.sin(
                                    - 2 * targetRange * (radarVals.wavenumber().v))  
        integrateTimeArr[cohereTrack] = time + returnTime

        if noiseOn == True:
            iIntegrate[cohereTrack] += expNoise(1, 2e-14)
            qIntegrate[cohereTrack] += expNoise(1, 2e-14)

        # Timing of gates relative to when pulse is sent (length =  #gates + 1)
        gateTimeArray = radarVals.gateTimeIntervalArray().v

        # Find gate indx of return pulse,
        gateLocation = (gateTimeArray >= returnTime) * \
                       (gateTimeArray < returnTime + radarVals.pulseWidth().v)
        gateIndx = np.where(gateLocation > 0)

        # Store return power in appropriate bins
        sampleArr[gateIndx[0] - 1, i] = returnPower
        snrArr[gateIndx[0] - 1, i] = finalValues.finalSNR(returnPower).v

        # Increase time by interpulse period
        time += radarVals.tIpp().v
        i += 1
        cohereTrack += 1

    iOut[ptNum] = np.mean(iIntegrate)
    qOut[ptNum] = np.mean(qIntegrate)
    integrateTime[ptNum] = np.mean(integrateTimeArr)

    ptNum += 1

#### ===== Make Plots for pt3 q2 ===== ####
print("Initial Position %f, %f m" %(targetVals.initPosition().v[0], targetVals.initPosition().v[1] ))
print("Final Position %f, %f m" %(targetLocation[0], targetLocation[1] ))
print("Time on Target: %f s" %time)
# helperFunctions.generatePt3q2Figures(integrateTime, iOut, qOut, returnPowerArr, 'approach')
# helperFunctions.generatePt3q3Figures(integrateTime, iOut, qOut, 1/radarVals.tIpp().v, radarVals.wavelength().v, 'approach')
helperFunctions.generatePt3q4Figures(integrateTime, iOut, qOut, 1/radarVals.tIpp().v, radarVals.wavelength().v, 'approach')

'''
#### ===== Print final values after simulation ===== ####
helperFunctions.printFinalValues(transmitArr[-1],
                                 db2linear(returnPowerArr[-1] - 30),
                                 arrivalTimeArr[-1], targetRangeArr[-1])

#### ===== Generate final figures ===== ####
helperFunctions.generateFigures(transmitArr, arrivalTimeArr,
                                returnPowerArr, sampleArr, snrArr)

#### =====  Calculate radial velocity and plot ===== ####

targRangeFromSNR = np.zeros((transLength))

for i in range(0, transLength):
    index = np.argmax(snrArr[:,i])
    targRangeFromSNR[i] = radarVals.gateTimeIntervalArray().v[index] * 3e8/2

timeArr = np.arange(0,  simulationAttributes.runtime + 2*radarVals.tIpp().v,
                    radarVals.tIpp().v)
windowSize = 20000
timeDiff = timeArr[windowSize] - timeArr[0]
timeDiff2 = timeArr[10] - timeArr[0]

drdtArr = np.zeros((transLength))
drdtTrue = np.zeros((transLength))

for i in range(10, transLength):
    drdtTrue[i] = (targetRangeArr[i] - targetRangeArr[i-10]) / timeDiff2

for i in range(windowSize, transLength):
    drdtArr[i] = (targRangeFromSNR[i] - targRangeFromSNR[i-windowSize]) / \
           timeDiff

#### ==== Plot radial velocity ===== ####
helperFunctions.radialVelocityPlot(timeArr, drdtArr, drdtTrue)
'''
