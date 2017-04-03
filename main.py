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

# Initialize arrays to store simulator data
transLength = int(2 + simulationAttributes.runtime / radarVals.tIpp().v)
transmitArr = np.zeros((transLength))
sampleArr = np.zeros((radarVals.totGates().v,
                      transLength))
snrArr = np.zeros((radarVals.totGates().v,
                      transLength))
arrivalTimeArr = np.zeros((transLength))
targetRangeArr = np.zeros((transLength))
returnPowerArr = np.zeros((transLength))
trueRadialVelArr   = np.zeros((transLength))
radialVelFromSNR   = np.zeros((transLength))

# Main method to perform simulation
while time <= simulationAttributes.runtime + radarVals.tIpp().v:
    
    # Update Target position, return time, and return power
    targetLocation = calculateValues.findPosition(targetLocation).v
    returnTime = calculateValues.expectedTime(targetLocation).v
    returnPower = calculateValues.findPower(targetLocation).v
    # Store transmit, arrival times, and power
    transmitArr[i] = time
    arrivalTimeArr[i] = returnTime
    targetRangeArr[i] = np.linalg.norm(targetLocation)
    returnPowerArr[i] = linear2db(returnPower) + 30  # dBm

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



