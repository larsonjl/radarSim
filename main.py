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
                            db2linear, finalValues, chiNoise, expNoise

# Create text file with initial conditions
# helperFunctions.initattributes2text(radarVals, targetVals, initialValues)
# Initialize first target location
def runSimulation(numIntegrate, attackAngle, noiseOn):
    # Initialize time, index
    time = 0
    i = 0
    nPts = 128
    ptNum = 0
    
    targetLocation = targetVals.initPosition().v
    gateArr = np.zeros((int(simulationAttributes.runtime/radarVals.tIpp().v),
                        radarVals.totGates().v))
    # Initialize arrays to store simulator data
    # transLength = int(2 + simulationAttributes.runtime / radarVals.tIpp().v)
    transLength = nPts * numIntegrate
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
    
    # Main method to perform simulation
    # while time <= simulationAttributes.runtime + radarVals.tIpp().v:
    while ptNum < nPts:
        cohereTrack = 0
        iIntegrate = np.zeros((numIntegrate))
        qIntegrate = np.zeros((numIntegrate))
        integrateTimeArr = np.zeros((numIntegrate))
        while cohereTrack < numIntegrate:
            # Update Target position, return time, and return power
            targetLocation = calculateValues.findPosition(targetLocation, 
                                                          attackAngle).v
            returnTime = calculateValues.expectedTime(targetLocation).v
            returnPower = calculateValues.findPower(targetLocation).v
            targetRange = np.linalg.norm(targetLocation)
            print(returnPower)
            if noiseOn == True:
                returnPower = chiNoise(1, returnPower, 4)
            # Store transmit, arrival times, and power
            transmitArr[i] = time
            arrivalTimeArr[i] = returnTime
            targetRangeArr[i] = targetRange
            returnPowerArr[i] = linear2db(returnPower) + 30  # dBm
    
                
            # Store coherent integration information
            iIntegrate[cohereTrack] = np.sqrt(returnPower) * np.cos(
                                        - 2 * targetRange * (
                                                radarVals.wavenumber().v))
            qIntegrate[cohereTrack] = np.sqrt(returnPower) * np.sin(
                                         - 2 * targetRange * (
                                                 radarVals.wavenumber().v))  
            integrateTimeArr[cohereTrack] = time + returnTime
    
            if noiseOn == True:
                iIntegrate[cohereTrack] = iIntegrate[cohereTrack] + \
                                                expNoise(1, 2e-14)
                qIntegrate[cohereTrack] = qIntegrate[cohereTrack] + \
                                                expNoise(1, 2e-14)
    
            # Timing of gates relative to when pulse is sent (length =  #gates + 1)
            gateTimeArray = radarVals.gateTimeIntervalArray().v
    
            # Find gate indx of return pulse,
            gateLocation = (gateTimeArray >= returnTime) * \
                           (gateTimeArray < 
                            returnTime + radarVals.pulseWidth().v)
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

    print("Initial Position %f, %f m" %(targetVals.initPosition().v[0], 
                                        targetVals.initPosition().v[1] ))
    print("Final Position %f, %f m" %(targetLocation[0], targetLocation[1] ))
    print("Time on Target: %f s" %time)

    return integrateTime, iOut, qOut, returnPowerArr

#### ===== Make Plots for pt3 q2 ===== ####
'''
# Run question 1 #
print("Running Question 1")
helperFunctions.generatePt3q1Figures()

# Run Question 2 #
print("Running Question 2 approach and recede")
integrateTime, iOut, qOut, returnPowerArr = runSimulation(1, 270, False)
helperFunctions.generatePt3q2Figures(integrateTime, iOut, qOut, 
                                         returnPowerArr, 'approach')
integrateTime, iOut, qOut, returnPowerArr = runSimulation(1, 90, False)
helperFunctions.generatePt3q2Figures(integrateTime, iOut, qOut, 
                                     returnPowerArr, 'recede')

# Run Question 3 #
print("# ==== Running Question 3 approach and recede ==== #")
integrateTime, iOut, qOut, returnPowerArr = runSimulation(1, 270, False)
helperFunctions.generatePt3q3Figures(integrateTime, 
                                     iOut, qOut, 1/radarVals.tIpp().v, 
                                     radarVals.wavelength().v, [-25, 25],
                                     1,'approach')
integrateTime, iOut, qOut, returnPowerArr = runSimulation(1, 90, False)
helperFunctions.generatePt3q3Figures(integrateTime,
                                     iOut, qOut, 1/radarVals.tIpp().v, 
                                     radarVals.wavelength().v, [-25,25],
                                     1, 'recede')
integrateTime, iOut, qOut, returnPowerArr = runSimulation(5, 90, False)
helperFunctions.generatePt3q3Figures(integrateTime, 
                                     iOut, qOut, 1/radarVals.tIpp().v,
                                     radarVals.wavelength().v, [-25, 25], 5, 
                                     'recede_ncoh5')

'''
# Run Question 4 #
print("# ==== Running Question 4 approach and recede ==== #")
integrateTime, iOut, qOut, returnPowerArr = runSimulation(1, 90, True)
helperFunctions.generatePt3q4Figures(integrateTime,
                                     iOut, qOut, 1/radarVals.tIpp().v, 
                                     radarVals.wavelength().v, 'recede_100less')