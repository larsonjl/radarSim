#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:23:12 2017

@author: jakelarson
"""
import inspect
import matplotlib.pyplot as plt
import seaborn as sns
import numpy  as np
from modelAttributes import radarVals, linear2db, finalValues


# Print  out initial attributes
def initattributes2text(radarVals, targetVals, initialValues):
    with open("outputData.txt", "w") as myfile:
        myfile.writelines("## ====  Initial Radar  Attributes ====  ##\n")
        for name, obj in inspect.getmembers(radarVals):
            if name[0:2] != '   ' and name != 'gateTimeIntervalArray' and \
                                      name != 'gateTimeRangeArray':

                myfile.writelines(obj().nameString + ' = ' + str(obj().v) + 
                                  ' ' + obj().units + "\n")

        myfile.writelines("\n")
        myfile.writelines("## ====  Initial Target  Attributes ====  ##\n")
        for name,  obj in inspect.getmembers(targetVals):
            if name[0:2] != '   ':
                myfile.writelines(obj().nameString + ' = ' + str(obj().v) +
                                  ' ' + obj().units + "\n")
        myfile.writelines("\n")
        myfile.writelines("## ====  Initial Values  ====  ##\n")
        for name,  obj in inspect.getmembers(initialValues):
            if name[0:2] != '   ':
                myfile.writelines(obj().nameString + ' = ' + str(obj().v) +
                                  ' ' + obj().units + "\n")


# Print  out final  values
def printFinalValues(finalTime, finalPower, finalTime2Target, finalRange):
    with open("outputData.txt", "a") as myfile:
        myfile.writelines("\n")
        myfile.writelines("## ====  Final  Values  ====  ##\n")
        for name,  obj in inspect.getmembers(finalValues):
            if name[0:2] != '   ':
                if name == 'timeofPulse':
                    myfile.writelines(obj(finalTime).nameString + ' = ' +
                                      str(obj(finalTime).v) +
                                      ' ' + obj(finalTime).units + "\n")
            if name == 'finalPowerLin':
                myfile.writelines(obj(finalPower).nameString + ' = ' +
                                  str(obj(finalPower).v) +
                                  ' ' + obj(finalPower).units + "\n")

            if name == 'finalPowerDB':
                myfile.writelines(obj(finalPower).nameString + ' = ' + 
                                  str(obj(finalPower).v) +
                                  ' ' + obj(finalPower).units + "\n")

            if name == 'finalPowerDBm':
                myfile.writelines(obj(finalPower).nameString + ' = ' +
                                  str(obj(finalPower).v) +
                                  ' ' + obj(finalPower).units + "\n")

            if name == 'finalTime2Target':
                myfile.writelines(obj(finalTime2Target).nameString+' = ' +
                                  str(obj(finalTime2Target).v) +
                                  ' ' + obj(finalTime2Target).units + "\n")

            if name == 'finalRange2Target':
                myfile.writelines(obj(finalRange).nameString + ' = ' +
                                  str(obj(finalRange).v) +
                                  ' ' + obj(finalRange).units + "\n")


# Plot  data
def generateFigures(transmitArr, arrivalTimeArr, returnPowerArr, sampleArr, snrArr):
    sns.set_context('talk')
    sns.set_style('ticks')

    # Generate time  series  plots
    plt.figure()
    ax = plt.subplot(211)
    plt.plot(transmitArr, 1e6 * arrivalTimeArr)
    plt.xlabel("Time of transmitted pulse  (s)")
    plt.ylabel("Arrival time  (microseconds)")
    plt.grid()
    ax.spines['top'].set_color('none')
    plt.xlim((0, 24))
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')

    ax = plt.subplot(212)
    plt.plot(transmitArr, returnPowerArr)
    plt.xlabel("Time of transmitted pulse  (s)")
    plt.ylabel("Received power  (dBm)")
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.xlim((0, 24))
    plt.tight_layout()
    plt.savefig('./Figures/timeSeries.png', dpi=350)
    plt.close()

    plt.figure()
    ax = plt.subplot(111)
    gateRanges = 1e-3 * (radarVals.gateTimeRangeArray().v)
    gateRanges = gateRanges - 0.5 * (gateRanges[1] - gateRanges[0])  # shift  to
    # sampleArr = np.ma.masked_where(sampleArr == 0, sampleArr)

    snrArr = np.ma.masked_where(snrArr < radarVals.snrDetect().v, snrArr)

    xx, yy = np.meshgrid(transmitArr, gateRanges)
    plt.pcolormesh(xx, yy, linear2db(snrArr) + 30, cmap=plt.cm.viridis)
    plt.ylim((0, 4))
    plt.xlim((0, 24))
    cb = plt.colorbar(label='dBm')
    cb.outline.set_visible(False)
    plt.xlabel('Time of transmitted pulse  (s)')
    plt.ylabel('Range to target  (km)')
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.savefig('./Figures/snr.png', dpi=350)
    plt.close()


def radialVelocityPlot(timeArr, drdtArr, drdtTrue):
    plt.figure()
    ax = plt.subplot(111)
    plt.plot(timeArr, drdtArr, label='Calculated')
    plt.plot(timeArr, drdtTrue, label = 'True')
    plt.legend()
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.xlim(0.1,)
    plt.xlabel('Time (s)')
    plt.ylabel('Target radial velocity (m/s)')
    plt.savefig('./Figures/radialVelocity.png', dpi=350)
    plt.close()
