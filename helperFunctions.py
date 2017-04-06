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
from scipy import signal
from modelAttributes import radarVals, linear2db, finalValues, expNoise,
                            chiNoise


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


def generatePt3q1Figures():
    sns.set_context('talk')
    sns.set_style('ticks')
    meanValue = 2e-14
    sampleNum = 10e3
    gSamples = expNoise(sampleNum, meanValue)
    
    # Create theoretical exponential
    xx = np.linspace(0, np.max(gSamples), 500)
    yy = (1 / meanValue) * np.exp(-xx / meanValue)

    # Plot
    plt.figure()
    ax = plt.subplot(111)
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.hist(gSamples, 50, normed=True)
    plt.plot(xx, yy)
    plt.ylabel('P(x)')
    plt.xlabel('x')
    plt.xlim(0,)
    plt.title('Exponential distribution theoretical and sampled (mean=2e-14, n=10000)')
    plt.tight_layout()
    plt.savefig('./Figures/hw3q1_fig1.png', dpi=350)
    plt.close()
    
    # Create chi-square
    chiSamples = chiNoise(10e3, 2e-14, 4)
    yyChi = ((4 * xx) / (2e-14**2)) * np.exp((-2 * xx) / 2e-14)
    plt.figure()
    ax = plt.subplot(111)
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.hist(chiSamples, 50, normed=True)
    plt.plot(xx, yyChi)
    plt.ylabel('P(x)')
    plt.xlabel('x')
    plt.xlim(0,1.2e-13)
    plt.title('Chi-squared theoretical and sampled (mean=2e-14, n=10e3, degrees=4)')
    plt.tight_layout()
    plt.savefig('./Figures/hw3q1_fig2.png', dpi=350)
    plt.close()


# Plotting Stuff
def generatePt3q2Figures(arrivalTimeArr, iOut, qOut, powerMeas, saveInfo):
    sns.set_context('talk')
    sns.set_style('ticks')
    # Q, I vs time 
    plt.figure(figsize=(10,5))
    ax = plt.subplot(111)
    plt.plot(arrivalTimeArr * 1e3, iOut, label ='i')
    plt.plot(arrivalTimeArr * 1e3, qOut, label = 'q')
    plt.xlabel("Arrival time (mS)")
    plt.ylabel("Amplitude")
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.tight_layout()
    plt.legend()
    plt.savefig('./Figures/hw3q2_fig1_%s.png'%saveInfo, dpi=350)
    plt.close()
    
    # Phasor diagram
    plt.figure(figsize=(5,5))
    sns.set_context('talk')
    sns.set_style('ticks')
    ax = plt.subplot(111)
    plt.plot(iOut[0:25], qOut[0:25])
    plt.scatter(iOut[0], qOut[0], c='g')
    plt.text(iOut[0], qOut[0], 'Start')
    
    plt.scatter(iOut[24], qOut[24], c='r')
    plt.text(iOut[24], qOut[24], 'End')

    plt.xlabel('I')
    plt.ylabel('Q')
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.tight_layout()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.set_aspect('equal', 'datalim')
    plt.title('Phasor diagram (bursts 1-25)')
    plt.savefig('./Figures/hw3q2_fig2_%s.png'%saveInfo, dpi=350)
    plt.close()
    
    
    # Verify return power = i^2 + q^2
    powerIQ = (iOut**2 + qOut**2)
    powerIQ = linear2db(powerIQ) + 30
    plt.figure(figsize=(10,10))
    ax = plt.subplot(111)
    plt.plot(powerIQ, powerMeas)
    plt.xlabel('I^2 + Q^2 (dbM)')
    plt.ylabel('Return Power (dbM)')
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.tight_layout()
    plt.title('I^2 + Q^2 versus returned power')
    plt.savefig('./Figures/hw3q2_fig3_%s.png'%saveInfo, dpi=350)
    plt.close()

def makePeriodogram(saveName, velocity, Pxx_den, transmitArr, iOut, 
                    qOut, xlims, plotQ, twinx):
    Pxx_den = (2/radarVals.wavelength().v) * Pxx_den

    # Generate figure
    sns.set_context('talk')
    sns.set_style('ticks')
    plt.figure()
    ax = plt.subplot(212)
    plt.scatter(velocity, Pxx_den)
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.xlim(xlims[0], xlims[1])
    print(np.max(Pxx_den))
    plt.ylim(np.min(Pxx_den), 1.1 * np.max(Pxx_den))
    plt.xlabel('Radial Velocity (m/s)')
    plt.ylabel('PSD [Watts/(m/s)]')
    ax = plt.subplot(211)
    plt.plot(transmitArr, iOut, label='i')
    if plotQ == True:
        plt.plot(transmitArr, qOut, label='q')
    plt.legend()
    plt.grid()
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.xlim(0,)
    plt.tight_layout()
    plt.savefig('./Figures/%s' % saveName, dpi=350)
    plt.close()
        
def makePeriodogramQ4(saveName, velocity, Pxx_den, transmitArr, iOut, 
                      qOut, xlims, plotQ, twinx):
        # TODO: Check this    
        Pxx_den = (2/radarVals.wavelength().v) * Pxx_den
    
        # Generate figure
        sns.set_context('talk')
        sns.set_style('ticks')
        plt.figure(figsize=(8,14))
        ax = plt.subplot(413)
        Pxx_den_db = np.zeros((len(Pxx_den)))
        for i in range(0,len(Pxx_den)):
            Pxx_den_db[i] = linear2db(Pxx_den[i])
        plt.scatter(velocity, Pxx_den_db, s=20)
        plt.grid()
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['right'].set_color('none')
        ax.yaxis.set_ticks_position('left')
        plt.xlim(xlims[0], xlims[1])
        plt.ylim(-210, 20 + np.max(Pxx_den_db))
        plt.xlabel('Radial Velocity (m/s)')
        plt.ylabel('PSD [db Watts/(m/s)]')
    

        ax = plt.subplot(412)
        plt.scatter(velocity, Pxx_den, s=20)
        plt.grid()
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['right'].set_color('none')
        ax.yaxis.set_ticks_position('left')
        plt.xlim(xlims[0], xlims[1])
        plt.ylim(np.min(Pxx_den), 1.1 * np.max(Pxx_den))
        plt.xlabel('Radial Velocity (m/s)')
        plt.ylabel('PSD [Watts/(m/s)]')

        ax = plt.subplot(411)
        plt.plot(transmitArr, iOut, label='i')
        if plotQ == True:
            plt.plot(transmitArr, qOut, label='q')
        plt.legend()
        plt.grid()
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['right'].set_color('none')
        ax.yaxis.set_ticks_position('left')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.xlim(0,)
        plt.tight_layout()
        
        ax = plt.subplot(412)
        # Calculate SNR
        maxSignal = np.max(Pxx_den)
        meanSignal = np.mean((Pxx_den[(Pxx_den < maxSignal) * (Pxx_den > 0)]))
        snr = linear2db(maxSignal/meanSignal)
        
        # Calculate target velocity
        velocityTarg = velocity[np.argmax(Pxx_den)]
        plt.figtext(0.2, 0.11, 'Velocity Estimate = %.2f m/s'%velocityTarg)
        plt.figtext(0.2, 0.08, 'Noise Power = %.2f dB'%linear2db(meanSignal))
        plt.figtext(0.2, 0.05, 'SNR = %.2f dB'%snr)
        plt.savefig('./Figures/%s' % saveName, dpi=350)
        plt.close()
        


def generatePt3q3Figures(transmitArr, iOut, qOut, samplingFreq, 
                         wavelength, xlims, numIntegrate, saveInfo):

    # Calculate periodogram
    realPlusI = iOut[:] + 1j * qOut[:]
    f1, Pxx_realPlusI = signal.periodogram(realPlusI,
                                           samplingFreq / numIntegrate)
    f2, Pxx_I = signal.periodogram(iOut, samplingFreq / numIntegrate )
    velocity1 = (f1 * wavelength) / 2  # radial velocity
    velocity2 = (f2 * wavelength) / 2  # radial velocity    
    makePeriodogram('hw3q3_fig1_%s'%saveInfo, velocity1,
                    Pxx_realPlusI, transmitArr, iOut, qOut, xlims,
                    plotQ= True, twinx = False)
    makePeriodogram('hw3q3_fig2_%s'%saveInfo, velocity2,
                    Pxx_I, transmitArr, iOut, qOut,
                    xlims, plotQ = False, twinx = False)


def generatePt3q4Figures(transmitArr, iOut, qOut, samplingFreq, 
                         wavelength, saveInfo):
    
    realPlusI = iOut[:] + 1j * qOut[:]
    f1, Pxx_realPlusI = signal.periodogram(realPlusI, samplingFreq)
    f2, Pxx_I = signal.periodogram(iOut, samplingFreq)
    velocity1 = (f1 * wavelength) / 2  # radial velocity
    velocity2 = (f2 * wavelength) / 2  # radial velocity    
    makePeriodogramQ4('hw3q4_fig1_%s'%saveInfo, velocity1,
                    Pxx_realPlusI, transmitArr, iOut, qOut,
                    [-50, 50], plotQ= True, twinx = True)
    makePeriodogramQ4('hw3q4_fig2_%s'%saveInfo, velocity2,
                    Pxx_I, transmitArr, iOut, qOut, [-50, 50],
                    plotQ = False, twinx = True)

    
    

def generateFigures(transmitArr, arrivalTimeArr, returnPowerArr, 
                    sampleArr, snrArr):

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
