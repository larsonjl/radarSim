#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 09:04:27 2017

@author: jakelarson
"""

import numpy as np

import scipy.stats as st

class my_pdf(st.rv_continuous):
    def _pdf(self,x):
        return 3*x**2  # Normalized over its range, in this case [0,1]

my_cv = my_pdf(a=0, b=1, name='my_pdf')


'''
def getRandomExp(N_samples, mu_exp):
    outQ = np.random.exponential(mu_exp, N_samples)
    outI = np.random.exponential(mu_exp, N_samples)
    return outI, outQ


def getRandomChi(N_samples, df):
'''