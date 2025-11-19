#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 15:00:56 2024

@author: mikkonea
"""

import numpy as np

def plancks_law(wl,T):
    # black-body emissivity (https://en.wikipedia.org/wiki/Planck%27s_law)
    # W m^-2 sr^-1 nm^-1
    c = 299792458 # m s^-1
    nu = c / (wl/1e9) # Hz
    h = 6.62607015e-34 # J Hz^-1
    k_B = 1.380649e-23 # J / K^-1
    coeff = 2 * h * nu ** 3 / (c ** 2)
    nu_rep = np.repeat(nu.reshape((wl.size,1)),T.size,axis=1)
    denom = np.exp(h * nu_rep / (k_B * T)) - 1
    B_nu = coeff / denom.T
    return B_nu