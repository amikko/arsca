#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:23:03 2021

@author: mikkonea
"""

class Medium:
    self.positions = [] # this will be np.array
    self.abs = []
    self.sca = []
    self.emi = []
    
    def __init__(self,absorbers,scatterers,emitters):
        self.abs = absorbers
        self.sca = scatterers
        self.emi = emitters
        
class Absorber:
    