#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:17:01 2022

@author: mikkonea
"""

from multiprocessing import Process, Manager
import numpy as np
from scipy.interpolate import interp1d

def f(d, l):
    d[1] = '1'
    x_out = np.max(d['x']) * np.random.random(10000)
    interpfun = interp1d(d['x'],d['data'])
    return interpfun(x_out)

if __name__ == '__main__':
    arr = np.random.random((100000))
    with Manager() as manager:
        d = manager.dict()
        d['data'] = arr
        d['x'] = np.linspace(0,100,arr.size)
        l = manager.list(range(10))

        p = Process(target=f, args=(d, l))
        p2 = Process(target=f, args=(d, l))
        p.start()
        p2.start()
        p.join()
        p2.join()

        print(d)
        print(l)
