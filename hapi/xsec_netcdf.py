#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:49:46 2019

@author: anttimi
"""


import netCDF4
import numpy as np
import os.path
import time

xsec_folder = './xsec/'
float_accuracy = 'f4' #this can also be set to 'f8'

def setup_file(basename,subst,wl):
    """
    This function sets up a netCDF4 datafile for cross-sections of substance
    subst for wavelengths wl. The file is called gas_subst.nc or aer_subst.nc.

    The spectra are static, which means that in your use case you need to
    store the reasonable maximum bandwidth.
    """

    if not os.path.isdir(xsec_folder):
        os.mkdir(xsec_folder)

    filename = xsec_folder + '%s_%s.nc' % (basename, subst)
    if os.path.isfile(filename):
        return
    print("Creating new file '%s'..." % filename)
    ds = netCDF4.Dataset(filename, 'w', format="NETCDF4")

    if basename == 'gas':
        subst_type = 'gas'
    elif basename == 'aer':
        subst_type = 'aerosol'
    else:
        subst_type = basename

    ds.description = "Cross-sections for %s '%s'." % (subst_type, subst)
    ds.history = "Created %s." % time.ctime(time.time())

    len_wl = wl.size

    #the dimensions are set to a variable for clarity of the nc data structure
    dim_tpps_idx = ds.createDimension('TPPs_index', 0)
    dim_wl = ds.createDimension('wavelength', len_wl)
    #dim_T = ds.createDimension('temperature', 1)
    #dim_p = ds.createDimension('pressure', 1)
    #dim_ps = ds.createDimension('partial_pressure', 1)
    #dim_wl = ds.createDimension('wavelength', 1)

    var_xsec = ds.createVariable('cross_section', float_accuracy, ('TPPs_index','wavelength',))
    var_T = ds.createVariable('temperature', float_accuracy, ('TPPs_index',))
    var_p = ds.createVariable('pressure', float_accuracy, ('TPPs_index',))
    var_ps = ds.createVariable('partial_pressure', float_accuracy, ('TPPs_index',))
    var_wl = ds.createVariable('wavelength', float_accuracy, ('wavelength',))

    var_wl[:] = wl

    var_xsec.units = '1/cm'
    var_T.units = 'K'
    var_p.units = 'atm'
    var_ps.units = 'atm'
    var_wl.units = 'nm'

    ds.close()

def input_xsec_data(basename,subst,T,P,Ps,wl,xsec):
    """
    Appends the data xsec(TPPs,wl) into the file basename_subst.nc.
    """
    filename = xsec_folder + '%s_%s.nc' % (basename, subst)
    if not os.path.isfile(filename):
        print("File '%s' not found! Halting..." % filename)
        return
    ds = netCDF4.Dataset(filename, 'a', format="NETCDF4")
    assert np.isclose(wl,ds['wavelength']).all(), ("The wavelengths 'wl' " +
    "aren't all equal with the ones in the file %s!" % filename)
    TPPs_len = T.size
    ds_len = ds['cross_section'].shape[0]

    start_idx = ds_len
    end_idx = start_idx + TPPs_len
    ds['temperature'][start_idx : end_idx] = T
    ds['pressure'][start_idx : end_idx] = P
    ds['partial_pressure'][start_idx : end_idx] = Ps
    ds['cross_section'][start_idx : end_idx, :] = xsec

    ds.close()

def load_xsec_data(basename,subst):
    """
    Loads all xsec-data from the file basename_subst.nc
    """
    filename = xsec_folder + '%s_%s.nc' % (basename, subst)
    if not os.path.isfile(filename):
        print("File '%s' not found! Halting..." % filename)
        return
    ds = netCDF4.Dataset(filename, 'r', format="NETCDF4")

    T = ds['temperature'][:].data
    P = ds['pressure'][:].data
    Ps = ds['partial_pressure'][:].data
    wl = ds['wavelength'][:].data
    xsec = ds['cross_section'][:].data

    ds.close()

    return (T,P,Ps,wl,xsec)

# -----------------------------------------
# The interpolation scheme for the database
# -----------------------------------------

def is_inside_box(x,a,b):
    """
    This function checks if x is within an n-box defined by a and b as its
    vertices.
    """
    ab = np.array([a,b])
    min_corner = ab.min(0)
    max_corner = ab.max(0)
    return np.all(x <= max_corner) and np.all(min_corner <= x)

def get_distances(x, p):
    n = p.shape[1]
    distance_weights = np.array([1,1,1])
    W = np.diag(distance_weights)
    dists = np.inf * np.ones(n)
    for i in range(0,n):
        dists[i] = (p[:,i] - x) @ W @ (p[:,i] - x)
    return dists

def find_closest_element(x, p):
    dists = get_distances(x, p)
    return np.argmin(dists)

def box_volume(a,b):
    measure_weights = np.array([1,1,1])
    return np.prod(np.multiply(np.absolute(a - b),measure_weights))

def find_smallest_box_bf(x, p):
    """
    This function finds the smallest enclosing box for the parameters by
    comparing all of them.

    Complexity: O(n^2)
    """
    n = p.shape[1]
    boxes = []
    measures = []
    for i in range(0,n):
        for j in range(0,i):
            if is_inside_box(x,p[:,i],p[:,j]):
                boxes.append((i,j))
                measures.append(box_volume(p[:,i],p[:,j]))
    if not boxes:
        return(-1,find_closest_element(x,p))
    else:
        #print("Bruteforce:")
        #print(np.min(measures))
        return(boxes[np.argmin(measures)])

def find_smallest_box(x, p, lim):
    """
    This function finds the smallest enclosing box for the parameters.

    Find the closest point and test other points if they enclose x. If they do,
    say that the box is smallest. If not, then select the second-closest point
    and try again. If a box isn't found, then the value is calculated.

    Complexity: O(n)
    """
    n = p.shape[1]
    if lim == -1:
        #in this case we interpolate from all the values in the database.
        lim = n
    dists = get_distances(x, p)
    dist_ord = dists.argsort()
    #print(dists[dist_ord[0]])
    #print(dists[dist_ord[1]])
    for l in range(min(lim,n)):
        for i in range(l + 1,n):
            if is_inside_box(x,p[:,dist_ord[l]],p[:,dist_ord[i]]):
                #print("heurestic:")
                #print(box_volume(p[:,dist_ord[0]],p[:,dist_ord[i]]))
                return (dist_ord[l],dist_ord[i])
    return(-1,dist_ord[0])

def interpolate(x,a,b,fa,fb):
    distance_weights = np.array([1,1,1])
    W = np.diag(distance_weights)
    da = a @ W @ a
    db = b @ W @ b
    t = da / (da + db)
    return fa * t + fb * (1 - t)

import xsec_netcdf

did_we_interp = []
interp_box_indices = []
interp_box_size = []

def generate_xsec(hitran_folder,wn,gas,T,P,Ps,cl_interp_p):

    global did_we_interp
    global interp_box_indices
    global interp_box_size
    print('Generating Voigt Shapes for: ' + gas)

    Tvals,Pvals,Psvals,wl,xsec = xsec_netcdf.load_xsec_data('gas',gas)

    TPPs = np.array([Tvals.ravel(), Pvals.ravel(), Psvals.ravel()])
    cs = np.zeros((T.size,wn.size))

    comp_xsec_idxs = []
    for alt_idx in range(0,alt.size):
        x = np.array([T[alt_idx],P_atm[alt_idx],Ps[alt_idx]])
        if TPPs.size > 0:
            for tpps_idx in range(TPPs.shape[1]):
                if np.all(np.isclose(x,TPPs[:,tpps_idx])):
                    #an almost exact same data point is found
                    #let's use that one straight away
                    i = tpps_idx
                    j = tpps_idx
                    break
            else:
                (i,j) = find_smallest_box(x,TPPs,cl_interp_p)
        else: #if there's no data yet
            i = -1

        if i == -1:
            #no box found -> compute new point
            comp_xsec_idxs.append(alt_idx)
        else:
            #box found -> interpolate
            cs[alt_idx,:] = interpolate(x,TPPs[:,i],TPPs[:,j],
            xsec[i][:],xsec[j][:])
        if i == -1:
            did_we_interp.append(False)
            interp_box_indices.append([-1,-1])
            interp_box_size.append(-1)
        else:
            did_we_interp.append(True)
            interp_box_indices.append([i,j])
            interp_box_size.append(box_volume(TPPs[:,i],TPPs[:,j]))

    print("Calculating %d / %d cross-sections..." %
          (len(comp_xsec_idxs), len(T)))

    cs[comp_xsec_idxs,:] = calc_voigt_shapes_allisotopes(hitran_folder,
          wn,gas,T[comp_xsec_idxs],P[comp_xsec_idxs],Ps[comp_xsec_idxs])

    if len(comp_xsec_idxs) > 0:
        xsec_netcdf.input_xsec_data('gas',gas,T[comp_xsec_idxs],
                P[comp_xsec_idxs],Ps[comp_xsec_idxs],wn,cs[comp_xsec_idxs,:])

    return cs
