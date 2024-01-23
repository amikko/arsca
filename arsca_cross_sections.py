#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:49:46 2019

@author: anttimi
"""

# TODO: Hide the internal HAPI number for gases, because they're not very 
# useful to use. 

import netCDF4
import numpy as np
import os.path
import time
import sys
path_to_folder_containing_this_file = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, path_to_folder_containing_this_file + '/hapi/')
import hapi

import arsca_utilities
case = 'foo'
configuration = 'bar'

xsec_db_folder = './csdb/xsec/'
hapi_db_folder = './csdb/hitran/'

float_accuracy = 'f4' #this can also be set to 'f8'

cl_interp_p = 16

global hapi_db_started, hapi_db_current_folder
hapi_db_started = False
hapi_db_current_folder = ""

def case_folder():
    return xsec_db_folder + case + '/'

def case_folder_hitran():
    return hapi_db_folder + case + '/'

def set_up_datafolders():
    if not os.path.isdir('./csdb'):
        os.mkdir('./csdb')
    for folder in [xsec_db_folder, hapi_db_folder]:
        if not os.path.isdir(folder):
            os.mkdir(folder)

def fetch_lines(wn_range,gas,gas_id):
    
    if os.path.isfile(case_folder_hitran() + '%s.data' % gas):
        print("Line data already downloaded for '%s'. Skipping download..." % gas)
    else:
        hapi.fetch(gas,
                    gas_id,
                    1, #the most common isotopologue is used
                    wn_range[0],
                    wn_range[1])

def calculate_xsec(gas,wn_range,wn_step,p,T,dil,func_selection='Lorentz'):
    function_listing = {'Lorentz' : hapi.absorptionCoefficient_Lorentz,
                        'Doppler' : hapi.absorptionCoefficient_Doppler,
                        'Voigt'   : hapi.absorptionCoefficient_Voigt,
                        'HT'      : hapi.absorptionCoefficient_HT,
                        'SDVoigt' : hapi.absorptionCoefficient_SDVoigt}
    line_broadening_func = function_listing[func_selection]
    nu, coef = line_broadening_func(
                        SourceTables = gas,
                        Diluent = {'air' : dil, 'self' : 1.0 - dil},
                        WavenumberStep = wn_step,
                        WavenumberRange = wn_range,
                        Environment = {'p' : p, 'T' : T})
    return nu, coef

# --------------------------------
# The cross-section database stuff
# --------------------------------

def setup_file(subst,wn):
    """
    This function sets up a netCDF4 datafile for cross-sections of substance
    subst for wavenumbers wn.

    The spectra are static, which means that in your use case you need to
    store the reasonable maximum bandwidth.
    """

    arsca_utilities.if_folder_not_exist_then_create(case_folder())

    filename = case_folder() + '%s.nc' % (subst)
    if os.path.isfile(filename):
        return
    print("Creating new file '%s'..." % filename)
    ds = netCDF4.Dataset(filename, 'w', format="NETCDF4")

    ds.description = "Cross-sections for substance '%s'." % (subst)
    ds.history = "Created %s." % time.ctime(time.time())

    len_wn = wn.size

    #the dimensions are set to a variable for clarity of the nc data structure
    dim_tpps_idx = ds.createDimension('TPPs_index', 0)
    dim_wl = ds.createDimension('wavenumber', len_wn)

    var_xsec = ds.createVariable('cross_section', float_accuracy, ('TPPs_index','wavenumber',))
    var_T = ds.createVariable('temperature', float_accuracy, ('TPPs_index',))
    var_p = ds.createVariable('pressure', float_accuracy, ('TPPs_index',))
    var_ps = ds.createVariable('partial_pressure', float_accuracy, ('TPPs_index',))
    var_wn = ds.createVariable('wavenumber', float_accuracy, ('wavenumber',))

    var_wn[:] = wn

    var_xsec.units = 'cm^2'
    var_T.units = 'K'
    var_p.units = 'atm'
    var_ps.units = 'atm'
    var_wn.units = '1/cm'

    ds.close()

def input_xsec_data(subst,T,P,Ps,wn,xsec):
    """
    Appends the data xsec(TPPs,wn) into the file basename_subst.nc.
    """
    filename = case_folder() + '%s.nc' % (subst)
    if not os.path.isfile(filename):
        print("File '%s' not found! Setting up..." % filename)
        setup_file(subst,wn)
    ds = netCDF4.Dataset(filename, 'a', format="NETCDF4")
    assert np.isclose(wn,ds['wavenumber']).all(), ("The wavenumbers 'wn' " +
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

def load_xsec_data(subst):
    """
    Loads all xsec-data from the file basename_subst.nc
    """
    filename = case_folder() + '%s.nc' % (subst)
    if not os.path.isfile(filename):
        print("File '%s' not found! Computing all the cross-sections..." % filename)
        return (None, None, None, None, None)
    ds = netCDF4.Dataset(filename, 'r', format="NETCDF4")

    T = ds['temperature'][:].data
    P = ds['pressure'][:].data
    Ps = ds['partial_pressure'][:].data
    wn = ds['wavenumber'][:].data
    xsec = ds['cross_section'][:].data

    ds.close()

    return (T,P,Ps,wn,xsec)

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

def generate_xsec(wn_range,wn_step,gas,gas_id,T,P,Ps,func_selection='Lorentz'):
    global hapi_db_started, hapi_db_current_folder
    if not hapi_db_started or case not in hapi_db_current_folder:
        #Either the HAPI db hasn't been stated or
        #we've changed the case but had not reimported the module.
        set_up_datafolders()
        hapi_db_current_folder = hapi_db_folder + case
        hapi.db_begin(hapi_db_current_folder)
        hapi_db_started = True
        
    fetch_lines(wn_range,gas,gas_id)
    
    did_we_interp = []
    interp_box_indices = []
    interp_box_size = []

    print('Generating line shapes for: ' + gas)

    Tvals,Pvals,Psvals,wn,xsec = load_xsec_data(gas)
    fresh_csdb = (type(Tvals) == type(None)) #couldn't find the nc-file

    comp_xsec_idxs = []
    if fresh_csdb:
        for alt_idx in range(0,T.size):
            comp_xsec_idxs.append(alt_idx)
        wnAmt = np.array(np.arange(wn_range[0],wn_range[1],wn_step)).size
        cs = np.zeros((T.size,wnAmt))
        #return 0, cs
    else:
        TPPs = np.array([Tvals.ravel(), Pvals.ravel(), Psvals.ravel()])
        cs = np.zeros((T.size,wn.size))
        for alt_idx in range(0,T.size):
            x = np.array([T[alt_idx],P[alt_idx],Ps[alt_idx]])
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

    #we asssume that the only interactions are between generic air and gas itself
    dil = 1.0 - Ps/P

    for comp_idx in comp_xsec_idxs:
        wn, cs[comp_idx,:] = calculate_xsec(gas,wn_range,wn_step,
    P[comp_idx],T[comp_idx],dil[comp_idx],func_selection)

    #store the computed data to the database
    if len(comp_xsec_idxs) > 0:
        input_xsec_data(gas,T[comp_xsec_idxs],
                P[comp_xsec_idxs],Ps[comp_xsec_idxs],wn,cs[comp_xsec_idxs,:])

    return wn, cs
