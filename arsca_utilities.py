#arsca-utilities
#some various functions, which
#   * don't really fit anywhere else
#   * don't have dependencies

import os
import uuid
import numpy as np

import netCDF4
from scipy.interpolate import interp1d

def unique_filename(folder,name,suffix,overwrite=False):
    filename_candidate = name + suffix
    file_exists = os.path.isfile(folder + filename_candidate)
    if file_exists and not overwrite:
        unique_part = str(uuid.uuid4()).split('-')[0]
        uniq_fname = name + '-' + unique_part + suffix
        print("File '%s' already exists in '%s'! Using a new filename '%s'." %
        (filename_candidate, folder, uniq_fname))
        file_name = uniq_fname
    else:
        file_name = filename_candidate
    file_path = folder + file_name
    return file_path

def if_folder_not_exist_then_create(folder_path):
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

def cartesian_product(*arrays):
    #The rightmost array is the "innermost", so that cartprod[i] and
    #cartprod[i+1] stand for [x_j,y_k,z_l] and [x_j,y_k,z_l+1] respectively

    #This is from https://stackoverflow.com/questions/11144513/cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def linear_to_multi_idx(idx,shape):
    #shape is a tuple with the matrix's dimensions
    mult_idx = [np.NaN] * len(shape)
    for dim_idx in range(len(shape)):
        dim_prod = int(np.prod(shape[(dim_idx+1):]))
        if idx // dim_prod > 0:
            mult_idx[dim_idx] = idx // dim_prod
        else:
            mult_idx[dim_idx] = 0
        idx = idx % dim_prod
    return tuple(mult_idx)

def _integrate_phase_function(costheta,ph):
    #this is the exact integral of the linearly interpolated function
    dct = np.diff(costheta)
    avgs = (0.5 * ph[1:] + ph[:-1])
    return np.dot(dct, avgs)

def aerosol_input_from_online_mie_calc(infilename,outfilename,old_normalization):
    # This creates the mie_IR.dat-like formatted phase matrix look-up tables
    # from a https://omlc.org/calc/mie_calc.html -styled input file
    miedata_in = np.genfromtxt(infilename)
    # Miedata_in has angles from -180 to 180 and we want from 0 to 180
    start = np.where(miedata_in[:,0] == 0)[0][0]
    miedata = miedata_in[start:,:] # Leave out everything before the 0 angle
    miedata = miedata[::-1,:] # flip it around, so it is from -180 to 0
    costheta = np.cos(np.pi / 180.0 * miedata[:,0]) # now the angles are from -1 to 1
    if old_normalization:
        # This is the original implementation
        I = (_integrate_phase_function(costheta,miedata[:,3])
            +_integrate_phase_function(costheta,miedata[:,2]))
        coeff = 2 / I
        M = np.zeros((costheta.size,17))
        M[:,0] = costheta
        M[:,1] = coeff * miedata[:,3]   # I_parallel to M_11
        M[:,6] = coeff * miedata[:,2]   # I_perp to M_22
        M[:,11] = miedata[:,8]  # S_33 to M_33
        M[:,12] = miedata[:,9]  # S_34 to M_34
        M[:,15] = -miedata[:,9] # -S_34 to M_43
        M[:,16] = miedata[:,8]  # S_33 to M_44, which it actually is
    else:
        # This is the new implementation to harmonize the matrices across 
        # different sources.
        I = 0.5 * np.trapz(miedata[:,5] + miedata[:,6],costheta)
        M = np.zeros((costheta.size,17))
        M[:,0] = costheta
        M[:,1] = miedata[:,6] / I  # I_parallel to M_11
        M[:,6] = miedata[:,5] / I  # I_perp to M_22
        M[:,11] = miedata[:,8] / I  # S_33 to M_33
        M[:,12] = miedata[:,9] / I  # S_34 to M_34
        M[:,15] = -miedata[:,9] / I # -S_34 to M_43
        M[:,16] = miedata[:,8] / I  # S_33 to M_44, which it actually is
    
    delim = '\t'
    np.savetxt(outfilename,M,delimiter=delim)

def convert_IQUV_muller_file_to_IIUV(infilename,outfilename):
    #When you have a usual muller matrix file, but it is in wrong basis
    
    IIUV2IQUV = np.array([[1.0, 1.0, 0.0, 0.0],
                          [1.0,-1.0, 0.0, 0.0],
                          [0.0, 0.0, 1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])

    IQUV2IIUV = np.linalg.inv(IIUV2IQUV)
    
    indata = np.genfromtxt(infilename)
    outdata = np.zeros_like(indata)
    outdata[:,0] = indata[:,0]
    for i in range(indata.shape[0]):
        outdata[i,1:] = (IQUV2IIUV @ indata[i,1:].reshape((4,4)) @ IQUV2IIUV).ravel()
        
    np.savetxt(outfilename,outdata,fmt='%.9f')