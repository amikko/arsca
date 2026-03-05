#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:35:11 2020

@author: mikkonea
"""

import sys
import numpy as np
import netCDF4
import arsca

step_siro = 1.0
#noph_siro = 100

try:
    casesel = sys.argv[1]
    noph_siro = int(sys.argv[2])
except:
    print("RUNNING DEFAULT SETTINGS")
    casesel = 'e1'
    noph_siro = 100

vzas = [0, 9, 18, 26, 34, 41, 48, 54, 60, 65, 70, 74, 78, 81, 84, 86, 88, 89, 90]
vaas = [0.001, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
szas = [30, 60, 80, 87, 90, 93, 96, 99]
#szas = [30, 60]
saas = [0]
solar_lon_e6 = [50,90,110,130]
salts= [0 + 1e-4, 120 - 1e-4] # sensor altitude, kilometres
#salts=[0 + 1e-4]
atmos_thickness = 120
R_earth = 6371

siro_custom_settings = {
    'noph' : noph_siro,
    'ratm' : R_earth + atmos_thickness,
    'step' : step_siro}

#def rotate_IQUV(vaa_angs,stokes):
#    for idx_vang,vang in enumerate(vaa_angs):
def old_rotate_IQUV(vza_angs,stokes):
    for idx_vang,vang in enumerate(vza_angs):
        ang = np.pi / 180.0 * (vang)
        ang = np.pi / 180.0 * 0
        rotmat = np.array([[1.0, 0.0,             0.0,             0.0],
                           [0.0, np.cos(2 * ang), np.sin(2 * ang), 0.0],
                           [0.0, np.sin(2 * ang), np.cos(2 * ang), 0.0],
                           [0.0, 0.0,             0.0,             1.0]])
        stoks = stokes[...,idx_vang,:,:]
        for i in range(stoks.shape[0]):
            for j in range(stoks.shape[1]):
                #print(stoks.shape)
                stokvek = stoks[i,j,:]
                stokvek[1] = stokvek[1]*(-1)
                stokvek[2] = stokvek[2]*(-1)
                stoks[i,j,:] = rotmat @ stokvek
        stokes[...,idx_vang,:,:] = stoks[...]
    return stokes

def rotate_IQUV(vza_angs,stokes):
    stokes[...,1] = -stokes[...,1]
    stokes[...,2] = -stokes[...,2]
    return stokes

def homogenous_layers(alttau):
    # nasty little layerses...
    # internally Siro linearly interpolates between levels, so more than one homogenous
    # layer is a pain to set up, but let's do it
    smaller_alts = alttau.copy()
    smaller_alts[:,0] = smaller_alts[:,0] - 1e-4 # 10 cm lower
    new_alttau = np.concatenate((alttau,smaller_alts))
    alttau = new_alttau[new_alttau[:,0].argsort()]
    alttau[:,1] = np.roll(alttau[:,1],1)
    alttau[-1,1] = alttau[-2,1]
    alttau = alttau[1:,:]
    return alttau

if casesel == 'd1':
    # single layer Rayleigh
    settings = {'n_sca' : 1,
                'n_abs' : 1, # needs to be one, but filled with zeroes
                'tau_sca' : 0.5,
                'wavelength' : 450,
                'albedo' : 0.0,
                'n_lay' : 1,                
                'rayleigh_depol' : 0.03} # currently hardcoded into Siro
elif casesel == 'd2':
    # single layer Rayleigh + albedo
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'tau_sca' : 0.1,
                'wavelength' : 450,
                'albedo' : 0.3,
                'n_lay' : 1,                
                'rayleigh_depol' : 0.03}

elif casesel == 'd3':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/ws_sphere_500nm.dat'"
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'tau_sca' : 0.2,
                'ssa' : 0.975683,
                'wavelength' : 500,
                'albedo' : 0.0,
                'n_lay' : 1}
elif casesel == 'd4':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/spheroid_350nm.dat'"
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'tau_sca' : 0.2,
                'ssa' : 0.787581,
                'wavelength' : 350,
                'albedo' : 0.0,
                'n_lay' : 1}
elif casesel == 'd5':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/water10mic_800nm.dat'"
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'tau_sca' : 5.0,
                'ssa' : 0.999979,
                'wavelength' : 800,
                'albedo' : 0.0,
                'n_lay' : 1}
elif casesel == 'd6':
    # single layer Rayleigh + albedo
    siro_custom_settings['BRF_FILENAME'] = "'input/brdf/mischenko_ocean_2ms.nc4'"
    siro_custom_settings['brdf_reflection'] = True
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'tau_sca' : 0.1,
                'wavelength' : 477,
                'albedo' : 0.3,
                'n_lay' : 1,
                'rayleigh_depol' : 0.03}
elif casesel == 'e1':
    alttau = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
    settings = {'n_sca' : 1,
                'n_abs' : 1, 
                'alttau' : homogenous_layers(alttau),
                'wavelength' : 450,
                'albedo' : 0.0}    
elif casesel == 'e2':
    alttausca = np.genfromtxt('datafiles/atmos/tau_rayleigh_320nm_usstd.dat',comments='#')
    alttauabs = np.genfromtxt('datafiles/atmos/tau_absorption_320nm_usstd.dat',comments='#')
   
    settings = {'n_sca' : 1,
                'n_abs' : 1, 
                'alttausca' : homogenous_layers(alttausca),
                'alttauabs' : homogenous_layers(alttauabs),
                'wavelength' : 320,
                'albedo' : 0.0}    
elif casesel == 'e3':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/desdust_450nm.dat'"
    alttausca = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
    alttauabs = np.genfromtxt('datafiles/atmos/tau_absorption_450nm_usstd.dat',comments='#')
    settings = {'n_sca' : 2,
                'n_abs' : 2, 
                'alttausca' : homogenous_layers(alttausca),
                'alttauabs' : homogenous_layers(alttauabs),
                'ssa' : 0.83447903,
                'tau_sca' : 0.5,
                'aer_lims' : [0.0, 3.0],
                'wavelength' : 450,
                'albedo' : 0.0}
elif casesel == 'e4':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/desdust_450nm.dat'"
    siro_custom_settings['AER2_FILENAME'] = "'input/miefiles/sulfate_450nm.dat'"
    alttausca = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
    alttauabs = np.genfromtxt('datafiles/atmos/tau_absorption_450nm_usstd.dat',comments='#')
    settings = {'n_sca' : 3,
                'n_abs' : 3,
                'alttausca' : homogenous_layers(alttausca),
                'alttauabs' : homogenous_layers(alttauabs),
                'ssa' : [0.83447903, 0.99999994],
                'tau_sca' : [0.5, 0.05],
                'aer_lims' : [[0.0, 3.0],[20.0, 21.0]],
                'wavelength' : 450,
                'albedo' : 0.0}
elif casesel == 'e5':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/cirrus_450nm.dat'"
    alttausca = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
    alttauabs = np.genfromtxt('datafiles/atmos/tau_absorption_450nm_usstd.dat',comments='#')
    settings = {'n_sca' : 2,
                'n_abs' : 2, 
                'alttausca' : homogenous_layers(alttausca),
                'alttauabs' : homogenous_layers(alttauabs),
                'ssa' : 1.000000,
                'tau_sca' : 1.0,
                'aer_lims' : [10.0, 11.0],
                'wavelength' : 450,
                'albedo' : 0.0}
elif casesel == 'e6':
    # single layer Rayleigh + albedo
    siro_custom_settings['BRF_FILENAME'] = "'input/brdf/mischenko_ocean_5ms.nc4'"
    siro_custom_settings['brdf_reflection'] = True
    alttausca = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
    alttauabs = np.genfromtxt('datafiles/atmos/tau_absorption_450nm_usstd.dat',comments='#')
    settings = {'n_sca' : 1,
                'n_abs' : 1,
                'alttausca' : homogenous_layers(alttausca),
                'alttauabs' : homogenous_layers(alttauabs),
                'wavelength' : 450,
                'albedo' : 0.3,
                'n_lay' : 1,
                'rayleigh_depol' : 0.03}
    szas = solar_lon_e6
else:
    print('bad case selection: %s' % casesel)
    halt
    
outfilename = 'iprt_phase3_Siro.nc'

def save_results_to_nc(casesel,idx_sza,stokes):
    var_names = ['%s_sza','%s_saa','%s_vza','%s_vaa','%s_zout','radiance_%s']
    var_names = [vname % casesel for vname in var_names]
    if not casesel == 'e6':
        stokes_in = np.reshape(stokes,(len(salts),1,1,len(vzas),len(vaas),4))
    else:
        stokes_in = np.reshape(stokes,(1,1,1,61,61,4))
    with netCDF4.Dataset(outfilename,'a') as ds:
        if not casesel == 'e6':
            ds[var_names[-1]][:,idx_sza,0,:,:,:] = stokes_in
        else:
            ds[var_names[-1]][:,0,idx_sza,:,:,:] = stokes_in

def transform_case_results(casesel,idx_sza):
    var_names = ['%s_sza','%s_saa','%s_vza','%s_vaa','%s_zout','radiance_%s']
    var_names = [vname % casesel for vname in var_names]
    with netCDF4.Dataset(outfilename,'a') as ds:
        vaa_angs = ds[var_names[-3]]
        if not casesel == 'e6':
            stokes = ds[var_names[-1]][:,idx_sza,0,:,:,:]
            stokes = rotate_IQUV(vaa_angs, stokes)
            ds[var_names[-1]][:,idx_sza,0,:,:,:] = stokes
        else:
            stokes = ds[var_names[-1]][:,0,idx_sza,:,:,:]
            stokes = rotate_IQUV(vaa_angs, stokes)
            ds[var_names[-1]][:,0,idx_sza,:,:,:] = stokes
def tau_to_xsec(tau,thickness):
    # converts optical thicknesses in km to xsec in cm2
    # assuming the number density to be 1/cm3
    return tau / thickness / 1e5
    
IIUV2IQUV = np.array([[1.0, 1.0, 0.0, 0.0],
                      [1.0,-1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0]])
# Siro uses IIUV-basis internally for the polarized RT calculations
IQUV2IIUV = np.linalg.inv(IIUV2IQUV)
schemes = [[1,2,0,0,2,1,0,0,0,0,3,4,0,0,-4,3],
           [1,2,0,0,2,5,0,0,0,0,3,4,0,0,-4,6],
           [1,2,0,0,2,1,0,0,0,0,3,4,0,0,-4,3],
           [1,2,0,0,2,5,0,0,0,0,3,4,0,0,-4,6],
           [1,2,0,0,2,5,0,0,0,0,3,4,0,0,-4,6],
           [1,2,0,0,2,5,0,0,0,0,3,4,0,0,-4,6]]

def aerosol_phase_matrix_file_generation():
    infiles = ['./datafiles/mie_results/waso.mie.dat',
               './datafiles/mie_results/sizedistr_spheroid.dat',
               './datafiles/mie_results/watercloud.mie.dat',
               './datafiles/mie_results/desert.cdf', # ssa: 0.83447903
               './datafiles/mie_results/sulfate.cdf', # ssa: 0.99999994
               './datafiles/mie_results/ic.ghm.baum.cdf'] # ssa: 1.000000
    outfiles = ['./rt_solvers/siro/input/miefiles/ws_sphere_500nm.dat',
                './rt_solvers/siro/input/miefiles/spheroid_350nm.dat',
                './rt_solvers/siro/input/miefiles/water10mic_800nm.dat',
                './rt_solvers/siro/input/miefiles/desdust_450nm.dat',
                './rt_solvers/siro/input/miefiles/sulfate_450nm.dat',
                './rt_solvers/siro/input/miefiles/cirrus_450nm.dat']
    mults = [1.0,
             1.0,
             1.0,
             1.0,
             1.0,
             0.5]
    for i in range(len(infiles)):
        infile = infiles[i]
        if '.cdf' in infile:
            try:
                with netCDF4.Dataset(infile) as ds:
                    wl_idx = np.argmin(np.abs(ds['wavelen'][:] - 0.45))
                    reff_idx = 18 if i == 5 else 0
                    angs = ds['theta'][wl_idx,reff_idx,0,:]
                    arr = np.zeros((angs.size,7))
                    arr[:,0] = angs
                    arr[:,1:] = mults[i] * ds['phase'][wl_idx,0,:,:].T
            except:
                # the ic.ghm.baum.cdf file is not included to the git repo, because
                # it is so large. The generated cirrus_450nm is included however.
                pass
        else:
            arr = np.genfromtxt(infile,comments='#')
            arr = arr.T
        rows = arr.shape[0]
        muller_data = np.zeros((rows,17))
        for r in range(rows):
            muller_data[r,0] = np.cos(arr[r,0] * np.pi / 180.0)
            for j in range(16):
                if schemes[i][j] != 0:
                    muller_data[r,j+1] = arr[r,schemes[i][j]] if schemes[i][j] > 0 else -arr[r,-schemes[i][j]]
            muller_data[r,1:] = (IQUV2IIUV @ muller_data[r,1:].reshape((4,4)) @ IQUV2IIUV).ravel()
        if 'desdust' in outfiles[i]:
            n_angle_pad = 10
            new_muller = np.zeros((rows + n_angle_pad,17))
            new_angs = np.linspace(muller_data[-2,0],muller_data[-1,0],n_angle_pad+1)
            new_muller[:rows,:] = muller_data[:,:]
            new_muller[-1,:] = muller_data[-1,:]
            new_muller[(rows-1):,0] = new_angs
            dist = muller_data[-1,0] - muller_data[-2,0]
            for r in range(n_angle_pad+1):
                for j in range(16):
                    new_muller[rows + r-1,j+1] = (new_muller[rows + r-1,0] - muller_data[-2,0]) / dist * new_muller[-1,j+1] + (muller_data[-1,0] - new_muller[rows + r-1,0]) / dist * new_muller[rows-2,j+1]
            new_muller = np.delete(new_muller,rows-1,axis=0)
            muller_data = new_muller
        np.savetxt(outfiles[i],muller_data)
            
aerosol_phase_matrix_file_generation()

#for idx_sza in range(len(szas)):
#    transform_case_results(casesel, idx_sza)
#halt

for idx_sza in range(len(szas)):
    #for idx_sza in [1]:
    #for idx_sza in [0]:    
    compute_siro = True
    
    arsca.simu.create_siro_settings(siro_custom_settings)
    
    arsca.set_case("iprt_%s" % casesel)
    arsca.set_configuration("base_conf_%s" % casesel)
    
    n_wn = 1
    wl = np.array([settings['wavelength']])
    
    n_wl = n_wn
    if casesel in ['d1', 'd2', 'd3', 'd4', 'd5', 'd6']:
        n_lev = settings['n_lay'] + 1
        altitudes = np.linspace(0,atmos_thickness,n_lev)
    elif casesel == 'e1':
        altitudes = settings['alttau'][:,0]
    elif casesel in ['e2','e3','e4','e5','e6']:
        altitudes = settings['alttausca'][:,0]
    else:
        print(settings.keys())
    n_medium_positions = altitudes.size
    n_coordinate = 3
    
    # MEDIUM DEFINITIONS
    
    medium = {}
    
    #medium['position'] is now a n_medium_positions x 3
    altitudes_ = altitudes.reshape((n_medium_positions, 1)) + R_earth
    position = np.zeros((n_medium_positions, n_coordinate))
    position[:,[0]] = altitudes_
    medium['position'] = position
    
    # Siro can be buggy, if used with non-scattering medium.
    n_scatterer = settings['n_sca']
    n_emitter = 0 #emission will be added later
    n_absorber = settings['n_abs']
    
    #scatterers and their cross-sections
    medium['scatterer'] = np.zeros((n_medium_positions,n_scatterer))
    medium['scattering_cross_section'] = np.zeros((n_medium_positions,n_wl,n_scatterer))
    
    medium['absorber'] = np.zeros((n_medium_positions,n_absorber))
    medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber))
    if casesel == 'd1' or casesel == 'd2' or casesel == 'd6':
        # d1 : single layer Rayleigh
        # d2 : single layer Rayleigh + albedo
        # d6 : single layer Rayleigh + sea surface
        medium['scatterer'][:,0] = 1.0 
        medium['scattering_cross_section'] = tau_to_xsec(settings['tau_sca'],atmos_thickness)
        
        # zero stands for rayleigh-scattering
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel_parameter'] = np.ones((1,n_scatterer)) * settings['rayleigh_depol']
        
    elif casesel in ['d3', 'd4', 'd5']:
        medium['scatterer_kernel'] = np.ones((n_scatterer,)) # this is to use the aerosols
        medium['scatterer'][:,0] = 1.0 
        medium['scattering_cross_section'] = tau_to_xsec(settings['tau_sca'] * settings['ssa'],atmos_thickness)
        medium['scatterer_kernel_parameter'] = np.zeros((1,n_scatterer))
        medium['absorber'][:,0] = 1.0
        medium['absorbing_cross_section'] = tau_to_xsec(settings['tau_sca'] * (1-settings['ssa']),atmos_thickness)
        
    elif casesel == 'e1':
        alttau = settings['alttau']
        medium['scatterer'][:,0] = 1.0 
        thicknesses = np.diff(alttau[::2,0])
        for i in range(0,98,2):
            thick_idx = i//2
            medium['scattering_cross_section'][i,0,0] = tau_to_xsec(alttau[i,1],thicknesses[thick_idx])
            medium['scattering_cross_section'][i+1,0,0] = medium['scattering_cross_section'][i,0,0]
        medium['scattering_cross_section'][-1,0,0] = medium['scattering_cross_section'][-2,0,0]
        
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel_parameter'] = np.ones((1,n_scatterer))# * settings['rayleigh_depol']
    elif casesel in ['e2','e6']:
        alttausca = settings['alttausca']
        alttauabs = settings['alttauabs']
        medium['scatterer'][:,0] = 1.0
        medium['absorber'][:,0] = 1.0
        thicknesses = np.diff(alttausca[::2,0])
        for i in range(0,98,2):
            thick_idx = i//2
            medium['scattering_cross_section'][i,0,0] = tau_to_xsec(alttausca[i,1],thicknesses[thick_idx])
            medium['scattering_cross_section'][i+1,0,0] = medium['scattering_cross_section'][i,0,0]
            medium['absorbing_cross_section'][i,0,0] = tau_to_xsec(alttauabs[i,1],thicknesses[thick_idx])
            medium['absorbing_cross_section'][i+1,0,0] = medium['absorbing_cross_section'][i,0,0]
        medium['scattering_cross_section'][-1,0,0] = medium['scattering_cross_section'][-2,0,0]
        medium['absorbing_cross_section'][-1,0,0] = medium['absorbing_cross_section'][-2,0,0]    
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel_parameter'] = np.ones((1,n_scatterer))# * settings['rayleigh_depol']
    elif casesel in ['e3','e4','e5']:
        # rayleigh scattter, molecular abs., aerosol
        alttausca = settings['alttausca']
        alttauabs = settings['alttauabs']
        medium['scatterer'][:,:] = 1.0
        medium['absorber'][:,:] = 1.0
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel'][1] = 1.0
        medium['scatterer_kernel_parameter'] = np.ones((1,n_scatterer))
        thicknesses = np.diff(alttausca[::2,0])
        
        #first the rayleigh
        for i in range(0,98,2):
            thick_idx = i//2
            medium['scattering_cross_section'][i,0,0] = tau_to_xsec(alttausca[i,1],thicknesses[thick_idx])
            medium['scattering_cross_section'][i+1,0,0] = medium['scattering_cross_section'][i,0,0]
            medium['absorbing_cross_section'][i,0,0] = tau_to_xsec(alttauabs[i,1],thicknesses[thick_idx])
            medium['absorbing_cross_section'][i+1,0,0] = medium['absorbing_cross_section'][i,0,0]
        medium['scattering_cross_section'][-1,0,0] = medium['scattering_cross_section'][-2,0,0]
        medium['absorbing_cross_section'][-1,0,0] = medium['absorbing_cross_section'][-2,0,0]    
        # then, the aerosol layer
        if not casesel == 'e4':
            aer_idx_start = np.argmin(np.abs(altitudes - settings['aer_lims'][0]))
            aer_idx_end = np.argmin(np.abs(altitudes - settings['aer_lims'][1]))
            layer_thickness = altitudes[aer_idx_end] - altitudes[aer_idx_start]
            medium['scattering_cross_section'][aer_idx_start:aer_idx_end,0,1] = tau_to_xsec(settings['ssa'] * settings['tau_sca'],layer_thickness)
            medium['absorbing_cross_section'][aer_idx_start:aer_idx_end,0,1] = tau_to_xsec((1-settings['ssa']) * settings['tau_sca'],layer_thickness)
        else:
            medium['scatterer_kernel'][2] = 2.0
            for j in range(1,3):
                aer_idx_start = np.argmin(np.abs(altitudes - settings['aer_lims'][j-1][0]))
                aer_idx_end = np.argmin(np.abs(altitudes - settings['aer_lims'][j-1][1]))
                layer_thickness = altitudes[aer_idx_end] - altitudes[aer_idx_start]
                medium['scattering_cross_section'][aer_idx_start:aer_idx_end,0,j] = tau_to_xsec(settings['ssa'][j-1] * settings['tau_sca'][j-1],layer_thickness)
                medium['absorbing_cross_section'][aer_idx_start:aer_idx_end,0,1] = tau_to_xsec((1-settings['ssa'][j-1]) * settings['tau_sca'][j-1],layer_thickness)
    #emitters and emissivities
    medium['emitter'] = np.zeros((n_medium_positions,n_emitter))
    medium['medium_emissivity'] = np.zeros((n_medium_positions,n_wl,4,n_emitter))
    medium['emitter_kernel'] = np.zeros((n_emitter,))
    medium['emitter_kernel_parameter'] = np.zeros((n_emitter,))
    
    #refractive index
    medium['refractive_index'] = np.ones((n_medium_positions,n_wl))
    
    # The interpolation functions
    # 1: linear interpolation
    medium['interpolation_function'] = np.ones((n_medium_positions,))
    medium['interpolation_parameter'] = np.zeros((n_medium_positions,))
    
    # INSTRUMENT DEFINITIONS
    instrument = {}
    
    n_los = len(salts) * len(vzas) * len(vaas)
    if casesel == 'e6':
        n_los = 61 * 61
    instrument['position'] = np.zeros((n_los,3))
    instrument['view_vector'] = np.zeros((n_los,3))
    
    z_ax = np.array([0.0, 0.0, 1.0])
    y_ax = np.array([0.0, 1.0, 0.0])
    x_ax = np.array([1.0, 0.0, 0.0])
    nadir = np.array([0.0,0.0,-1.0])
    idx = 0
    vza_list = []
    vaa_list = []
    if casesel != 'e6':
        sza = szas[idx_sza]
        saa = saas[0]
        # note: the minus sign in the y_ax on the following line redefines the polarization directions
        # to match MYSTIC.
        solar_dir = arsca.tf.solar_direction(z_ax,x_ax,sza,saa)
        for idx_salt, sensor_alt in enumerate(salts):
            sat_pos = np.array([0.0, 0.0, sensor_alt + R_earth])
            if sensor_alt > 60: # the top of atmosphere case
                view_vec_0 = nadir
                zen_ang_dir = -1.0
            else: # the bottom of atmos. case
                view_vec_0 = -nadir
                zen_ang_dir = 1.0
            for idx_vza, vza in enumerate(vzas):
                for idx_vaa, vaa in enumerate(vaas):
                    view_vec_zen = arsca.tf.arb_rotation(view_vec_0,zen_ang_dir * vza * np.pi/180,y_ax)
                    view_vec_zenazi = arsca.tf.arb_rotation(view_vec_zen,-1 * vaa * np.pi/180,z_ax)
                    #print(vza,vaa,sensor_alt,view_vec_zenazi) #ok!
                    if idx_vza == 8 and idx_sza == 1 and idx_vaa == 18:
                        print(sat_pos,view_vec_zenazi,solar_dir)
                    instrument['position'][idx,:] = sat_pos
                    instrument['view_vector'][idx,:] = view_vec_zenazi
                    vza_list.append(vza)
                    vaa_list.append(vaa)
                    idx += 1
    else:
        # the geostationary view of Earth's disk
        sza = szas[idx_sza]
        saa = saas[0]
        solar_dir = arsca.tf.solar_direction(z_ax,x_ax,sza,saa)
        sensor_alt = 3e5 # km
        sat_pos = np.array([[0.0, 0.0, sensor_alt + R_earth]])
        view_vecs = arsca.inst.camera_fov([61,61], [2.4*np.pi/180, 2.4*np.pi/180], -z_ax, y_ax)
        instrument['view_vector'][:,:] = view_vecs
        instrument['position'][:,:] = np.repeat(sat_pos,n_los,axis=0)
    # BOUNDARY DEFINITIONS
    boundary = {}
    boundary['shape'] = np.array([1,1]) #both are spherical surfaces
    boundary['parameter'] = np.array([R_earth,R_earth + atmos_thickness]) #the radius of the spheres
    
    #reflections from the boundary (0=pass through,1=lambertian,2=semispecular,3=brdf)
    # TODO: Hide this to the background
    #boundary['reflection_kernel'] = np.array([3,0])
    boundary['reflection_kernel'] = np.array([1,0])
    boundary['reflection_kernel_parameter'] = np.array([[settings['albedo'],0.0]])
    
    source = {}
    
    source['input_wavelength'] = wl
    source['output_wavelength'] = wl
    source['type'] = np.array([0])
    source['incident_direction'] = np.array([
            [np.nan, np.nan, np.nan]])
    source['incident_direction'][0,:] = -solar_dir
    source['position'] = np.nan * np.ones_like(source['incident_direction'])
    
    inc_stokes = np.reshape(np.array([1.0,0.0,0.0,0.0]),(1,1,4))
    inc_stokes = np.repeat(inc_stokes,n_wl,axis=1)
    source['incident_stokes'] = np.repeat(inc_stokes,1,axis=0)
    source['parameter'] = np.array([0.0])
    source['source_angular_radius'] = np.array([1.0])
    
    import time
    input_fname = arsca.io.create_simulator_input(medium,instrument,source,boundary)
    
    siro_start = time.time()
    radiance_siro = arsca.simu.run('siro',input_fname)
    radiance_siro_sum = np.sum(radiance_siro[:,:,:],axis=2)
    siro_time = time.time() - siro_start
    if not casesel == 'e6':
        stokes = radiance_siro_sum.reshape((len(salts),1,len(vzas),len(vaas),4))
    else:
        stokes = radiance_siro_sum.reshape((1,1,61,61,4))
    #print(stokes.shape)
    #stokes = rotate_IQUV(vaas, stokes[:,0,:,:,:])
    stokes = rotate_IQUV(vzas, stokes[:,0,:,:,:])
    save_results_to_nc(casesel, idx_sza, stokes)
    