#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:35:11 2020

@author: mikkonea
"""

import sys
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
import scipy.optimize
import arsca
import pickle

step_siro = 1.0
#noph_siro = 100

casesel = sys.argv[1]
noph_siro = int(sys.argv[2])
#casesel = 'd1'

vzas = [0, 9, 18, 26, 34, 41, 48, 54, 60, 65, 70, 74, 78, 81, 84, 86, 88, 89, 90]
vaas = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
szas = [30, 60, 80, 87, 90, 93, 96, 99]
saas = [0]
salts= [0 + 1e-4, 120 - 1e-4] # sensor altitude, kilometres
atmos_thickness = 120
R_earth = 6371

siro_custom_settings = {
    'noph' : noph_siro,
    'ratm' : R_earth + atmos_thickness,
    'step' : step_siro}

if casesel == 'd1':
    # single layer Rayleigh
    settings = {'n_sca' : 1,
                'n_abs' : 1, # needs to be one, but filled with zeroes
                'tau_sca' : 0.5,
                'wavelength' : 450,
                'albedo' : 0.0,
                'n_lay' : 1,                
                'rayleigh_depol' : 0.03} # currently hardcoded into Siro
elif casesel == 'e1':
    alttau = np.genfromtxt('datafiles/atmos/tau_rayleigh_450nm_usstd.dat',comments='#')
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
    settings = {'n_sca' : 1,
                'n_abs' : 1, 
                'alttau' : alttau,
                'wavelength' : 450,
                'albedo' : 0.0}
elif casesel == 'd3':
    siro_custom_settings['AER_FILENAME'] = "'input/miefiles/%s'"
    
else:
    print('bad case selection: %s' % casesel)
    halt
    
outfilename = 'iprt_phase3_Siro.nc'

def save_results_to_nc(casesel,idx_sza,stokes):
    var_names = ['%s_sza','%s_saa','%s_vza','%s_vaa','%s_zout','radiance_%s']
    var_names = [vname % casesel for vname in var_names]
    stokes_in = np.reshape(stokes,(len(salts),1,1,len(vzas),len(vaas),4))
    with netCDF4.Dataset(outfilename,'a') as ds:
        ds[var_names[-1]][:,idx_sza,0,:,:,:] = stokes_in
        
def tau_to_xsec(tau,thickness):
    # converts optical thicknesses in km to xsec in cm2
    # assuming the number density to be 1/cm3
    return tau / thickness / 1e5
    
for idx_sza in range(len(szas)):
    #for idx_sza in [0]:    
    compute_siro = True
    
    arsca.simu.create_siro_settings(siro_custom_settings)
    
    arsca.set_case("iprt_%s" % casesel)
    arsca.set_configuration("base_conf_%s" % casesel)
    
    n_wn = 1
    wl = np.array([settings['wavelength']])
    
    n_wl = n_wn
    if casesel == 'd1':
        n_lev = settings['n_lay'] + 1
        altitudes = np.linspace(0,atmos_thickness,n_lev)
    elif casesel == 'e1':
        altitudes = settings['alttau'][:,0]
    
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
    
    if casesel == 'd1':
        # single layer Rayleigh
        
        medium['scatterer'][:,0] = 1.0 
        medium['scattering_cross_section'] = tau_to_xsec(settings['tau_sca'],atmos_thickness)
        
        # zero stands for rayleigh-scattering
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel_parameter'] = np.ones((1,n_scatterer)) * settings['rayleigh_depol']
        medium['absorber'] = np.zeros((n_medium_positions,n_absorber))
        medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber))
        
    elif casesel == 'd3':
        medium['scatterer_kernel'][1] = 1 # this is to use the aerosols
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
        medium['absorber'] = np.zeros((n_medium_positions,n_absorber))
        medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber))
    
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
    instrument['position'] = np.zeros((n_los,3))
    instrument['view_vector'] = np.zeros((n_los,3))
    
    z_ax = np.array([0.0, 0.0, 1.0])
    y_ax = np.array([0.0, 1.0, 0.0])
    x_ax = np.array([1.0, 0.0, 0.0])
    nadir = np.array([-1.0,0.0,0.0])
    idx = 0
    if casesel != 'e6':
        sza = szas[idx_sza]
        saa = saas[0]
        solar_dir = arsca.tf.solar_direction(x_ax,y_ax,sza,saa)
        
        
        for idx_salt, sensor_alt in enumerate(salts):
            sat_pos = np.array([sensor_alt + R_earth, 0.0, 0.0])
            if sensor_alt > 60: # the top of atmosphere case
                view_vec_0 = nadir
                zen_ang_dir = -1.0
            else: # the bottom of atmos. case
                view_vec_0 = -nadir
                zen_ang_dir = 1.0
            for idx_vza, vza in enumerate(vzas):
                for idx_vaa, vaa in enumerate(vaas):
                    view_vec_zen = arsca.tf.arb_rotation(view_vec_0,zen_ang_dir * vza * np.pi/180,z_ax)
                    view_vec_zenazi = arsca.tf.arb_rotation(view_vec_zen,-1 * vaa * np.pi/180,x_ax)
                    #print(vza,vaa,sensor_alt,view_vec_zenazi) #ok!
                    instrument['position'][idx,:] = sat_pos
                    instrument['view_vector'][idx,:] = view_vec_zenazi
                    idx += 1
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
    radiance_siro_sum = np.sum(radiance_siro,axis=2)
    siro_time = time.time() - siro_start

    stokes = radiance_siro_sum.reshape((len(salts),1,len(vzas),len(vaas),4))
    save_results_to_nc(casesel, idx_sza, stokes)