#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:35:11 2020

@author: mikkonea
"""



import numpy as np
import netCDF4
from scipy.interpolate import interp1d
import scipy.optimize
import arsca
import os.path
from multiprocessing import Pool
import os
img_sizes = [10,20,50,100,2,5,30,40]
anim_vzas = np.linspace(-80,80,33)
#anim_vzas = [-0.5]
i_img = 0
i_vza = 16

if True:
    plume_choice = 1
    source_choice = 1
    plume_end_spreads = [0.2, 0.4]
    plume_end_spread = plume_end_spreads[plume_choice]
    ground_pixel_size = 0.06
    image_size = np.array([img_sizes[i_img],img_sizes[i_img]])
    #image_size = np.array([15,15])
    sat_altitude = 417.0
    # ISS: 413 - 422 km
    half_image = np.array([100,100]) / 2 * ground_pixel_size
    image_angle = 2 * np.arctan(half_image / sat_altitude)# * 180 / np.pi
    n_los = np.prod(image_size)
    
    
    #source_choice = 3
    emission_prop = 0.02 #
    sources = np.array([1,13,25]) # Mt/a
    
    sources = sources * 1e6 * 1e3 * 1e3 / (365 * 24 * 3600) #g/s
    
    r_g = 0.07e-6 * 1e2# fine mode strongly absorbing in cm
    bc_dens = 2.0 #g/cm^3
    bc_vol = r_g ** 3 * np.pi * (4/3)
    bc_mass = bc_dens * bc_vol
    
    co2_molar_mass = 44.01 # g/mol
    N_avo = 6.022e23
    co2_mass = co2_molar_mass / N_avo # g
    co2_sources = (1 - emission_prop) * sources / co2_mass # 1/s
    bc_sources = emission_prop * sources / bc_mass
    
    
    plume_size = 30
    wind_speed = 5 # m/s
    wind_direction = np.array([-np.sqrt(0.5), -np.sqrt(0.5)]) #normalize this!
    plume_elevation = 1.0 # km in the total length
    #plume_end_spread = 0.4 # km at the total length
    plume_start_spread = 0.1
    plume_length = 3 # km
    s_in_full_plume = plume_length / (1e-3 * wind_speed)
    co2_in_full_plume = co2_sources[source_choice] * s_in_full_plume
    bc_in_full_plume = bc_sources[source_choice] * s_in_full_plume
    def approx_plume_vol(plume_length,plume_start_spread,plume_end_spread):
        # volume of a truncated cone
        r = plume_start_spread
        R = plume_end_spread
        height = plume_length
        vol = 1/3 * np.pi * height * (r**2 + r*R + R**2)
        return vol
    plume_vol = approx_plume_vol(plume_length, plume_start_spread, plume_end_spread)
    avg_co2_dens_in_full_plume = 0.65 * co2_in_full_plume / plume_vol
    avg_bc_dens_in_full_plume = 0.65 * bc_in_full_plume / plume_vol
    s_of_emission_in_one_element = s_in_full_plume / plume_size
    #s_of_emission_in_one_pixel = ground_pixel_size / (1e-3 * wind_speed)
    #s_of_emission_in_cm = 1e-5 / (1e-3 * wind_speed)
    
    #overlap = 0.5
    
    plume_argument = np.zeros((plume_size,)) # the sigma in km
    #seconds_per_element = np.zeros((plume_size,))
    def plume_spread_fun(i):
        return plume_start_spread + i/plume_size * (plume_end_spread - plume_start_spread)
    for i in range(plume_size):
        plume_argument[i] = plume_spread_fun(i)
        #seconds_per_element[i] = s_in_full_plume * (plume_argument[i] / plume_length)
    km3_to_cm3 = 1e15
    
    R_earth = 6371.0 # km
    plume_start = np.array([R_earth,0,0])
    plume_step = np.array([plume_elevation, plume_length * wind_direction[0], plume_length * wind_direction[1]])
    #plume_step = plume_step / np.linalg.norm(plume_step)
    plume_step = plume_step / plume_size
    plume_position = np.zeros((plume_size,3)) 
    plume_inclusion = np.zeros((plume_size,)) # the center point in ppm (roughly)
    plume_aerosols = np.zeros((plume_size,)) 
    for i in range(plume_size):
        #plume_position[i,:] = i * (plume_step * overlap * plume_argument[i]) + plume_start
        plume_position[i,:] = i * plume_step + plume_start
        plume_inclusion[i] = avg_co2_dens_in_full_plume / km3_to_cm3
        plume_aerosols[i] = avg_bc_dens_in_full_plume / km3_to_cm3
        #plume_inclusion[i] = s_of_emission_in_one_element * co2_sources[source_choice] / km3_to_cm3
        #plume_aerosols[i] = s_of_emission_in_one_element * bc_sources[source_choice] / km3_to_cm3
        #plume_inclusion[i] = s_of_emission_in_cm * co2_sources[source_choice] / km_to_cm
        #plume_aerosols[i] = s_of_emission_in_cm * bc_sources[source_choice] / km_to_cm
        #plume_inclusion[i] = seconds_per_element[i] * co2_sources[source_choice] / km3_to_cm3
        #plume_aerosols[i] = seconds_per_element[i] * bc_sources[source_choice] / km3_to_cm3
    
    solar_zen = 77
    solar_azi = 124
    
    band_names = ['o2a','wco2','sco2']
    #oco2 lims: 0.758-0.772, 1.594-1.619, 2.042-2.082
    import pickle

    if True:
        #band_choice = int(sys.argv[1])
        band_choice = 2
        #band_lims = [[747.0,773.0],[1590.0,1675.0],[1990.0, 2095.0]]
        band_lims = [[747.0,773.0],[1590.0,1675.0],[1990.0, 2000.0]]
        band_wls = [217, 283, 300] # about 0.3 nm
        band_wls = [217, 283, 14] # 
        aerosol_names = ['','SWIR1-TSA.dat','SWIR2-TSA.dat']
        arsca.simu.change_raysca_settings('aerosol_file',aerosol_names[band_choice])
        arsca.simu.change_raysca_settings('single_scattering','False')
        arsca.simu.change_raysca_settings('parallelization','False')
        #arsca.simu.change_raysca_settings('brdf_file','reflectance_secunda_09092022.nc4')
        atm_choice = 0
        atmos = arsca.io.read_nc_atmos('./datafiles/atmos/snowite_atmos_v3.nc4',atm_choice,'ggg')
        
        arsca.set_case("plume_test_tursca_release_small_pixel_%s" % (band_names[band_choice]))
        arsca.set_configuration("initial_%dx%d_vza%d" % (img_sizes[i_img],img_sizes[i_img],anim_vzas[i_vza]))
        
        wl_range = band_lims[band_choice] # nm
        n_wl_coarse = band_wls[band_choice]
        wl_coarse = np.linspace(wl_range[0],wl_range[1],n_wl_coarse)
        wn_coarse = arsca.tf.wl2wn(wl_coarse)
        wn_reso = np.max(np.abs(np.diff(wn_coarse)))
        
        n_sel_wl = 10
        R_earth = 6371.0 # km
        altitudes = atmos['Height']
        altitudes = np.arange(altitudes.size)
        
        n_medium_positions = altitudes.size + plume_size
        n_coordinate = 3
        
        # MEDIUM DEFINITIONS
        
        medium = {}
        
        #medium['position'] is now a n_medium_positions x 3
        altitudes_ = altitudes.reshape((altitudes.size, 1)) + R_earth
        n_alt = altitudes_.size
        position = np.zeros((n_medium_positions, n_coordinate))
        position[:n_alt,[0]] = altitudes_
        position[n_alt:,:] = plume_position
        
        medium['position'] = position
        
        n_emitter = 0 #emission will be added later
        if band_choice == 0:
            n_absorber = 1
        elif band_choice == 1:
            n_absorber = 3
        else:
            n_absorber = 2
            
        if band_choice == 0: # O2A band
            gases_mav = ['o2']
            gases_hapi = ['O2']
            gases_hapi_id = [7]
        elif band_choice == 1: # WCO2 band
            gases_mav = ['co2','h2o','ch4']
            gases_hapi = ['CO2','H2O','CH4'] #these are the local file names for the absorption lines
            gases_hapi_id = [2,1,6] #these are the values used by HITRAN
        elif band_choice == 2: # SCO2 band
            #gases_mav = ['co2','h2o','n2o']
            #gases_hapi = ['CO2','H2O','N2O']
            #gases_hapi_id = [2,1,4]
            gases_mav = ['co2','h2o']
            gases_hapi = ['CO2','H2O'] #these are the local file names for the absorption lines
            gases_hapi_id = [2,1] #these are the values used by HITRAN
        
        #atmos = arsca.io.read_atmos('./files_used_in_examples/scout/so20140319.mav', gases_mav)
        
        #atmos = arsca.io.read_nc_atmos('./snowite_atmos_v3.nc4',atm_choice)
            
        n_scatterer = 2
        
        n_absorber_tot = n_absorber + 1
        # zero stands for rayleigh-scattering
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel'][1] = 1
        medium['scatterer_kernel_parameter'] = np.zeros((1,n_scatterer))
        
        medium['absorber'] = np.zeros((n_medium_positions,n_absorber_tot))
        
        wn_step = 0.02
        wn_range = arsca.tf.wl2wn(wl_range)
        wn_range[0] = wn_range[0] - wn_step #widen it a bit
        wn_range[1] = wn_range[1] + wn_step
        n_wn = np.array(np.arange(wn_range[0],wn_range[1],wn_step)).size
        
        xsec_abs = np.zeros((altitudes.size,n_wn,n_absorber))
        
        pres_interpfun = interp1d(atmos['Height'],atmos['Pres'],fill_value="extrapolate")
        temp_interpfun = interp1d(atmos['Height'],atmos['Temp'],fill_value="extrapolate")
        pres = pres_interpfun(altitudes)
        temp = temp_interpfun(altitudes)
        
        for abs_idx in range(n_absorber):
            gas_interpfun = interp1d(atmos['Height'],atmos[gases_mav[abs_idx]],fill_value="extrapolate")
            gas_column = gas_interpfun(altitudes)
            p_s = pres * gas_column
            wn, xsec = arsca.xsec.generate_xsec(wn_range,wn_step,gases_hapi[abs_idx],
                                                gases_hapi_id[abs_idx],temp,pres,p_s,
                                                func_selection='Voigt')
            n_wl = wn.size
            #n_wl = n_sel_wl
            
            wl = arsca.tf.wn2wl(wn)
            if 'absorbing_cross_section' not in medium.keys():
                medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber_tot))
            if 'scattering_cross_section' not in medium.keys():
                #scatterers and their cross-sections
                medium['scatterer'] = np.zeros((n_medium_positions,n_scatterer))
                medium['scattering_cross_section'] = np.zeros((n_medium_positions,n_wl,n_scatterer))
                xsec_sca = np.zeros((n_medium_positions,n_wl,n_scatterer))
                     
                sca_interpfun = interp1d(atmos['Height'],atmos['Density'])
                medium['scatterer'][:n_alt,0] = sca_interpfun(altitudes).reshape((altitudes.size,))
                rayleigh_xsec = arsca.scatter.rayleigh_xsec(wl,extrapolate=True)
                xsec_sca[:n_alt,:,0] = np.repeat(rayleigh_xsec.reshape((1,n_wl)),n_alt,axis=0)
                    
                medium['scattering_cross_section'] = xsec_sca
            
            
            # here the cross-section gets flipped to wavelengths
            #xsec_abs[:,:,abs_idx] = xsec_conv[:,::-1]
            xsec_abs[:,:,abs_idx] = xsec[:,::-1]
            #wl_xsec_dense = arsca.tf.wn2wl(wn_conv)
            #xsec_interpfun = interp1d(wl_xsec_dense,xsec_abs[:,:,abs_idx],axis=1)
            
            xsec_column = xsec_abs[:,:,abs_idx]
            #xsec_column = xsec_interpfun(wl)
        
            # The reference to scatterers on the row below is merely because gas_column 
            # is fraction of all particles and absorber needs to be numbers density
            medium['absorber'][:n_alt,abs_idx] = gas_column * medium['scatterer'][:n_alt,0].ravel()
            
            medium['absorbing_cross_section'][:n_alt,:,abs_idx] = xsec_column
        
        ## THE AEROSOLS
        
        type_choice = 1
        ice_input_folder = './datafiles/aerosols/'
        aer_file = ['extssa-twa.dat','extssa-tsa.dat'][type_choice]
        
        aer_data = np.genfromtxt(ice_input_folder + aer_file)
        aer_wl = aer_data[:,0]
        aer_xsec_sca = aer_data[:,1]
        aer_xsec_abs = aer_data[:,2]
        
        xsec_aer_sca_interpfun = interp1d(aer_wl,aer_xsec_sca,fill_value='extrapolate')
        xsec_aer_abs_interpfun = interp1d(aer_wl,aer_xsec_abs,fill_value='extrapolate')
        
        aod_wl = 765
        avg_ext = xsec_aer_sca_interpfun(aod_wl) + xsec_aer_abs_interpfun(aod_wl)
        
        
        
        medium['absorber'][:n_alt,-1] = np.zeros_like(altitudes)
        medium['scatterer'][:n_alt,1] = np.zeros_like(altitudes)
        
        aer_sca_xsec = xsec_aer_sca_interpfun(wl)
        aer_abs_xsec = xsec_aer_abs_interpfun(wl)
        medium['absorbing_cross_section'][:,:,-1] = np.repeat(aer_abs_xsec.reshape((1,n_wl)),n_medium_positions,axis=0)
        medium['scattering_cross_section'][:,:,1] = np.repeat(aer_sca_xsec.reshape((1,n_wl)),n_medium_positions,axis=0)
        
        # The interpolation functions
        medium['interpolation_function'] = np.ones((n_medium_positions,))
        medium['interpolation_parameter'] = np.zeros((n_medium_positions,))
        
        #handle plume
        for i in range(-plume_size,0):
            medium['absorber'][i,0] = plume_inclusion[i] * np.ones_like(medium['scatterer'][0,0])
            medium['absorbing_cross_section'][i,:,0] = medium['absorbing_cross_section'][0,:,0]
        
            medium['scatterer'][i,0] = np.zeros_like(plume_inclusion[i] / 1e6 * medium['scatterer'][0,0])
            medium['scattering_cross_section'][i,:,0] = np.zeros_like(medium['absorbing_cross_section'][0,:,0])
    
            medium['scatterer'][i,1] = plume_aerosols[i] * np.ones_like(plume_inclusion[i] / 1e6 * medium['scatterer'][0,0])
            medium['absorber'][i,-1] = plume_aerosols[i] * np.ones_like(plume_inclusion[i] / 1e6 * medium['scatterer'][0,0])
            #medium['scattering_cross_section'][i,:,1] = np.zeros_like(medium['absorbing_cross_section'][0,:,0])
    
            medium['interpolation_function'][i] = 2
            medium['interpolation_parameter'][i] = plume_argument[i]
        
        
        #emitters and emissivities
        medium['emitter'] = np.zeros((n_medium_positions,n_emitter))
        medium['medium_emissivity'] = np.zeros((n_medium_positions,n_wl,4,n_emitter))
        medium['emitter_kernel'] = np.zeros((n_emitter,))
        medium['emitter_kernel_parameter'] = np.zeros((n_emitter,))
        
        #refractive index
        medium['refractive_index'] = np.ones((n_medium_positions,n_wl))
        
        
        # INSTRUMENT DEFINITIONS
        
        instrument = {}
        #sat_altitude = 1705.0 # Some height for good measure.
        #print("SAT ALTITUDE TEST!")
        #sat_altitude = 193.0
        instrument['position'] = np.zeros((n_los,3))
        instrument['view_vector'] = np.zeros((n_los,3))
        
        z_ax = np.array([0.0, 0.0, 1.0])
        y_ax = np.array([0.0, 1.0, 0.0])
        x_ax = np.array([1.0, 0.0, 0.0])
        nadir = np.array([-1.0,0.0,0.0])
        
        solar_zen_degs = solar_zen
        solar_azi_degs = solar_azi
        
        
        vza = anim_vzas[i_vza]
        solar_dir = arsca.tf.solar_direction(x_ax,y_ax,solar_zen_degs,solar_azi_degs)
        nadir_mode = False
        if not nadir_mode:
            sat_dir = arsca.tf.arb_rotation(x_ax,vza * np.pi/180,z_ax)
            up_dir = arsca.tf.arb_rotation(y_ax,vza * np.pi/180,z_ax)
        else:
            sat_dir = x_ax
        sat_centerview_dir = -sat_dir
        for i_los in range(n_los):
            sat_altitude_mod = sat_altitude / np.cos(vza / 180 * np.pi)
            sat_altitude_mod = sat_altitude
            instrument['position'][i_los,:] = sat_altitude_mod * sat_dir + np.array([R_earth,0.0,0.0])
            instrument['position'][i_los,:] = instrument['position'][i_los,:] / np.linalg.norm(instrument['position'][i_los,:]) * (sat_altitude + R_earth)
            
        if not nadir_mode: # in glint, this gets rotated
            sat_centerview_dir = arsca.tf.arb_rotation(sat_centerview_dir,vza * np.pi / 180,z_ax)
            
        sat_centerview_dir = np.array([R_earth,0.0,0.0]) - instrument['position'][0,:]
        sat_centerview_dir = sat_centerview_dir / np.linalg.norm(sat_centerview_dir)
        up_dir = arsca.tf.perp_norm_vec(sat_centerview_dir,z_ax)
        instr_dirs = arsca.inst.camera_fov(image_size,
                                           image_angle,
                                           view_center=sat_centerview_dir,
                                           up_vec=up_dir)
        for i_los in range(n_los):
            #instrument['position'][i_los,:] = 
            instrument['view_vector'][i_los,:] = instr_dirs[i_los,:]
        
        # BOUNDARY DEFINITIONS
        boundary = {}
        boundary['shape'] = np.array([1,1]) #both are spherical surfaces
        boundary['parameter'] = np.array([6371.0,6371.0 + 50.0]) #the radius of the spheres
        #reflections from the boundary (0=pass through,1=lambertian,2=semispecular,3=brdf)
        # TODO: Hide this to the background
        #boundary['reflection_kernel'] = np.array([3,0])
        #print("LAMBERTIAN ENABLED")
        boundary['reflection_kernel'] = np.array([2,0])
        boundary['reflection_kernel_parameter'] = np.array([[0.1,0.0]])
        
        source = {}
        
        #Just the glint source
        source['input_wavelength'] = wl#[:n_sel_wl]
        source['output_wavelength'] = wl#[:n_sel_wl]
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
        
        input_fname = arsca.io.create_simulator_input(medium,instrument,source,boundary)
print('Input file created:',input_fname)
