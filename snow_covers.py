#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:35:11 2020

@author: mikkonea
"""

# Tässä studyssä on tärkeää todeta, että vaikka Siro on validoitu UV-aallonpituuksilla
# niin ei voida suoraan sanoa, että ne matchaa VIS-NIR-SWIR alueilla, vaikka siellä
# ilmiöt ovatkin yksinkertaisempia. Vai voidaanko, ehkä. Sironnat on pienempiä.
# En tiiä tarviikohan limb-mittausta tehdä. 

# Nadir/Glint/Limb
# No scattering/Rayleigh/Rayleigh+Aerosol
# Kaikki kolme bändiä
# Sopivat vakio-albedot

# 250 nm, (sulfaatti), n_real = 1.61, 1.58, 1.57 in each bands respectively
# 

# pitääkö plotata esittelydataa, esim. että tämmönen oli nyt se käytetty kaasuprofiili
# 

import numpy as np
import netCDF4
from scipy.interpolate import interp1d
import scipy.optimize
import arsca
import pickle

#validation_fname = 'snow_surf_comparison_new_art.pick'
#validation_fname = 'snow_surf_testing.pick'
validation_fname = 'snow_surf_v5.pick'
#validation_fname = 'water_surf_simulation.pick'

# joitain juttuja tehdään jokaisessa ongelmanmäärittelyssä
# oisko niihin jotain työkaluja nopeuttamiseksi?

# pitää olla jokin layout, jonka pohjalta voi rakentaa inputin

#oco2 lims: 0.758-0.772, 1.594-1.619, 2.042-2.082

# Noin 15 viewing anglea, muutama solar zenith angle, yhdenlainen ilmakehä, 4 eri lumipintaa

# TODO: Laskeppa AOD!
# TODO: Troposfääriaerosoleja nadiirin ja glinttiin (n 2 km korkeudessa)
# TODO: OCO-2:n luminen targetti voisi olla dataa parhaaseen malliin! 
import time
total_start = time.time()
band_names = ['o2a','wco2','sco2']
surf_names = ['art','loose_snow','new_snow','rough_snow','lambertian']
surf_filenames = ['art_snow_msbrdf.nc4', 'snortex_snow1.nc4', 'new_snow_v0.nc4','greenland_snow6.nc4','greenland_snow6.nc4']
compute_siro = False
for band_choice in range(3):
    for surf_choice in range(5):
        solar_zens = [36,54,69]
        solar_azis = [0,0,0]
        arsca.simu.change_raysca_settings('brdf_file',surf_filenames[surf_choice])
        
        try:
            with open(validation_fname,'rb') as fhandle:
                radiance_dict = pickle.load(fhandle)
        except FileNotFoundError:
            print("No file found. Setting up a new one...")
            radiance_dict = {}
            for sname in surf_names:
                radiance_dict[sname] = {}
                for bname in band_names:
                    radiance_dict[sname][bname] = np.array([np.nan])
        
        if not np.isnan(radiance_dict[surf_names[surf_choice]][band_names[band_choice]][0]).all():
            #raise ValueError("This configuration already computed!") # this is to "return" or halt
            print("This configuration already computed!")
            continue
        
        band_lims = [[755.0,775.0],[1590.0,1620.0],[2040.0, 2080.0]]
        #band_albs = [0.9,0.1,0.02] # solar zenith 55 from the picture
        #band_albs = [0.87651711, 0.17815365, 0.05095223] # from the model
        band_albs = [0.88, 0.18, 0.05] # from the model
        
        arsca.set_case("snow_surf_simu_%s" % band_names[band_choice])
        arsca.set_configuration("initial_snow_simus_%s" % band_names[band_choice])
        
        wl_range = band_lims[band_choice] # nm
        wl_range = band_lims[band_choice] # nm
            
        wn_step = 0.01
        wn_range = arsca.tf.wl2wn(wl_range)
        wn_range[0] = wn_range[0] - wn_step #widen it a bit
        wn_range[1] = wn_range[1] + wn_step
        n_wn = np.array(np.arange(wn_range[0],wn_range[1],wn_step)).size
        n_wl = n_wn
        wl = np.linspace(wl_range[0],wl_range[1],n_wl)
        wn_coarse = arsca.tf.wl2wn(wl)
        wn_reso = np.max(np.abs(np.diff(wn_coarse)))
            
        R_earth = 6371.0 # km
        
        altitudes = np.linspace(0,70,20)
        
        n_medium_positions = altitudes.size
        n_coordinate = 3
        
        # MEDIUM DEFINITIONS
        
        medium = {}
        
        #medium['position'] is now a n_medium_positions x 3
        altitudes_ = altitudes.reshape((n_medium_positions, 1)) + R_earth
        position = np.zeros((n_medium_positions, n_coordinate))
        position[:,[0]] = altitudes_
        medium['position'] = position
        
        n_emitter = 0 #emission will be added later
        if band_choice == 0:
            n_absorber = 1
        else:
            n_absorber = 3
            
        if band_choice == 0: # O2A band
            gases_mav = ['1o2']
            gases_hapi = ['O2']
            gases_hapi_id = [7]
        elif band_choice == 1: # WCO2 band
            gases_mav = ['1co2','1h2o','1ch4']
            gases_hapi = ['CO2','H2O','CH4'] #these are the local file names for the absorption lines
            gases_hapi_id = [2,1,6] #these are the values used by HITRAN
        elif band_choice == 2: # SCO2 band
            gases_mav = ['1co2','1h2o','1n2o']
            gases_hapi = ['CO2','H2O','N2O']
            gases_hapi_id = [2,1,4]
        
        atmos = arsca.io.read_atmos('./datafiles/atmos/so20140319.mav', gases_mav)
        
        n_scatterer = 1
        #scatterers and their cross-sections
        medium['scatterer'] = np.zeros((n_medium_positions,n_scatterer))
        medium['scattering_cross_section'] = np.zeros((n_medium_positions,n_wl,n_scatterer))
        xsec_sca = np.zeros((n_medium_positions,n_wl,n_scatterer))
             
        sca_interpfun = interp1d(atmos['Height'],atmos['Density'])
        medium['scatterer'][:,0] = sca_interpfun(altitudes).reshape((altitudes.size,))
        rayleigh_xsec = arsca.scatter.rayleigh_xsec(wl,extrapolate=True)
        xsec_sca[:,:,0] = np.repeat(rayleigh_xsec.reshape((1,n_wl)),n_medium_positions,axis=0)
            
        medium['scattering_cross_section'] = xsec_sca
        
        # zero stands for rayleigh-scattering
        medium['scatterer_kernel'] = np.zeros((n_scatterer,))
        medium['scatterer_kernel_parameter'] = np.zeros((1,n_scatterer))
        
        #generation of the absorbtion cross-sections
        
        medium['absorber'] = np.zeros((n_medium_positions,n_absorber))
        medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber))
        
        wn_step = 0.005
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
            
            xsec_conv = xsec
            #for i in range(xsec.shape[0]):
            #    wn_conv,xsec_conv[i,:],_,_,_ = arsca.xsec.hapi.convolveSpectrumSame(wn,xsec[i,:],Resolution=wn_reso,SlitFunction=arsca.xsec.hapi.SLIT_GAUSSIAN)
            
            # here the cross-section gets flipped to wavelengths
            xsec_abs[:,:,abs_idx] = xsec_conv[:,::-1]
            wl_xsec_dense = arsca.tf.wn2wl(wn)
            
            xsec_interpfun = interp1d(wl_xsec_dense,xsec_abs[:,:,abs_idx],axis=1)
            
            xsec_column = xsec_interpfun(wl)
        
            # The reference to scatterers on the row below is merely because gas_column 
            # is fraction of all particles and absorber needs to be numbers density
            medium['absorber'][:,abs_idx] = gas_column * medium['scatterer'][:,0].ravel()
            
            medium['absorbing_cross_section'][:,:,abs_idx] = xsec_column
        
        
        #emitters and emissivities
        medium['emitter'] = np.zeros((n_medium_positions,n_emitter))
        medium['medium_emissivity'] = np.zeros((n_medium_positions,n_wl,4,n_emitter))
        medium['emitter_kernel'] = np.zeros((n_emitter,))
        medium['emitter_kernel_parameter'] = np.zeros((n_emitter,))
        
        #refractive index
        medium['refractive_index'] = np.ones((n_medium_positions,n_wl))
        
        # The interpolation functions
        medium['interpolation_function'] = np.ones((n_medium_positions,))
        medium['interpolation_parameter'] = np.zeros((n_medium_positions,))
        
        # INSTRUMENT DEFINITIONS
        
        instrument = {}
        # The altitude is necessary, but in this case it doesn't matter as long as
        # we're above the atmosphere. In the case of 90 degrees of viewing zenith
        # angle, we'll need to have a minimum "altitude" (in quotes because it gets 
        # rotated from zenith) of about 950 km to stay above the atmosphere of 
        # 70 km.
        sat_altitude = 1705.0 # Some height for good measure.
        
        #los_angs = np.linspace(0,85,18) # 0, 5, 10, ..., 80, 85.
        los_angs = np.linspace(0,85,18) # 0, 5, 10, ..., 80, 85.
        #los_angs = np.linspace(0,85,3)
        n_los = los_angs.size
        
        instrument['position'] = np.zeros((n_los,3))
        instrument['view_vector'] = np.zeros((n_los,3))
        
        z_ax = np.array([0.0, 0.0, 1.0])
        y_ax = np.array([0.0, 1.0, 0.0])
        x_ax = np.array([1.0, 0.0, 0.0])
        nadir = np.array([-1.0,0.0,0.0])
        
        for i_los in range(n_los):
            view_degs = los_angs[i_los]
            sat_view_dir = arsca.tf.arb_rotation(x_ax,-view_degs * np.pi/180,z_ax)
            instrument['position'][i_los,:] = sat_altitude * sat_view_dir + np.array([R_earth,0.0,0.0])
            # in this case, the satellite is just at 428 km altitude
            instrument['view_vector'][i_los,:] = -sat_view_dir
                
        
        # BOUNDARY DEFINITIONS
        boundary = {}
        boundary['shape'] = np.array([1,1]) #both are spherical surfaces
        boundary['parameter'] = np.array([6371.0,6371.0 + 70.0]) #the radius of the spheres
        #reflections from the boundary (0=pass through,1=lambertian,2=semispecular,3=brdf)
        # TODO: Hide this to the background
        boundary['reflection_kernel'] = np.array([3,0])
        if surf_choice == 4:   
            boundary['reflection_kernel'] = np.array([1,0])
        boundary['reflection_kernel_parameter'] = np.array([[band_albs[band_choice],0.0]])
        
        source = {}
        
        n_source = len(solar_zens)
        source['input_wavelength'] = wl
        source['output_wavelength'] = wl
        source['type'] = np.zeros((n_source,))
        source['incident_direction'] = np.nan * np.ones((n_source,3))
        for i_source in range(n_source):
            solar_zen_degs = solar_zens[i_source]
            solar_azi_degs = solar_azis[i_source]
            solar_dir = arsca.tf.solar_direction(x_ax,y_ax,solar_zen_degs,solar_azi_degs)
            source['incident_direction'][i_source,:] = -solar_dir
            
        source['position'] = np.nan * np.ones_like(source['incident_direction'])
        
        inc_stokes = np.reshape(np.array([1.0,0.0,0.0,0.0]),(1,1,4))
        inc_stokes = np.repeat(inc_stokes,n_wl,axis=1)
        source['incident_stokes'] = np.repeat(inc_stokes,n_source,axis=0)
        source['parameter'] = 0.0 * np.ones((1,n_source))
        source['source_angular_radius'] = 1.0 * np.ones((n_source,))
        

        input_fname = arsca.io.create_simulator_input(medium,instrument,source,boundary)
        
        if compute_siro:
            siro_start = time.time()
            radiance_siro = arsca.simu.run('siro',input_fname)
            siro_time = time.time() - siro_start
        #    rads_siro.append(radiance_siro)
        #else:
        #    rads_siro.append(np.nan)
        
        raysca_start = time.time()
        radiance_raysca = arsca.simu.run('raysca',input_fname)
        raysca_time = time.time() - raysca_start
    
        #rads_raysca.append(radiance_raysca)
        #print(siro_time)
        print(raysca_time)
    
        radiance_dict[surf_names[surf_choice]][band_names[band_choice]] = [radiance_raysca]
        
        with open(validation_fname,'wb') as fhandle:
            # Writing results to the file.
            pickle.dump(radiance_dict,fhandle)
            
