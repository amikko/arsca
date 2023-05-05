#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:35:11 2020

@author: mikkonea
"""


import sys
sys.path.insert(0,'/home/mikkonea/Projects/arsca')
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
import scipy.optimize
import arsca
import pickle

validation_fname = 'validation_study_radiances_oops_all_siros_noeps.pick'
validation_fname = 'validation_study_radiances_oops_all_siros.pick'
validation_fname = 'validation_study_radiances_low_alt_aero_high_alb.pick'
#validation_fname = 'validation_study_radiances_low_alt_aero_zero_alb_many_photons.pick'
#validation_fname = 'vasdflidation_study_radiances_low_alt_aero_zero_alb_many_photons.pick'
validation_fname = 'validation_study_radiances_final_v2.pick'

#oco2 lims: 0.758-0.772, 1.594-1.619, 2.042-2.082

# aerosolit
# sulfaatit: r = 0.125 µm
# 0.765µm: 1.61 - i0.0
# 1.61µm: 1.58 - i0.0
# 2.06µm: 1.57 - i0.0
# black carbon (diesel soot) r = 0.05 µm
# 0.765µm: 1.6835 - i0.30050 
# 1.61µm: 1.8625 - i0.34050
# 2.06µm: 1.9210 - i0.34800


def create_phase_matrices():
    arsca.util.aerosol_input_from_online_mie_calc('sulfate_mie_O2A.txt','sulf_O2A.dat',False)
    arsca.util.aerosol_input_from_online_mie_calc('sulfate_mie_WCO2.txt','sulf_WCO2.dat',False)
    arsca.util.aerosol_input_from_online_mie_calc('sulfate_mie_SCO2.txt','sulf_SCO2.dat',False)
    arsca.util.aerosol_input_from_online_mie_calc('soot_mie_O2A.txt','soot_O2A.dat',False)
    arsca.util.aerosol_input_from_online_mie_calc('soot_mie_WCO2.txt','soot_WCO2.dat',False)
    arsca.util.aerosol_input_from_online_mie_calc('soot_mie_SCO2.txt','soot_SCO2.dat',False)

aerosol_names = ['sulf_O2A.dat', 'sulf_WCO2.dat', 'sulf_SCO2.dat']
noph_siro = 10000
for band_choice in range(3):
    for scatterers in range(1,3):
        for geom_choice in ['limb','nadir','glint']:
            #for geom_choice in ['nadir','glint']:
            #for geom_choice in ['nadir','limb']:
            print(band_choice,scatterers,geom_choice)
            if geom_choice == 'limb':
                arsca.simu.change_raysca_settings('main_beam_step_length',1.0)
                arsca.simu.change_raysca_settings('scattering_step_length',1.0)
            else:
                arsca.simu.change_raysca_settings('main_beam_step_length',0.25)
                arsca.simu.change_raysca_settings('scattering_step_length',1.0)
            solar_zens = [60,35,50]
            #solar_zens = [0,0,60]
            #solar_zens = [35,0,60]
            solar_azis = [120,0,0]
            
            compute_siro = True
            if compute_siro:
                siro_custom_settings = {
                    'AER_FILENAME' : "'input/miefiles/%s'" % aerosol_names[band_choice],
                    'noph' : noph_siro,
                    'ratm' : 6371.0 + 70.0,
                    'step' : 1.0}
                arsca.simu.create_siro_settings(siro_custom_settings)
            arsca.simu.change_raysca_settings('aerosol_file',aerosol_names[band_choice])
            arsca.simu.change_raysca_settings('parallelization',False)
            band_names = ['o2a','wco2','sco2']
            scatter_names = ['no_sca','rayleigh','rayleigh+mie']
            #scatter_names = ['no_sca']
            geom_names = ['limb','nadir','glint']
            
            try:
                with open(validation_fname,'rb') as fhandle:
                    radiance_dict = pickle.load(fhandle)
            except FileNotFoundError:
                print("No file found. Setting up a new one...")
                radiance_dict = {}
                for gname in geom_names:
                    radiance_dict[gname] = {}
                    for bname in band_names:
                        radiance_dict[gname][bname] = {}
                        for sname in scatter_names:
                            radiance_dict[gname][bname][sname] = np.array([np.nan])
            
            # this -v is to enable recomputation of mie scattering things
            #if scatterers != 2 and not np.isnan(radiance_dict[geom_choice][band_names[band_choice]][scatter_names[scatterers]][0]).all():
            #if not np.isnan(radiance_dict[geom_choice][band_names[band_choice]][scatter_names[scatterers]][0]).all():
            if False:
                #raise ValueError("This configuration already computed!") # this is to "return" or halt
                print("This configuration already computed!")
                continue
            
            
            band_lims = [[755.0,775.0],[1590.0,1620.0],[2040.0, 2080.0]]
            #band_albs = [0.9,0.1,0.02] # solar zenith 55 from the picture
            #band_albs = [0.87651711, 0.17815365, 0.05095223] # from the model
            band_albs = [0.88, 0.18, 0.05] # from the model
            band_albs = [0.8, 0.8, 0.8] # from the model
            band_albs = [0.15, 0.15, 0.15]
            #band_albs = [0.0, 0.0, 0.0]
            arsca.set_case("raysca_validation_study_newspec_%s" % band_names[band_choice])
            arsca.set_configuration("validation_study_%s" % band_names[band_choice])
            
            wl_range = band_lims[band_choice] # nm
            
            
            wn_step = 0.01
            wn_range = arsca.tf.wl2wn(wl_range)
            wn_range[0] = wn_range[0] - wn_step #widen it a bit
            wn_range[1] = wn_range[1] + wn_step
            n_wn = np.array(np.arange(wn_range[0],wn_range[1],wn_step)).size
            #NOTE: Here we need to reconfigure the spectra computations for
            # larger spectra.
            
            n_wl = n_wn
            #n_wl_coarse = 200
            n_wl_coarse = n_wl
            wl = np.linspace(wl_range[0],wl_range[1],n_wl_coarse)
            wn_coarse = arsca.tf.wl2wn(wl)
            wn_reso = np.max(np.abs(np.diff(wn_coarse)))
            
            R_earth = 6371.0 # km
            altitudes = np.linspace(0,70,100)
            
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
            n_scatterer = scatterers if scatterers > 0 else 1
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
            
            #scatterers and their cross-sections
            medium['scatterer'] = np.zeros((n_medium_positions,n_scatterer))
            medium['scattering_cross_section'] = np.zeros((n_medium_positions,n_wl,n_scatterer))
            xsec_sca = np.zeros((n_medium_positions,n_wl,n_scatterer))
                 
            if scatterers == 0: # no scattering
                #for proper interaction with siro this probably needs to be 
                #so that the scatterer amount is one, but cross-section is zero.
                sca_interpfun = interp1d(atmos['Height'],atmos['Density'],fill_value="extrapolate")
                medium['scatterer'][:,0] = sca_interpfun(altitudes).reshape((altitudes.size,))
                # With zero multiplication, so that the scattering doesn't happen!
                rayleigh_xsec = 0 * arsca.scatter.rayleigh_xsec(wl,extrapolate=True)
                xsec_sca[:,:,0] = np.repeat(rayleigh_xsec.reshape((1,n_wl)),n_medium_positions,axis=0)
            
            elif scatterers in [1, 2]: # rayleigh stuffs
                sca_interpfun = interp1d(atmos['Height'],atmos['Density'])
                medium['scatterer'][:,0] = sca_interpfun(altitudes).reshape((altitudes.size,))
                rayleigh_xsec = arsca.scatter.rayleigh_xsec(wl,extrapolate=True)
                xsec_sca[:,:,0] = np.repeat(rayleigh_xsec.reshape((1,n_wl)),n_medium_positions,axis=0)
            
            if scatterers == 2: # mie aerosols
                # it is constant through the whole band; from the original mie file
                max_conc = 500
                if geom_choice == 'limb':
                    mean_conc_alt = 15
                else:
                    mean_conc_alt = 2.5
                std_conc_alt = 2.5
                
                aerosols = max_conc * np.exp(-((altitudes - mean_conc_alt) / std_conc_alt) ** 2 )
                aerosol_xsec = np.array([0.017075, 0.00084129, 0.00030242]) # micron^2
                aerosol_xsec = aerosol_xsec * (1e-6 / 1e-2) ** 2 #cm^2
                aerosol_xsec_band = aerosol_xsec[band_choice] * np.ones((1,n_wl))
                xsec_sca[:,:,1] = np.repeat(aerosol_xsec_band.reshape((1,n_wl)),n_medium_positions,axis=0)
                current_aod_765 = np.sum(np.diff(altitudes) * aerosols[:-1] * aerosol_xsec[0] * 1e5)
                aod_765 = 0.05 # at nadir
                aerosols = aerosols * aod_765 / current_aod_765
                medium['scatterer'][:,1] = aerosols
            
            medium['scattering_cross_section'] = xsec_sca
            
            # zero stands for rayleigh-scattering
            medium['scatterer_kernel'] = np.zeros((n_scatterer,))
            if scatterers == 2:
                medium['scatterer_kernel'][1] = 1 # this is to use the aerosols
            medium['scatterer_kernel_parameter'] = np.zeros((1,n_scatterer))
            
            #generation of the absorbtion cross-sections
            
            medium['absorber'] = np.zeros((n_medium_positions,n_absorber))
            medium['absorbing_cross_section'] = np.zeros((n_medium_positions,n_wl,n_absorber))
            
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
                                                    gases_hapi_id[abs_idx],temp,pres,p_s)
                
                #xsec_conv = np.zeros_like(xsec)
                #for i in range(xsec.shape[0]):
                #    wn_conv,xsec_conv[i,:],_,_,_ = arsca.xsec.hapi.convolveSpectrumSame(wn,xsec[i,:],Resolution=wn_reso,SlitFunction=arsca.xsec.hapi.SLIT_GAUSSIAN)
                xsec_conv = xsec
                
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
            # three geomtries: limb, nadir, glint
            instrument = {}
            sat_altitude = 705.0
            
            n_los = 1
            instrument['position'] = np.zeros((n_los,3))
            instrument['view_vector'] = np.zeros((n_los,3))
            
            z_ax = np.array([0.0, 0.0, 1.0])
            y_ax = np.array([0.0, 1.0, 0.0])
            x_ax = np.array([1.0, 0.0, 0.0])
            nadir = np.array([-1.0,0.0,0.0])
            
            if geom_choice == 'limb':
                tan_altitude = 15
                solar_zen_degs = solar_zens[0]
                solar_azi_degs = solar_azis[0]
                (sat_pos,sat_view_vec,_,solar_dir,_) = arsca.tf.limb_geometries(sat_altitude,np.array([tan_altitude]),np.array([solar_zen_degs]),np.array([solar_azi_degs]))
                instrument['position'][0,:] = sat_pos
                instrument['view_vector'][0,:] = sat_view_vec
            elif geom_choice == 'nadir':
                solar_zen_degs = solar_zens[1]
                solar_azi_degs = solar_azis[1]
                solar_dir = arsca.tf.solar_direction(x_ax,y_ax,solar_zen_degs,solar_azi_degs)
                sat_position = np.array([R_earth + sat_altitude, 0.0, 0.0])
                instrument['position'][0,:] = sat_position
                instrument['view_vector'][0,:] = nadir
            elif geom_choice == 'glint':
                solar_zen_degs = solar_zens[2]
                solar_azi_degs = solar_azis[2]
                solar_dir = arsca.tf.solar_direction(x_ax,y_ax,solar_zen_degs,solar_azi_degs)
                sat_glint_dir = arsca.tf.arb_rotation(x_ax,-solar_zen_degs * np.pi/180,z_ax)
                instrument['position'][0,:] = sat_altitude * sat_glint_dir + np.array([R_earth,0.0,0.0])
                # in this case, the satellite is just at 428 km altitude
                instrument['view_vector'][0,:] = -sat_glint_dir
                    
            
            # BOUNDARY DEFINITIONS
            boundary = {}
            boundary['shape'] = np.array([1,1]) #both are spherical surfaces
            boundary['parameter'] = np.array([6371.0,6371.0 + 70]) #the radius of the spheres
            
            #reflections from the boundary (0=pass through,1=lambertian,2=semispecular,3=brdf)
            # TODO: Hide this to the background
            #boundary['reflection_kernel'] = np.array([3,0])
            boundary['reflection_kernel'] = np.array([1,0])
            boundary['reflection_kernel_parameter'] = np.array([[band_albs[band_choice],0.0]])
            
            source = {}
            
            #Just the glint source
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
            
            
            raysca_start = time.time()
            radiance_raysca = arsca.simu.run('raysca',input_fname)
            #radiance_raysca = [0]
            
            raysca_time = time.time() - raysca_start
            if compute_siro:
                siro_start = time.time()
                radiance_siro = arsca.simu.run('siro',input_fname)
                siro_time = time.time() - siro_start
            #    rads_siro.append(radiance_siro)
            #else:
            #    rads_siro.append(np.nan)
            
            
    
            #rads_raysca.append(radiance_raysca)
            print(str(scatterers) + geom_choice)
            
            print(raysca_time)
            if not compute_siro:
                radiance_siro = radiance_dict[geom_choice][band_names[band_choice]][scatter_names[scatterers]][1]
            else:
                print(siro_time)
            radiance_dict[geom_choice][band_names[band_choice]][scatter_names[scatterers]] = [radiance_raysca, radiance_siro]
            
            with open(validation_fname,'wb') as fhandle:
                # Writing results to the file.
                pickle.dump(radiance_dict,fhandle)
