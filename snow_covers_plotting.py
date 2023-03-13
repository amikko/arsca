#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 02:11:56 2020

@author: mikkonea
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

img_folder_name = "./images/snow_surf_comparison_v5/"
#img_folder_name = "./images/water_surf_simu/"

validation_fname = 'snow_surf_comparison.pick'
validation_fname = 'snow_surf_comparison_new_art.pick'
validation_fname = 'snow_surf_v5.pick'
with open(validation_fname,'rb') as fhandle:
    radiance_dict = pickle.load(fhandle)

sp = np.genfromtxt('sun_fixed.spec')
#pols = ['H','V']
from scipy.interpolate import interp1d
ifun = interp1d(sp[:,0],sp[:,1])
band_lims = [[755.0,775.0],[1590.0,1620.0],[2040.0, 2080.0]]
band_names = ['O2A','10 x WCO2','100 x SCO2']
band_keys = ['o2a','wco2','sco2']

surf_names = ['art','loose_snow','new_snow','rough_snow','lambertian']
surf_names_ = ['ART', 'Loose snow', 'New snow', 'Rough snow','Lambertian']

geom_names = surf_names

title = []

siro_line = ['b-','b:']
#raysca_line = ['r-','r:']
#raysca_line = ['r-','r:']
raysca_line = ['rx-','bx-','kx-']

siro_ss_line = ['k:','k:']
los_angs = np.linspace(0,85,18)
#ang_choices = [0,1,2,3,4,5,6,7,8]
#los_angs = np.array([0,30,70])
#ang_choices = [0,1,2]

zen_ang_choices = [[0,1,2],[2],[0,1],[1],[0,1,2]]

mult_coeffs= [1,10,100]

# summakuvat eri zeniiteille ja geometrioille: 12 kpl
# kaikki bandit samaan kuvaan? Yes


solar_zens = [36,54,69]
errlistI = []
errlistQ = []

siro_plotting = False
plot_spectra = True

for g_idx,gkey in enumerate(geom_names):
    #for zen_i in range(len(solar_zens)):
    for zen_i in zen_ang_choices[g_idx]:
        f, (axI, axQ) = plt.subplots(1, 2, figsize=(14,4))
        max_I = -1
        min_I = 1
        max_Q = -1
        min_Q = 1
        for b_idx, bkey in enumerate(band_keys):
            n_wl = radiance_dict['lambertian'][bkey][0][0,0,:,0].size
            wl_range = band_lims[b_idx] # nm
            wl = np.linspace(wl_range[0],wl_range[1],n_wl)
            radi = ifun(wl)
            rads = radiance_dict[gkey][bkey]
            I_raysca = np.mean(mult_coeffs[b_idx] * radi * rads[0][zen_i,:,:,0],axis=1)
            Q_raysca = np.mean(mult_coeffs[b_idx] * (radi * rads[0][zen_i,:,:,1]),axis=1)
            axI.plot(los_angs,I_raysca,raysca_line[b_idx],label='%s' % band_names[b_idx])
            axQ.plot(los_angs,Q_raysca,raysca_line[b_idx],label='%s' % band_names[b_idx])
            max_I = np.max(I_raysca) if np.max(I_raysca) > max_I else max_I
            min_I = np.min(I_raysca) if np.min(I_raysca) < min_I else min_I
            max_Q = np.max(Q_raysca) if np.max(Q_raysca) > max_Q else max_Q
            min_Q = np.min(Q_raysca) if np.min(Q_raysca) < min_Q else min_Q
        axI.plot([solar_zens[zen_i],solar_zens[zen_i]],[min_I,max_I],'k:')
        axQ.plot([solar_zens[zen_i],solar_zens[zen_i]],[min_Q,max_Q],'k:')
        axQ.set_xlabel('Viewing zenith angle (degrees)')
        axI.set_xlabel('Viewing zenith angle (degrees)')
        axI.set_ylabel('Average radiance (mW/m^2/sr/nm)')
        axQ.set_ylabel('Average radiance (mW/m^2/sr/nm)')
        f.suptitle('%s, SZA %d degrees' % (surf_names_[g_idx], solar_zens[zen_i]))
        axI.set_title('Stokes I')
        axQ.set_title('Stokes Q')
        axI.legend()
        axQ.legend()
        #plt.show()
                
        img_name = "angs_%s_%s_%d.png" % (gkey,bkey,solar_zens[zen_i])
        plt.savefig(img_folder_name + img_name)


raysca_line = ['r-','k-','b-']
band_names = ['O2A','WCO2','SCO2']
#indices are 4, 6 and 7
ang_indices = [7,11,14]
if plot_spectra:
    for b_idx, bkey in enumerate(band_keys):
        n_wl = radiance_dict['lambertian'][bkey][0][0,0,:,0].size
        for g_idx,gkey in enumerate(geom_names):
            for zen_i in zen_ang_choices[g_idx]:
                ang_choices = [0,ang_indices[zen_i]]
                f, (axI, axQ) = plt.subplots(1, 2, figsize=(14,4))
                for choice_idx,ang_idx in enumerate(ang_choices):
                    i = ang_idx
                    rads = radiance_dict[gkey][bkey]
                    I_raysca = rads[0][zen_i,i,:,0]
                    Q_raysca = rads[0][zen_i,i,:,1]
                    if siro_plotting:
                        I_siro = np.sum(np.nan_to_num(rads[1][:,i,:,0]),axis=1)
                        Q_siro = np.sum(np.nan_to_num(rads[1][:,i,:,1]),axis=1)
                        I_siro_ss = np.sum(np.nan_to_num(rads[1][:,i,:2,0]),axis=1)
                        Q_siro_ss = np.sum(np.nan_to_num(rads[1][:,i,:2,1]),axis=1)
                    #
                    #plt.figure(figsize=(8,6))
                    #f, (axI, axQ, axErr) = plt.subplots(1, 3, sharey=True, figsize=(8,6))
                    wl_range = band_lims[b_idx] # nm
                    wl = np.linspace(wl_range[0],wl_range[1],n_wl)
                    radi = ifun(wl)
                    angstr = "%1.1f" % los_angs[ang_idx]
                    axI.plot(wl,radi * I_raysca,raysca_line[choice_idx % 3],label='%s' % angstr)
                    #axI.plot(wl,radi * I_siro_ss,siro_line[1],label='Siro (I)')
                    axQ.plot(wl,radi * Q_raysca,raysca_line[choice_idx % 3],label='%s' % angstr)
                    #axQ.plot(wl,radi * Q_siro_ss,siro_line[1],label='Siro (Q)')
                    
                    plt.plot()
                #axErr.plot(wl, 100 * (I_raysca - I_siro_ss) / I_siro_ss, siro_line[0],label='Single-scattering (I)')
                #axErr.plot(wl, 100 * (I_raysca - I_siro) / I_siro, siro_line[1],label='Full radiance (I)')
                #axErr.plot(wl, 100 * (Q_raysca - Q_siro_ss) / Q_siro_ss, raysca_line[0],label='Single-scattering (Q)')
                #axErr.plot(wl, 100 * (Q_raysca - Q_siro) / Q_siro_ss, raysca_line[1],label='Full radiance (Q)')
                #axErr.plot(wl,radi *  (I_raysca - I_siro_ss), siro_line[0],label='RaySca - Single-scattering Siro (I)')
                #axErr.plot(wl,radi *  (I_raysca - I_siro), siro_line[1],label='RaySca - Siro (I)')
                #axErr.plot(wl,radi *  (Q_raysca - Q_siro_ss), raysca_line[0],label='RaySca - Single-scattering Siro (Q)')
                #axErr.plot(wl,radi *  (Q_raysca - Q_siro), raysca_line[1],label='RaySca - Siro (Q)')
                
                axQ.set_xlabel('wavelength (nm)')
                axI.set_ylabel('radiance (mW/m^2/sr/nm)')
                f.suptitle('%s, %s band, %d degrees' % (surf_names_[g_idx],band_names[b_idx], solar_zens[zen_i]))
                axI.set_title('Stokes I')
                axQ.set_title('Stokes Q')
                #axErr.set_title('Relative error (%)')
                #axErr.set_title('Absolute error (RaySca - Siro)')
                #axErr.set_ylim([-1.0,1.0])
                axI.legend()
                axQ.legend()
                #axErr.legend()
                img_name = "%s_%s_%d.png" % (gkey,bkey,solar_zens[zen_i])
                plt.savefig(img_folder_name + img_name)
                #plt.show()
                
                if siro_plotting:
                    err_eps = 1e-6
                    err_maskI = np.abs(I_siro) > err_eps
                    err_maskQ = np.abs(Q_siro) > err_eps
                    #err_maskIss = np.abs(I_siro_ss) > err_eps
                    #err_maskQss = np.abs(Q_siro_ss) > err_eps
                    errspec_I = ( I_raysca[err_maskI] - I_siro[err_maskI]) / I_siro[err_maskI]
                    errspec_Q = ( Q_raysca[err_maskQ] - Q_siro[err_maskQ]) / Q_siro[err_maskQ]
                #errspec_Iss = ( I_raysca[err_maskIss] - I_siro_ss[err_maskIss]) / I_siro_ss[err_maskIss]
                #errspec_Qss = ( Q_raysca[err_maskQss] - Q_siro_ss[err_maskQss]) / Q_siro_ss[err_maskQss]
