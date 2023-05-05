#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 02:11:56 2020

@author: mikkonea
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

img_folder_name = "./images/validation_study_zero_alb_many_photons/"
img_folder_name = "./images/validation_study_final_v2/"

validation_fname = 'validation_study_radiances.pick'
validation_fname = 'validation_study_radiances_v11.pick'
validation_fname = 'validation_study_radiances_low_alt_aero.pick'
validation_fname = "validation_study_radiances_low_alt_aero_high_alb.pick"
validation_fname = 'validation_study_radiances_low_alt_aero_zero_alb.pick'
validation_fname = 'validation_study_radiances_low_alt_aero_zero_alb_many_photons.pick'
validation_fname = 'validation_study_radiances_final_v2.pick'

def combine_validations():
    newfname = 'validation_study_radiances_comb.pick'
    validation_fname_raysca = 'validation_study_radiances_v11.pick'
    validation_fname_siro = 'validation_study_radiances_oops_all_siros.pick'
    with open(validation_fname_raysca,'rb') as fhandle:
        radiance_dict_raysca = pickle.load(fhandle)
    with open(validation_fname_siro,'rb') as fhandle:
        radiance_dict_siro = pickle.load(fhandle)
    radiance_dict = {}
    for k1 in radiance_dict_raysca.keys():
        radiance_dict[k1] = {}
        for k2 in radiance_dict_raysca[k1].keys():
            radiance_dict[k1][k2] = {}
            for k3 in radiance_dict_raysca[k1][k2].keys():
                radiance_dict[k1][k2][k3] = [radiance_dict_raysca[k1][k2][k3][0],radiance_dict_siro[k1][k2][k3][1]]
                
    with open(newfname,'wb') as fhandle:
        # Writing results to the file.
        pickle.dump(radiance_dict,fhandle)
    #validation_fname = 'validation_study_radiances_oops_all_siros.pick'

with open(validation_fname,'rb') as fhandle:
    radiance_dict = pickle.load(fhandle)

sp = np.genfromtxt('sun_fixed.spec')
#pols = ['H','V']
from scipy.interpolate import interp1d
ifun = interp1d(sp[:,0],sp[:,1])
band_lims = [[755.0,775.0],[1590.0,1620.0],[2040.0, 2080.0]]
band_names = ['O2A','WCO2','SCO2']
#band_names = ['NIR','SWIR-1','SWIR-2']
band_keys = ['o2a','wco2','sco2']

scatter_names = ['no_sca','rayleigh','rayleigh+mie']
#scatter_names = ['no_sca']
scattering_title = ['no scattering', 'Rayleigh scattering', 'Rayleigh + Mie']
mode_names = ['limb','nadir', 'glint']
#mode_names = ['nadir']
title = []

siro_line = ['b-','b:']
#raysca_line = ['r-','r:']
raysca_line = ['r-','r:']
siro_ss_line = ['k-','k:']

# 27 spektrikuvaa ja suhteellista virhettä
# I ja Q erikseen
# RaySca, Siro full, Siro SS

errlistI = []
errlistQ = []
errlistT = []
#band_plot_lims = {'nadir':([-12,240],[][])} 

# Oisko 1-to-1?
for b_idx, bkey in enumerate(band_keys):
    #n_wl = 200
    try:
        n_wl = radiance_dict['limb'][bkey]['rayleigh+mie'][0][0,0,:,0].size
    except IndexError:
        continue
    for s_idx, skey in enumerate(scatter_names):
        for gkey in mode_names:
            rads = radiance_dict[gkey][bkey][skey]
            try:
                I_raysca = rads[0][0,0,:,0]
            except IndexError:
                # no computation done yet
                continue
            Q_raysca = rads[0][0,0,:,1]
            I_siro = np.sum(np.nan_to_num(rads[1][:,0,:,0]),axis=1)
            Q_siro = np.sum(np.nan_to_num(rads[1][:,0,:,1]),axis=1)
            I_siro_ss = np.sum(np.nan_to_num(rads[1][:,0,:2,0]),axis=1)
            Q_siro_ss = np.sum(np.nan_to_num(rads[1][:,0,:2,1]),axis=1)
            
            #plt.figure(figsize=(8,6))
            #f, (axI, axQ, axErr) = plt.subplots(1, 3, sharey=True, figsize=(8,6))
            f, temptup = plt.subplots(nrows=2, ncols=2,sharex=True, figsize=(9,5), constrained_layout=True, dpi=300)
            (axI, axQ) = temptup[0]
            (axErrI, axErrQ) = temptup[1]
            wl_range = band_lims[b_idx] # nm
            wl = np.linspace(wl_range[0],wl_range[1],n_wl)
            radi = ifun(wl)
            
            axI.plot(wl,radi * I_raysca,raysca_line[0],label='RaySca (I)',alpha=0.8)
            #axI.plot(wl,radi * I_siro_ss,siro_ss_line[0],label='Siro single-scattering (I)',alpha=0.5)
            #axI.set_ylim([150,180])
            axI.plot(wl,radi * I_siro,siro_line[1],label='Siro (I)',alpha=0.5)
            axQ.plot(wl,radi * Q_raysca,raysca_line[0],label='RaySca (Q)',alpha=0.8)
            #axQ.plot(wl,radi * Q_siro_ss,siro_ss_line[0],label='Siro single-scattering (Q)',alpha=0.5)
            axQ.plot(wl,radi * Q_siro,siro_line[1],label='Siro (Q)',alpha=0.5)
            #axErr.plot(wl, 100 * (I_raysca - I_siro_ss) / I_siro_ss, siro_line[0],label='Single-scattering (I)')
            #axErr.plot(wl, 100 * (I_raysca - I_siro) / I_siro, siro_line[1],label='Full radiance (I)')
            #axErr.plot(wl, 100 * (Q_raysca - Q_siro_ss) / Q_siro_ss, raysca_line[0],label='Single-scattering (Q)')
            #axErr.plot(wl, 100 * (Q_raysca - Q_siro) / Q_siro_ss, raysca_line[1],label='Full radiance (Q)')
            axErrI.plot(wl,radi *  (I_raysca - I_siro_ss), raysca_line[0],label='RaySca - Single-scattering Siro (I)',alpha=0.5)
            axErrI.plot(wl,radi *  (I_raysca - I_siro), siro_line[1],label='RaySca - Siro (I)',alpha=0.5)
            axErrQ.plot(wl,radi *  (Q_raysca - Q_siro_ss), raysca_line[0],label='RaySca - Single-scattering Siro (Q)',alpha=0.5)
            axErrQ.plot(wl,radi *  (Q_raysca - Q_siro), siro_line[1],label='RaySca - Siro (Q)',alpha=0.5)
            
            axErrQ.set_xlabel('wavelength (nm)')
            axErrI.set_xlabel('wavelength (nm)')
            axI.set_ylabel('radiance (mW/m^2/sr/nm)')
            axErrI.set_ylabel('radiance difference')
            f.suptitle('%s band, %s mode, %s' % (band_names[b_idx], gkey, scattering_title[s_idx]))
            axI.set_title('Stokes I')
            axQ.set_title('Stokes Q')
            #axErr.set_title('Relative error (%)')
            axErrI.set_title('Absolute error (Raysca - Siro) in I')
            axErrQ.set_title('Absolute error (Raysca - Siro) in Q')
            #axErr.set_ylim([-1.0,1.0])
            #axI.legend()
            #axQ.legend()
            #axErr.legend()
            img_name = "%s_%s_%s.jpg" % (gkey,bkey,skey)
            plt.savefig(img_folder_name + img_name)
            
            #plt.show()
            
            err_eps = 1e-6
            err_maskI = np.abs(I_siro) > err_eps
            err_maskQ = np.abs(Q_siro) > err_eps
            err_maskIss = np.abs(I_siro_ss) > err_eps
            err_maskQss = np.abs(Q_siro_ss) > err_eps
            errspec_I = ( I_raysca[err_maskI] - I_siro[err_maskI]) / I_siro[err_maskI]
            errspec_Q = ( Q_raysca[err_maskQ] - Q_siro[err_maskQ]) / Q_siro[err_maskQ]
            errspec_Iss = ( I_raysca[err_maskIss] - I_siro_ss[err_maskIss]) / I_siro_ss[err_maskIss]
            errspec_Qss = ( Q_raysca[err_maskQss] - Q_siro_ss[err_maskQss]) / Q_siro_ss[err_maskQss]
            
            if skey == 'no_sca':
                continue
                #we'll skip all the non-scattering segments
            
            def find_outlier(spec):
                if spec.size == 0:
                    return 0
                abs_errs = np.abs(np.array([np.min(spec),np.max(spec)]))
                if abs_errs[0] > abs_errs[1]:
                    return np.min(spec)
                else:
                    return np.max(spec)
            
            bias_I = np.mean(errspec_I)
            bias_Q = np.mean(errspec_Q)
            bias_Iss = np.mean(errspec_Iss)
            bias_Qss = np.mean(errspec_Qss)
            errstd_I = np.abs(np.std(errspec_I))
            errstd_Q = np.abs(np.std(errspec_Q))
            errstd_Iss = np.abs(np.std(errspec_Iss))
            errstd_Qss = np.abs(np.std(errspec_Qss))
            errmax_I = find_outlier(errspec_I)
            errmax_Q = find_outlier(errspec_Q)
            errmax_Iss = find_outlier(errspec_Iss)
            errmax_Qss = find_outlier(errspec_Qss)
            
            print_bands = {'o2a' : 'O$_2$ A ', 'wco2': 'Weak CO$_2$', 'sco2' : 'Strong CO$_2$'}
            print_scatt = {'no_sca' : 'no scattering', 'rayleigh': 'Rayleigh', 'rayleigh+mie' : 'R. \\& Aerosol'}
            print_geom = {'limb' : 'Limb', 'nadir': 'Nadir', 'glint' : 'Glint'}
            
            decimf = '1.1f'

            err_string = "%s & %s & %s" % (print_bands[bkey], print_scatt[skey], print_geom[gkey])
            
            #err_stringI = err_string + f"& %{decimf} (%{decimf}) & %{decimf} (%{decimf}) & %{decimf} (%{decimf}) \\\\" % (100 * bias_Iss, 100 * bias_I, 100 * errstd_Iss, 100 * errstd_I, 100 * errmax_Iss, 100 * errmax_I)
            #err_stringQ = err_string + f"& %{decimf} (%{decimf}) & %{decimf} (%{decimf}) & %{decimf} (%{decimf}) \\\\" % (100 * bias_Qss, 100 * bias_Q, 100 * errstd_Qss, 100 * errstd_Q, 100 * errmax_Qss, 100 * errmax_Q)
            err_string_total = err_string + f"& %{decimf} (%{decimf}) & %{decimf} (%{decimf}) & %{decimf} (%{decimf}) & %{decimf} (%{decimf}) \\\\" % (100 * bias_Iss, 100 * bias_I, 100 * errstd_Iss, 100 * errstd_I, 100 * bias_Qss, 100 * bias_Q, 100 * errstd_Qss, 100 * errstd_Q)
            
            #errlistI.append(err_stringI)
            #errlistQ.append(err_stringQ)
            errlistT.append(err_string_total)
            

scatter_names = ['no_sca','rayleigh','rayleigh+mie']
mode_color = ['b','r','y']

binlims1 = {'o2a' : [[-20,5],[-10,10]],
           'wco2' : [[-15,4],[-7.5,10]],
           'sco2' : [[-15,5],[-7.5,10]]}

binlims2 = {'o2a' : [[-40,8],[-35,15]],
           'wco2' : [[-20,7],[-10,15]],
           'sco2' : [[-17,6],[-10,15]]}

binlims = {'rayleigh' : binlims1,
         'rayleigh+mie' : binlims2}
for b_idx, bkey in enumerate(band_keys):
    
    try:
        n_wl = radiance_dict['limb'][bkey]['rayleigh+mie'][0][0,0,:,0].size
    except IndexError:
        continue
    
    for s_idx, skey in enumerate(scatter_names):
        f, temptup = plt.subplots(nrows=2, ncols=2, figsize=(9,5), constrained_layout=True, dpi=300)
        (axIss, axQss) = temptup[0]
        (axI, axQ) = temptup[1]
        axQ.set_title('Q, multiple scattering')
        axI.set_title('I, multiple scattering')
        axQss.set_title('Q, single scattering')
        axIss.set_title('I, single scattering')
            
        for g_idx,gkey in enumerate(mode_names):
            rads = radiance_dict[gkey][bkey][skey]
            try:
                I_raysca = rads[0][0,0,:,0]
            except IndexError:
                # no computation done yet
                continue
            Q_raysca = rads[0][0,0,:,1]
            I_siro = np.sum(np.nan_to_num(rads[1][:,0,:,0]),axis=1)
            Q_siro = np.sum(np.nan_to_num(rads[1][:,0,:,1]),axis=1)
            I_siro_ss = np.sum(np.nan_to_num(rads[1][:,0,:2,0]),axis=1)
            Q_siro_ss = np.sum(np.nan_to_num(rads[1][:,0,:2,1]),axis=1)
            
            err_eps = 1e-6
            err_maskI = np.abs(I_siro) > err_eps
            err_maskQ = np.abs(Q_siro) > err_eps
            err_maskIss = np.abs(I_siro_ss) > err_eps
            err_maskQss = np.abs(Q_siro_ss) > err_eps
            errspec_I = ( I_raysca[err_maskI] - I_siro[err_maskI]) / I_siro[err_maskI]
            errspec_Q = ( Q_raysca[err_maskQ] - Q_siro[err_maskQ]) / Q_siro[err_maskQ]
            errspec_Iss = ( I_raysca[err_maskIss] - I_siro_ss[err_maskIss]) / I_siro_ss[err_maskIss]
            errspec_Qss = ( Q_raysca[err_maskQss] - Q_siro_ss[err_maskQss]) / Q_siro_ss[err_maskQss]
            # Note: First time escaping %-sign in %-formatted string. Why is it %% and not \%?!
            binnum = 100
            f.suptitle('Relative error (%%) (Raysca - Siro) / Siro \n  %d bins, %d wavelengths, %s band, %s' % (binnum,I_raysca.size,band_names[b_idx], scattering_title[s_idx]))
            binlim = binlims[skey][bkey]
            binsI = np.linspace(binlim[0][0],binlim[0][1],binnum+1)
            binsQ = np.linspace(binlim[1][0],binlim[1][1],binnum+1)
            axI.hist(100*errspec_I,binsI,alpha=0.6,color=mode_color[g_idx])
            axQ.hist(100*errspec_Q,binsQ,alpha=0.6,color=mode_color[g_idx])
            axIss.hist(100*errspec_Iss,binsI,alpha=0.6,label=gkey,color=mode_color[g_idx])
            axQss.hist(100*errspec_Qss,binsQ,alpha=0.6,color=mode_color[g_idx])
        img_name = "errorhist_%s_%s.jpg" % (bkey,skey)
        axI.set_xlabel("Relative error (%)")
        axQ.set_xlabel("Relative error (%)")
        axIss.set_ylabel('# of wavelengths')
        axI.set_ylabel('# of wavelengths')
        axIss.legend()
        plt.savefig(img_folder_name + img_name)
            
        #f.show()

            
for e in errlistT:
    print(e)
for e in errlistI:
    print(e)
for e in errlistQ:
    print(e)