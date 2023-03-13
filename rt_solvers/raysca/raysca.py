#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 03:57:36 2019

@author: taintur
"""

"""
23.11.2022
Mikkonen, Antti

Setting up the full 3D atmosphere can be a bit cumbersome at this point.
Maybe we'll do the research by creating some kind of geometry which we'll vary
iteratively while doing the inversion!

Could there be an optimization based on the medium points we're tracing through?
In the neearest neighbour setup our steps would be simplified considerably in
some larger structures?
"""

"""
11.3.2022
Mikkonen, Antti
Okay, now we need to make some dazzling optimizations and fast. Plan is as follows:

    The path tracing is done first and sequentially.
        The path contains the nodes on how radiance changes directions
        For each scatterer a new path is created
        path segment transmittances are computed beforehand
        path transmittances are then combined from scattering phase matrices and
        path segment transmittances

    The transmissivity of each of the paths is computed afterwards using fancy
parallel methods. The parallelization comes into account in the integration part.
        The transmissivity computations initially check which interpolations
        are needed along each of the paths, depending on the integration length.
        DONE! extinction is precomputed!
        The large path systems are


"""
def keyboard(banner=None):
    import sys
    import code
    ''' Function that mimics the matlab keyboard command '''
    # use exception trick to pick up the current frame
    # os.system("""osascript -e 'display notification "Keyboard ready" with title "Ready"'""")
    try:
        raise None
    except Exception:
        frame = sys.exc_info()[2].tb_frame.f_back
    print("# Use quit() to exit :) Happy debugging!")
    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner, local=namespace)
    except SystemExit:
        return

# TODO: This needs to be refactored into smaller modules
# TODO: Before modifying stuff too much, compute a test case against which it
# is easy to test if things work correctly.

# TODO: Verify that surface reflectances are handled in correct order, i.e.
# inc_zen is in the direction of the source. This shouldn't be an issue if the
# brdf is symmetrical with respect to the zenith angles.

# TODO: Wavelength-dependent muller matrices for scattering?

# TODO: ARSCA-function to change RaySca settings for automating settings changes

# TODO: If multiple-scattering with Monte Carlo is implemented, then the
# polaris-function from Siro is needed

# TODO: Should scattering and absorbing elements be separated?

# TODO: Non-far-field sources 1/r radiance diminishing information can be stored
# in the norm of the source_direction function.

# TODO: To make this more accessible and nicer, the lambda interpolators should
# be implemented using classes! This would eliminate the need for 'pathos'
# package. This should be done in the refactoring phase

# TODO: The parallelization over wavelengths should be smart, i.e. so the maximum
# number of cores is utilized

# TODO: What happens, when a beam is transmitted and reflected at the same time?
# Happens also with refraction.

# TODO: Add logic to trace to multiple sources at the same time.

# TODO: Optimize the computation by precomputing the cross-section and medium
# products as they're used everywhere. This will be extremely useful when computing
# camera-like images

import netCDF4

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# This will just get the processor counts
import multiprocessing

# yaml is for reading the settings
import yaml

# tqdm enables the progress bars
import tqdm
from .reflection import set_up_reflection
from .memory import total_size
# These matrices are here to transform stokes between two basises. The input
# and the output of raysca is in IQUV, but the internals are in IIUV, since
# they're copied from Siro.
IIUV2IQUV = np.array([[1.0, 1.0, 0.0, 0.0],
                      [1.0,-1.0, 0.0, 0.0],
                      [0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0]])

IQUV2IIUV = np.linalg.inv(IIUV2IQUV)

cm_in_km = 1e5

settings_file = "raysca_settings.yaml"

#phase folder contains all the phase matrices used
phase_folder = "phase_matrices/"

paths_set = False

def set_paths(raysca_path):
    global settings_file, phase_folder, paths_set
    if not paths_set:
        basepath = raysca_path[:-len('raysca.py')]
        settings_file = basepath + settings_file
        phase_folder = basepath + phase_folder
        paths_set = True

def read_settings():
    # These are the parameters loaded from the configuration files
    global main_beam_step_length, scattering_step_length
    global max_distance, single_scattering, extinction_threshold, monte_carlo_mode
    global single_trace_for_all_wl, parallelization, process_count, show_progress
    global source_los_cartesian_product, aerosol_file, brdf_file, optimized_mode
    global stokes_output_format, super_accurate, interpolator
    global gaussian_drop_off, gaussian_sum

    with open(settings_file, 'r') as stream:
        try:
            settings = yaml.safe_load(stream)
            main_beam_step_length = float(settings['main_beam_step_length'])
            scattering_step_length = float(settings['scattering_step_length'])
            max_distance = float(settings['max_distance'])
            single_scattering = bool(settings['single_scattering'])
            monte_carlo_mode = bool(settings['monte_carlo_mode'])
            optimized_mode = bool(settings['optimized_mode'])
            extinction_threshold = float(settings['extinction_threshold'])
            single_trace_for_all_wl = bool(settings['single_trace_for_all_wl'])
            parallelization = bool(settings['parallelization'])
            process_count = int(settings['process_count'])
            show_progress = bool(settings['show_progress'])
            source_los_cartesian_product = bool(settings['source_los_cartesian_product'])
            aerosol_file = str(settings['aerosol_file'])
            brdf_file = str(settings['brdf_file'])
            stokes_output_format = str(settings['stokes_output_format'])
            super_accurate = bool(settings['super_accurate'])
            interpolator = str(settings['interpolator'])
            gaussian_sum = bool(settings['gaussian_sum'])
            if stokes_output_format not in ['IIUV','IQUV']:
                raise ValueError("Invalid Stokes output format type: %s" % stokes_output_format)
        except yaml.YAMLError as exc:
            print(exc)
            raise ValueError("RaySca couldn't read settings!")
    """"
    if parallelization:
        # pathos enables the parallelization of lambda functions (NOTE: This isn't in
        # the usual Anaconda distribution, but it is obtainable using the command)
        # pip install git+https://github.com/uqfoundation/pathos.git@master
        import pathos
    """

def read_input_nc(filename):
    """
    This reads the input filename and populates the data dictionaries.
    """
    groups = ['medium','instrument','source','boundary']
    medium = {}
    instrument = {}
    source = {}
    boundary = {}
    with netCDF4.Dataset(filename) as ds:
        for group in groups:
            variables = list(ds[group].variables.keys())
            for var in variables:
                exec(group + "['%s'] = ds['%s']['%s'][:].data" % (var,group,var))
    return (medium, instrument, source, boundary)

def interpolate_wavelengths(cross_section, wavelengths_in, wavelengths_out):
    # This part will interpolate the input wavelengths into output ones
    # NOTE: This is highly illegal. Use only when testing out numerics!
    # This should be done in ARSCA with corresponding instrument function.
    f = interp1d(wavelengths_in, cross_section, axis=1)
    return f(wavelengths_out)

# Suunnitelma tuolle taichille:
# Luo extinktiofunktiot ti.funktioina
# Kaikkien ekstintioiden akkumulaatio on ti.kernel!


def gauss_fun(pos_in,pos_base,std_base,val_base):
    # 3-dimensional isotropic Gaussian
    # Source: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

    std_max = 4 # Further away than this, the value is zero.
    # Value of 4 corresponds to 2sigma
    var = std_base ** 2
    diff = pos_in - pos_base
    diffTdiff = np.diag(diff @ diff.T)
    if len(val_base.shape) > 1:
        value = np.zeros((val_base.shape[1],))
    else:
        value = np.zeros((1,))

    for i in range(diffTdiff.shape[0]):
        if diffTdiff[i] > std_max * var[i]:
            continue
        elif not gaussian_sum:
            exp = - diffTdiff / var
            idx = np.argmax(exp)
            div = np.sqrt((2 * np.pi) ** 3 * var[idx])
            value = val_base[idx] / div * np.exp(exp[idx])
            break
        else:
            exp = - diffTdiff / var
            div = np.sqrt((2 * np.pi) ** 3 * var)
            if len(val_base.shape) > 1:
                value = np.sum(val_base / div.reshape((div.size,1)) * np.exp(exp.reshape((exp.size,1))),axis=0)
            else:
                value = np.sum(val_base / div * np.exp(exp))
            break
    return value


def set_up_interpolation(medium, wavelengths_in, wavelengths_out):
    """
    Sets up the interpolation functions for the medium positions.
    Also pre-interpolates the cross-section input wavelengths to output
    wavelengths.

    Returns a list of all interpolation functions used in the medium.
    """

    # TODO: Set up classes instead of lambda functions for the function templates
    # Change lambdas to classes in other functions also. The lambdas are the
    # reason why multiprocessing requires an external dependency for serialization
    # NOTE: Functools.partial is supposedly a nicer way to do this.

    n_abs = medium['absorber'].shape[1]
    n_sca = medium['scatterer'].shape[1]
    n_emi = medium['emitter'].shape[1]

    intp_fun_selection = medium['interpolation_function'][:]

    unique_funs = np.unique(intp_fun_selection)

    interp_funs = {'absorber' : [],
                   'scatterer' : [],
                   'emitter' : [],
                   'absorbing_cross_section' : [],
                   'scattering_cross_section' : [],
                   'medium_emissivity' : [],
                   'extinction' : []}

    for current_fun_selection in unique_funs:
        mask = (intp_fun_selection == current_fun_selection)
        if current_fun_selection == 1:
            # spherical shells interpolation
            # This will assume that the spherical elements are in given in order,
            # which you really should do!
            base_interp_fun = lambda pos : np.linalg.norm(pos)
            R = np.array(list(map(np.linalg.norm, medium['position'][mask])))
            # You'd might think that the lambda functions here would be unnecessarily complex, and you'd be right!
            # It appears that if you construct lambda functions in a loop with different variables, then all of them would use the last one
            # Who knows what's going on with this!

            for i_abs in range(n_abs):
                interp_funs['absorber'].append(lambda pos, medium_vals=medium['absorber'][mask,i_abs] : interp1d(R, medium_vals, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
                xsec_abs = interpolate_wavelengths(medium['absorbing_cross_section'][mask,:,i_abs], wavelengths_in, wavelengths_out)
                interp_funs['absorbing_cross_section'].append(lambda pos, medium_vals=xsec_abs : interp1d(R, medium_vals, axis=0, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
            for i_sca in range(n_sca):
                interp_funs['scatterer'].append(lambda pos, medium_vals=medium['scatterer'][mask,i_sca] : interp1d(R, medium_vals, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
                xsec_sca = interpolate_wavelengths(medium['scattering_cross_section'][mask,:,i_sca], wavelengths_in, wavelengths_out)
                interp_funs['scattering_cross_section'].append(lambda pos, medium_vals=xsec_sca : interp1d(R, medium_vals, axis=0, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
            for i_emi in range(n_emi):
                interp_funs['emitter'].append(lambda pos, medium_vals=medium['emitter'][mask,i_emi] : interp1d(R, medium_vals, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
                emiss = interpolate_wavelengths(medium['medium_emissivity'][mask,:,:,i_emi], wavelengths_in, wavelengths_out)
                interp_funs['medium_emissivity'].append(lambda pos, medium_vals=emiss : interp1d(R, medium_vals, axis=0, kind=interpolator,assume_sorted=True)(base_interp_fun(pos)))
            if optimized_mode:
                R_poss = medium['position'][mask]
                tau_abs = np.zeros((R_poss.shape[0],wavelengths_in.size))
                tau_sca = np.zeros((R_poss.shape[0],wavelengths_in.size))

                for R_idx in range(R_poss.shape[0]):
                    R_pos = R_poss[R_idx]
                    for idx_abs in range(-n_abs,0):
                        tau_abs[R_idx,:] = tau_abs[R_idx,:] + (
                          interp_funs['absorber'][idx_abs](R_pos)
                        * interp_funs['absorbing_cross_section'][idx_abs](R_pos) )
                    for idx_sca in range(-n_sca,0):
                        tau_sca[R_idx,:] = tau_sca[R_idx,:] + (
                          interp_funs['scatterer'][idx_sca](R_pos)
                        * interp_funs['scattering_cross_section'][idx_sca](R_pos) )
                tau_ext = cm_in_km * (tau_abs + tau_sca)
                interp_funs['extinction'].append(lambda pos, medium_vals=tau_ext : interp1d(R, medium_vals, axis=0, kind=interpolator,copy=False,assume_sorted=True)(base_interp_fun(pos)))
        elif current_fun_selection == 2:
            #2: Interpolation is radially from the medium position.
            #    There's no dependence between different medium positions. The basis
            #    function is a 3-dimensional Gaussian hump. The interpolation parameter
            #    is the standard deviation of the Gaussian.

            # Inversiossa sitten pitääpi tehä vaan siroamattomalla
            medium_std = medium['interpolation_parameter'][mask]
            medium_pos = medium['position'][mask]

            for i_abs in range(n_abs):
                interp_funs['absorber'].append(lambda pos, medium_vals=medium['absorber'][mask,i_abs] : gauss_fun(pos,medium_pos,medium_std,medium_vals))
                xsec_abs = interpolate_wavelengths(medium['absorbing_cross_section'][mask,:,i_abs], wavelengths_in, wavelengths_out)
                interp_funs['absorbing_cross_section'].append(lambda pos, medium_vals=xsec_abs : gauss_fun(pos,medium_pos,medium_std,medium_vals))
            for i_sca in range(n_sca):
                interp_funs['scatterer'].append(lambda pos, medium_vals=medium['scatterer'][mask,i_sca] : gauss_fun(pos,medium_pos,medium_std,medium_vals))
                xsec_sca = interpolate_wavelengths(medium['scattering_cross_section'][mask,:,i_sca], wavelengths_in, wavelengths_out)
                interp_funs['scattering_cross_section'].append(lambda pos, medium_vals=xsec_sca : gauss_fun(pos,medium_pos,medium_std,medium_vals))
            for i_emi in range(n_emi):
                interp_funs['emitter'].append(lambda pos, medium_vals=medium['emitter'][mask,i_emi] : gauss_fun(pos,medium_pos,medium_std,medium_vals))
                emiss = interpolate_wavelengths(medium['medium_emissivity'][mask,:,:,i_emi], wavelengths_in, wavelengths_out)
                interp_funs['medium_emissivity'].append(lambda pos, medium_vals=emiss : gauss_fun(pos,medium_pos,medium_std,medium_vals))
            if optimized_mode:
                R_poss = medium['position'][mask]
                tau_abs = np.zeros((R_poss.shape[0],wavelengths_in.size))
                tau_sca = np.zeros((R_poss.shape[0],wavelengths_in.size))

                for R_idx in range(R_poss.shape[0]):
                    R_pos = R_poss[R_idx]
                    for idx_abs in range(-n_abs,0):
                        tau_abs[R_idx,:] = tau_abs[R_idx,:] + (
                          interp_funs['absorber'][idx_abs](R_pos)
                        * interp_funs['absorbing_cross_section'][idx_abs](R_pos) )
                    for idx_sca in range(-n_sca,0):
                        tau_sca[R_idx,:] = tau_sca[R_idx,:] + (
                          interp_funs['scatterer'][idx_sca](R_pos)
                        * interp_funs['scattering_cross_section'][idx_sca](R_pos) )
                tau_ext = cm_in_km * (tau_abs + tau_sca)
                interp_funs['extinction'].append(lambda pos, medium_vals=tau_ext : gauss_fun(pos,medium_pos,medium_std,medium_vals))

        elif current_fun_selection == 3:
            pass
            #nearest neighbour blocks
            #absorber and scatterer vary in 3D
            #absorbing and scattering cross-sections vary in 1D
        else:
            raise ValueError("No interpolation function '%d' defined in raysca!" % current_fun_selection)
    return interp_funs

def muller_rayleigh(theta, depolarization):
    """
    Rayleigh scattering with depolarization factor.
    Chandrasekhar p. 49
    Originally in Siro by Liisa Oikarinen on 8.12.1998
    """
    if type(depolarization) == type(np.array([1])):
        # NOTE: This is a safeguard due to some curious circumstances and
        # it shouldn't slow the code too much down.
        depolarization = depolarization[0]
    gamma = depolarization / (2.0 - depolarization)
    coeff = 3.0 / (8.0 * np.pi) * (1.0 / (1.0 + 2.0 * depolarization))
    M11 = (np.cos(theta) ** 2 * (1.0 - gamma) + gamma)
    M33 = coeff * np.cos(theta) * (1.0 - gamma)
    M44 = coeff * np.cos(theta) * (1.0 - 3.0 * gamma)
    S = np.array([[coeff * M11  , gamma * coeff, 0.0        , 0.0        ],
                  [gamma * coeff, coeff        , 0.0        , 0.0        ],
                  [0.0          , 0.0          , coeff * M33, 0.0        ],
                  [0.0          , 0.0          , 0.0        , coeff * M44]])
    return S

def phase_function_rayleigh(theta):
    return 0.75 / (4.0 * np.pi) * (1 + np.cos(theta) ** 2)

def phase_rayleigh_wrapper(dir_in,dir_out,wl,sca_param):
    #sca_param here doesn't do anything yet
    theta = np.arccos(np.dot(dir_in,dir_out))
    ph = phase_function_rayleigh(theta)
    return ph

def muller_rayleigh_wrapper(dir_in,dir_out,wl,sca_param):
    theta = np.arccos(np.dot(dir_in,dir_out))
    S = muller_rayleigh(theta, sca_param)
    return (S, dir_out)

def muller_mie(theta, muller_data):
    costheta = np.cos(theta)
    costheta_idx = np.argmax(muller_data[:,0] >= costheta)
    Mu = muller_data[costheta_idx,:]
    # TODO: interpolation here, if needed
    S = np.array([[Mu[1],  Mu[2],  Mu[3],  Mu[4]],
                  [Mu[5],  Mu[6],  Mu[7],  Mu[8]],
                  [Mu[9],  Mu[10], Mu[11], Mu[12]],
                  [Mu[13], Mu[14], Mu[15], Mu[16]]])
    return 1 / (2 * np.pi) * S

def phase_function_mie(theta, muller_data):
    costheta = np.cos(theta)
    costheta_idx = np.argmax(muller_data[:,0] >= costheta)
    # The muller data ought to be in IIUV basis, so the phase function is
    # given by the sum of M11 and M22.
    phase = 1 / (4 * np.pi) * muller_data(costheta_idx,1) + muller_data(costheta_idx,6)
    return phase

#NOTE: I don't recall what was the design reason for the wrappers for these
# functions. Refactor 'em out from 2020.2.1 onwards, if this comment hadn't
# changed.

def phase_mie_wrapper(dir_in,dir_out,wl,muller_data):
    theta = np.arccos(np.dot(dir_in,dir_out))
    ph = phase_function_mie(theta, muller_data)
    return ph

def muller_mie_wrapper(dir_in,dir_out,wl,muller_data):
    theta = np.arccos(np.dot(dir_in,dir_out))
    S = muller_mie(theta, muller_data)
    return (S, dir_out)
#TODO: actually the phase function ought to return the dir_out

def rotation_matrix(cosphi,sinphi):
    cosphi2 = cosphi ** 2
    sinphi2 = sinphi ** 2
    cos2phi = 2.0 * cosphi2 - 1.0
    sin2phi = 2.0 * sinphi * cosphi

    R_mat = np.array(
            [[ cosphi2,sinphi2, 0.5 * sin2phi,0.0],
             [ sinphi2,cosphi2,-0.5 * sin2phi,0.0],
             [-sin2phi,sin2phi,       cos2phi,0.0],
             [     0.0,    0.0,           0.0,1.0]])

    return R_mat

def forward_rotation(dir_in,dir_out):
    # This rotates the polarization vectors so that they are properly aligned

    #TODO: Define eps and numacc in the configuration file?
    # Should the be the same?
    numacc = 5e-7

    # dir_in is parallel to the z-axis
    cond_parallel_z = np.abs(np.abs(dir_in[2]) - 1.0) < numacc

    # dir_in and dir_out are parallel
    cond_parallel = np.dot(dir_in - dir_out, dir_in - dir_out) < numacc
    cond_antiparallel = np.dot(dir_in + dir_out, dir_in + dir_out) < numacc

    if cond_parallel_z or cond_parallel or cond_antiparallel:
        return np.eye(4)
    else:
        # TODO: Figure out what is going on in here
        dotprod = np.dot(dir_in, dir_out)
        proj1 = 1.0 / np.sqrt(1.0 - dotprod ** 2)
        proj2 = 1.0 / np.sqrt(dir_in[0] ** 2 + dir_in[1] ** 2)
        norma = proj1 * proj2
        cosphi = norma * (dir_out[0] * dir_in[1] - dir_out[1] * dir_in[0])
        sinphi = norma * (-dir_out[2] + dir_in[2] * dotprod)
        return rotation_matrix(cosphi,sinphi)

def rotate_polarization(dir_in,dir_out,S,phasef,theta,cosalpha,sinalpha):
    # Currently this isn't used anywhere, it is here just for future multiple
    # scattering implementations.
    # This is the function polaris from polarisation.f90
    eps = 1e-12
    costheta = np.cos(theta)
    cond_forward_sca = np.dot(dir_in - dir_out, dir_in - dir_out) < eps
    cond_backward_sca = np.dot(dir_in + dir_out, dir_in + dir_out) < eps

    if cond_forward_sca or cond_backward_sca:
        rotated_S = S
    else:
        if np.abs(np.abs(dir_in[2]) - 1.0) < eps:
            # The incoming direction is parallel to the z-axis
            cosbeta = 0.0
            sinbeta = 1.0
        elif np.abs(dir_out[2] - 1.0) < eps:
            # The outgoing radiation is toward the positive z-axis
            norm = np.sqrt(1.0 - dir_in[2] ** 2)
            cosbeta = -dir_in[0] / norm
            sinbeta =  dir_in[1] / norm
        elif np.abs(dir_out[2] + 1.0) < eps:
            # The outgoing radiation is toward the negative z-axis
            norm = np.sqrt(1.0 - dir_in[2] ** 2)
            cosbeta = -dir_in[0] / norm
            sinbeta = -dir_in[1] / norm
        else:
            determinant = dir_in[0] * dir_out[1] - dir_in[1] * dir_out[0]
            if np.abs(determinant) < eps:
                cosbeta = 0.0
                sinbeta = 1.0 # TODO: In the original Siro implementation
                # this was commented with "should this be -1?". Research why
                # it should be so.
            else:
                cosbeta = -cosalpha * np.sqrt(dir_in[0] ** 2 + dir_in[1] ** 2) / np.sqrt(dir_out[0] ** 2 + dir_out[1] ** 2)
                sinbeta = cosbeta * (dir_in[2] - dir_out[2] * costheta) / determinant

        R_alpha = rotation_matrix(cosalpha,sinalpha)
        R_beta = rotation_matrix(cosbeta,sinbeta)
        rotated_S = R_alpha @ S @ R_beta / phasef(theta)
    return rotated_S

def set_up_scattering(medium):
    """
    Sets up scattering functions for each of the scatterers. This could include
    setting up look-up-tables from files for example.

    Returns a list of all scattering functions used in the computation.
    """

    sca_fun_idxs = medium['scatterer_kernel'][:]
    sca_fun_params = medium['scatterer_kernel_parameter'][:,:]
    sca_funs = []

    for (fun_idx,current_fun_selection) in enumerate(sca_fun_idxs):
        if current_fun_selection == 0:
            # Rayleigh scattering
            rayleigh_muller = lambda dir_in, dir_out, wl, sca_param=sca_fun_params[:,fun_idx] : muller_rayleigh_wrapper(dir_in,dir_out,wl,sca_param)
            rayleigh_phase = lambda dir_in, dir_out, wl, sca_param=sca_fun_params[:,fun_idx] : phase_rayleigh_wrapper(dir_in,dir_out,wl,sca_param)
            sca_funs.append([rayleigh_muller, rayleigh_phase])
        elif current_fun_selection == 1:
            # Phase matrix of the aerosol
            muller_data_in = np.genfromtxt(phase_folder + aerosol_file)
            mie_muller = lambda dir_in, dir_out, wl, muller_data=muller_data_in : muller_mie_wrapper(dir_in,dir_out,wl,muller_data)
            mie_phase = lambda dir_in, dir_out, wl, muller_data=muller_data_in : phase_mie_wrapper(dir_in,dir_out,wl,muller_data)
            sca_funs.append([mie_muller, mie_phase])
        else:
            raise ValueError("No scattering function '%d' defined in raysca!" % current_fun_selection)

    return sca_funs

# =====================
# Monte Carlo functions
# =====================

def sample_scattering_direction():
    #this will select the scatterer and sample a new scattering direction
    pass

def get_scattering_mc_table():
    pass

def sample_reflection_direction():
    pass

def get_reflection_mc_table():
    pass

def sample_scattering_position(path_dict):
    if path_dict['cumulate']:
        pass
    else:
        pass


def get_path(geom, interp_funs, refl_funs, photon_position, photon_direction, photon_in_domain, fast_trace):
    #Generates a path array which stores the photon positions, directions,
    #transmissivities (weights) and scattering probabilities
    path_dict = {'position' : [],
                 'direction' : [],
                 'transmissivty' : [],
                 'cumulate' : not fast_trace, # if this is the main beam, then we'll cumulate the transmissivty beforehand
                 'scattering_efficiency' : [],
                 'cumulative_scattering_probability' : [],
                 'reflect_at_end' : False,
                 'boundary_at_end' : 0}
    # if reflect_at_end is False, then the final element of the lists is
    # an intramedium scattering event (presumably right before a beam exits
    # the medium). boundary_at_end stores the boundary index
    transmissivity = np.ones_like(geom['incident_stokes'])
    idx_wl = geom['wavelength_index']
    tau_size = idx_wl.size
    n_sca = len(interp_funs['scatterer'])

    while True:
        (step_transmissivity, next_position, photon_direction) = trace_photon(
                    geom, interp_funs, refl_funs,
                    photon_position, photon_direction,
                    photon_in_domain, single_step=True,
                    fast_trace=fast_trace)
        if not photon_in_domain(photon_position):
            b_idx = find_crossed_boundary(geom['boundary'], photon_position, next_position)
            path_dict['boundary_at_end'] = b_idx
            if refl_funs[b_idx].refl_type != 0:
                path_dict['reflect_at_end'] = True
            break
        else:
            photon_position = next_position
        path_dict['position'].append(photon_position)
        # for example refraction may change the photon direction in
        # trace_photon
        path_dict['direction'].append(photon_direction)
        transmissivity = transmissivity * step_transmissivity
        path_dict['transmissivity'].append(transmissivity)
        tau_sca = np.zeros((n_sca,tau_size))
        for idx_sca in range(n_sca):
            tau_sca[idx_sca,:] = (interp_funs['scatterer'][idx_sca](photon_position)
            * interp_funs['scattering_cross_section'][idx_sca](photon_position)[idx_wl])
        path_dict['scattering_efficiency'].append(tau_sca)
        probs = []
        for idx_sca in range(n_sca):
            # NOTE: The mean of the scattering is considered here
            prob = np.mean(interp_funs['scatterer'][idx_sca](photon_position)
            * interp_funs['scattering_cross_section'][idx_sca](photon_position)[idx_wl])
            probs.append(prob)
        path_dict['cumulative_scattering_probability'].append(probs)
    cum_sca_prob = np.array(path_dict['cumulative_scattering_probability'])
    cum_sca_prob_sum = np.sum(cum_sca_prob,axis=1)

    if path_dict['reflect_at_end']:
        step_len = scattering_step_length if fast_trace else main_beam_step_length
        steps = step_len * np.ones((len(path_dict['position']),))
        medium_scatter_probability = 1 - np.exp(-np.sum(cum_sca_prob_sum * steps))

        # Append the final step for surface reflectance
        path_dict['cumulative_scattering_probability'].append(medium_scatter_probability)
        (step_transmissivity, last_position, last_direction) = trace_photon(
                    geom, interp_funs, refl_funs,
                    photon_position, photon_direction,
                    photon_in_domain, single_step=True,
                    fast_trace=fast_trace)
        path_dict['transmissivity'].append(step_transmissivity)
        path_dict['position'].append(last_position)
        path_dict['direction'].append(last_direction)
    else:
        # The boundary at the end doesn't reflect the photon, so we'll force a
        # scattering event in the atmosphere.
        medium_scatter_probability = 1

    for k in ['position','direction','transmissivity','scattering_efficiency','cumulative_scattering_probability']:
        path_dict[k] = np.array(path_dict[k])

    path_dict['cumulative_scattering_probability'] = path_dict['cumulative_scattering_probability'] / medium_scatter_probability

    if path_dict['cumulate']:
        #TODO: cumsum probability and then scale!
        #TODO: cumprod transmissivity!
        #tsekkaa että mitkä dimensiot on kyseessä
        pass
    return path_dict

def set_up_geometry(source_idx,los_idx,wl_idx,instrument,source,boundary):
    """
    Sets up all the needed geometry of the ray-tracer. The instrument locations
    and viewing beams, source directions and locations and the computation domain
    boundaries.
    """

    geom = {}

    shapes = boundary['shape'][:]
    shape_params = boundary['parameter'][:]

    boundary_funs = []
    boundary_normal_funs = []

    # NOTE: This is done now in every step, but in actuality, this should be done
    # just once. Boundaries are the same for each of the tasks.
    for shape_idx, shape in enumerate(shapes):
        if shape == 1:
            # spherical surface
            # The function will return a number larger than zero if the point is
            # on the outside and a negative number when we're inside the sphere.
            # The surface normal should point always INSIDE the medium.
            shape_param = shape_params[shape_idx]
            boundary_funs.append(lambda r, shape_param=shape_param : np.dot(r,r) - shape_param ** 2)

            if shape_param == np.max(shape_params):
                # This is the case for the outermost spherical boundary
                boundary_normal_funs.append(lambda r : -r / np.linalg.norm(r))
            else:
                boundary_normal_funs.append(lambda r : r / np.linalg.norm(r))

        else:
            raise ValueError("No shape function '%d' defined in raysca!" % shape)

    geom['boundary'] = boundary_funs
    geom['boundary_normal'] = boundary_normal_funs

    source_type = source['type'][source_idx]
    if source_type == 0:
        #far-field source
        #the source is in the same direction everywhere
        source_dir_fun = lambda r, src_idx=source_idx : -source['incident_direction'][src_idx]
        geom['source_dir'] = source_dir_fun
    else:
        raise ValueError("No source function '%d' defined in raysca!" % source_type)

    # The idea of this function is to determine if the ray hits directly to the
    # source when tracing it from point r into the direction dir
    radius_s = source['source_angular_radius'][source_idx] # this is in degrees, so it gets converted to radians here --v
    source_direct_illumination = lambda pho_dir, r, radius_s=radius_s : np.arccos(np.dot(pho_dir,source_dir_fun(r))) < (radius_s * 2 * np.pi / 360)
    geom['source_direct_illumination'] = source_direct_illumination

    if single_trace_for_all_wl:
        geom['incident_stokes'] = source['incident_stokes'][source_idx,:,:]
        geom['wavelength'] = source['output_wavelength'][:]
        geom['wavelength_index'] = np.arange(source['output_wavelength'][:].size)
    else:
        geom['incident_stokes'] = source['incident_stokes'][source_idx,wl_idx,:]
        geom['wavelength'] = source['output_wavelength'][wl_idx]
        geom['wavelength_index'] = wl_idx

    geom['incident_stokes'] = geom['incident_stokes'] @ IQUV2IIUV

    geom['view_vector'] = instrument['view_vector'][los_idx,:]
    geom['instrument_position'] = instrument['position'][los_idx,:]

    return geom

def find_crossed_boundary(boundary_funs, current_position, next_position):
    for bfun_idx, bfun in enumerate(boundary_funs):
        if bfun(current_position) * bfun(next_position) <= 0:
            return bfun_idx
    else:
        # This means that the we are not in the domain anymore, but we didn't
        # find the boundary we crossed. Let's raise some errors.
        # Disable the error for a while t. Antti 7.11.2021
        raise ValueError("Couldn't find the crossed boundary! Ray is lost!")
        return -1 # This is to assure that something gets returned if no boundary
                  # gets crossed.

def move_photon_to_boundary(b_fun, bn_fun, current_position, next_position):
    dirvec = next_position - current_position
    seg_fun = lambda t: b_fun(t * dirvec + current_position) ** 2 # the second power, because it is minimization
    t0 = minimize_scalar(seg_fun,bounds=[0,1])
    current_position = t0.x * dirvec + current_position
    # current_position is now exactly at the boundary
    # Let's move it inside just a little bit
    bound_tol = 1e-6
    return bound_tol * bn_fun(current_position) + current_position


def cartesian_product(*arrays):
    #NOTE: This is identical definition to ARSCA's utilities module
    #The rightmost array is the "innermost", so that cartprod[i] and
    #cartprod[i+1] stand for [x_j,y_k,z_l] and [x_j,y_k,z_l+1] respectively

    #This is from https://stackoverflow.com/questions/11144513/cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def stokes_times_muller(stokes,M):
    # stokes.shape = (n_wl,n_stokes)
    # M.shape = (n_wl,n_stokes,n_stokes)
    # NOTE: This is an ugly hack, and a desperate fix of some problems in other
    # places.
    if len(stokes.shape) < 2:
        stokes = stokes.reshape((1,stokes.size))
    out_stokes = np.nan * np.ones_like(stokes)
    for i in range(stokes.shape[0]):
        out_stokes[i,:] = stokes[i,:] @ M[i,:,:].T
    return out_stokes


def set_up_RT_tasks(n_source, n_los, n_wavelength, n_compute_cores):
    """
    Divides the computation into independent tasks. These tasks may then
    be divided onto different threads. Solar geometries, lines of sight and
    wavelengths are independent

    Returns a list of lists of tuples of (solar_geom_idx, los_idx, wl_idx).
    The outermost list index stands for the index of the computation core.
    The second list index stands for the task number of the current core.
    The tuple contains the geometry and wavelength indices for the current
    task.

    If the setting single_trace_for_all_wl is enabled, wavelengths aren't used
    in the task generation.

    If the setting source_los_cartesian_product is enabled, then the each of
    different source and line of sight pairs are computed separately.
    """

    if single_trace_for_all_wl:
        wl_range = np.array([-1])
    else:
        wl_range = np.arange(n_wavelength)

    all_tasks = cartesian_product(np.arange(n_source),
                                  np.arange(n_los),
                                  wl_range)

    if not source_los_cartesian_product:
        # we'll slice the non-necessary elements out
        all_tasks = all_tasks[np.where(all_tasks[:,0] == all_tasks[:,1])]
    n_tasks = len(all_tasks)
    compute_core_ids = np.array(range(n_tasks)).reshape((n_tasks,1)) % n_compute_cores
    #The permutation one will distribute the tasks randomly onto different
    #threads
    #work_equalizer_idx = np.random.permutation(n_tasks)
    work_equalizer_idx = np.arange(n_tasks)
    all_tasks_mixed = all_tasks[work_equalizer_idx,:]
    all_tasks_array = np.hstack((np.array(all_tasks_mixed),compute_core_ids))

    return all_tasks_array

def trace_photon(geom, interp_funs, refl_funs,
                 photon_position, photon_direction,
                 photon_in_domain, single_step=False,
                 fast_trace=False):

    if not fast_trace:
        step_length = main_beam_step_length
    else:
        step_length = scattering_step_length

    idx_wl = geom['wavelength_index']

    transmissivity = np.ones_like(geom['incident_stokes'])
    while photon_in_domain(photon_position):
        current_position = photon_position
        next_position = photon_position + step_length * photon_direction
        if not photon_in_domain(next_position):
            # We need to check if the boundary blocks radiation and to do that,
            # we have to figure out which boundary was crossed.
            b_idx = find_crossed_boundary(geom['boundary'], current_position, next_position)
            # NOTE: boundary_normal on the next line should probably be called with an
            # exact point of crossing, which could be computed. If this turns out
            # to be a problem, let's fix it.

            photon_exiting_position = move_photon_to_boundary(geom['boundary'][b_idx], geom['boundary_normal'][b_idx], photon_position, next_position)

            wl = geom['wavelength'] if geom['wavelength'].shape == () else geom['wavelength'][idx_wl]
            _, T, _ = refl_funs[b_idx].reflection(photon_direction, photon_direction, geom['boundary_normal'][b_idx](photon_exiting_position), wl)
            try:
                _ = T[0,0,0]
            except IndexError:
                print("bloody hell")
            transmissivity = update_transmissivity(transmissivity,interp_funs,geom,current_position,photon_exiting_position)
            transmissivity = stokes_times_muller(transmissivity, T)
            photon_position = next_position # this is to ensure that the compute_stokes knows to exit too
            break

        transmissivity = update_transmissivity(transmissivity,interp_funs,geom,current_position,next_position)
        photon_position = next_position
        if single_step:
            break
    return (transmissivity, photon_position, photon_direction)

def update_transmissivity(transmissivity,interp_funs,geom,current_position,next_position):
    step_length = np.linalg.norm(next_position - current_position)

    n_abs = len(interp_funs['absorber'])
    n_sca = len(interp_funs['scatterer'])
    n_emi = len(interp_funs['emitter'])

    idx_wl = geom['wavelength_index']
    if single_trace_for_all_wl:
        tau_size = idx_wl.size
        tau_idx = idx_wl
    else:
        tau_size = 1
        tau_idx = 0

    # The extinction
    # The absorption and scattering are averaged over the distance end-points.
    # This is equivalent of doing a linear interpolation over the distance
    # and integrating over it.
    tau_abs = np.zeros((1,tau_size))
    for idx_abs in range(n_abs):
        tau_abs[0,tau_idx] = tau_abs[0,tau_idx] + 0.5 * (
          interp_funs['absorber'][idx_abs](current_position)
        * interp_funs['absorbing_cross_section'][idx_abs](current_position)[idx_wl]
        + interp_funs['absorber'][idx_abs](next_position)
        * interp_funs['absorbing_cross_section'][idx_abs](next_position)[idx_wl] )
    tau_sca = np.zeros((n_sca,tau_size))
    for idx_sca in range(n_sca):
        tau_sca[idx_sca,:] = 0.5 * (
          interp_funs['scatterer'][idx_sca](current_position)
        * interp_funs['scattering_cross_section'][idx_sca](current_position)[idx_wl]
        + interp_funs['scatterer'][idx_sca](next_position)
        * interp_funs['scattering_cross_section'][idx_sca](next_position)[idx_wl] )
    tau_abs = cm_in_km * tau_abs
    tau_sca_tot = cm_in_km * np.sum(tau_sca,axis=0)
    tau_ext = tau_abs + tau_sca_tot

    transmissivity = np.exp(-tau_ext * step_length).repeat(4,axis=0).T * transmissivity
    return transmissivity

def compute_stokes(geom, interp_funs, sca_funs, refl_funs):
    """
    This is it. The main man. Head honcho. This guy'll handle all your problems,
    but only if they're RT problems.

    Computes single-los, single-source, single-wavelength RT problem in the
    geometry and medium defined in geom and interp_funs. Whatever changes the
    direction of the radiation is defined in sca_funs and refl_funs.
    """
    n_sca = len(interp_funs['scatterer'])
    idx_wl = geom['wavelength_index']

    wl = geom['wavelength'] if geom['wavelength'].shape == () else geom['wavelength'][idx_wl]

    stokes_in = geom['incident_stokes']
    stokes_out = np.zeros_like(stokes_in)

    photon_position = geom['instrument_position']
    photon_direction = geom['view_vector']

    # For the sake of simplicity, let's assume that all the medium is inside the
    # defined domain, which is bounded by the boundaries. That is, it doesn't get
    # scattered or attenuated anywhere else.

    # NOTE: This only works in cases with two boundary functions. That might cause
    # problems later on.
    photon_in_domain = lambda r : geom['boundary'][0](r) * geom['boundary'][1](r) < 0

    distance_travelled = 0.0

    step_length = main_beam_step_length

    next_position = np.nan * np.ones((3,))

    while not photon_in_domain(photon_position):
        # TODO: In this part, we should check also the exact crossing point
        # AND any effect which might come from crossing some certain boundary.
        next_position = photon_position + step_length * photon_direction
        distance_travelled = distance_travelled + step_length
        if distance_travelled > max_distance:
            #The photon missed the domain completely. Who aimed this ray?!
            if geom['source_direct_illumination'](photon_direction, photon_position):
                # if the photon hits the source without passing through the medium
                return stokes_in
            else:
                return stokes_out

        else:
            prev_position = photon_position
            photon_position = next_position

    b_idx = find_crossed_boundary(geom['boundary'], prev_position, photon_position)
    photon_position = move_photon_to_boundary(geom['boundary'][b_idx], geom['boundary_normal'][b_idx], prev_position, photon_position)

    if monte_carlo_mode:
        print("Create the initial Monte Carlo table here")
        raise Exception

    # Now the photon should be in the domain. Let's continue propagating it
    transmissivity_main_beam = np.ones_like(stokes_in)

    photon_exiting_domain = False
    while np.linalg.norm(transmissivity_main_beam) > extinction_threshold and not photon_exiting_domain:
        (transmissivity_step, next_position,_ ) = trace_photon(
                                               geom, interp_funs, refl_funs,
                                               photon_position, photon_direction,
                                               photon_in_domain, single_step=True,
                                               fast_trace=False)
        if not photon_in_domain(next_position):
            # If the photon escapes the domain, we'll stop tracing it
            # The photon has hit a boundary and next_position is outside the domain, while
            # photon_position is still on the inside.

            b_idx = find_crossed_boundary(geom['boundary'], photon_position, next_position)
            next_position = move_photon_to_boundary(geom['boundary'][b_idx], geom['boundary_normal'][b_idx], photon_position, next_position)
            if single_scattering:
                scattering_step_length = np.linalg.norm(next_position - photon_position)
            transmissivity_main_beam = update_transmissivity(transmissivity_main_beam,interp_funs,geom,photon_position,next_position)
            photon_position = next_position
            photon_exiting_domain = True

        else:
            if single_scattering:
                scattering_step_length = step_length
            transmissivity_main_beam = transmissivity_main_beam * transmissivity_step
            photon_position = next_position

        if single_scattering:
            # scattering peel-off
            source_dir = geom['source_dir'](photon_position)
            (transmissivity_peel_off,_ ,_ ) = trace_photon(
                                                   geom, interp_funs, refl_funs,
                                                   photon_position, source_dir,
                                                   photon_in_domain, single_step=False,
                                                   fast_trace=True)
            for idx_sca in range(n_sca):
                # Propagation direction might get changed in the scattering.
                # Now it is just discarded, since we do the peel-off
                (S, _) = sca_funs[idx_sca][0](photon_direction, source_dir, wl)

                scattering_efficiency = (interp_funs['scattering_cross_section'][idx_sca](photon_position)[idx_wl]
                                       * interp_funs['scatterer'][idx_sca](photon_position)
                                       * cm_in_km
                                       * scattering_step_length)

                #This is to enable the multiplication below when handling several wavelengths at once.
                scattering_efficiency = scattering_efficiency.reshape((scattering_efficiency.size,1)).repeat(4,axis=1)

                R_fwd = forward_rotation(photon_direction, source_dir)
                stokes_scattered = transmissivity_main_beam * transmissivity_peel_off * stokes_in @ S.T @ R_fwd.T
                stokes_out = stokes_out + scattering_efficiency * stokes_scattered


    else: # while condition not longer satisfied
        if not photon_exiting_domain:
            #in this case, the beam is fully extincted while in the medium
            return stokes_out


    # In single-scatterings, we'll check the contribution of surface reflection
    # to the transmission, or if the beam exits toward the source.
    source_dir = geom['source_dir'](photon_position)

    R, T, _ = refl_funs[b_idx].reflection(photon_direction, source_dir, geom['boundary_normal'][b_idx](photon_position), wl)

    R_fwd = forward_rotation(photon_direction, source_dir)

    if geom['source_direct_illumination'](photon_direction, photon_position):
        stokes_pass_through = transmissivity_main_beam * stokes_in @ T.T
        stokes_reflected = np.zeros_like(stokes_in)
    else:
        stokes_pass_through = np.zeros_like(stokes_in)
        (transmissivity_reflection,_ ,_ ) = trace_photon(
                                               geom, interp_funs, refl_funs,
                                               photon_position, source_dir,
                                               photon_in_domain, single_step=False,
                                               fast_trace=False)
        stokes_reflected = stokes_times_muller(transmissivity_main_beam * transmissivity_reflection * stokes_in, R) @ R_fwd.T

    stokes_out = stokes_out + stokes_reflected + stokes_pass_through
    return stokes_out

class Path:
    def __init__(self,segment_ids,scatter_event):
        #scatter_event:
        # dict containning scatterer/reflector choice, dir_in, dir_out, wl, surf_norm, param
        # wl dependent phase matrix should be doable, no?
        # sehän on sama juttu kuin BRDF:ssä!
        self.computed = False
        self.segment_ids = segment_ids # Only two segments in a RaySca path
        self.scatter_event = scatter_event # The scatter event joins the segments

    def compute_transmissivity_direct(self,segments,subsegments,geom,interp_funs,sca_funs,refl_funs):
        print("Implement this!")
        #lol, eipä tarvinnukkaa

    def compute_transmissivity(self,mb_segments,sca_segments,subsegments,geom,interp_funs,sca_funs,refl_funs):
        if len(self.segment_ids) != 2:
            # In this case, we should have only one segment, which means
            # we're doing an limb occlusion or direct sunlight measurement
            #self.compute_transmissivity_direct(segments,subsegments,geom,interp_funs,sca_funs,refl_funs)
            print("This shouldn't happen.")
            return
        seg0 = mb_segments[self.segment_ids[0]]
        seg0.compute_transmissivity(geom,subsegments,mb_segments,interp_funs)
        seg1 = sca_segments[self.segment_ids[1]]
        seg1.compute_transmissivity(geom,subsegments,sca_segments,interp_funs)

        stokes_in = geom['incident_stokes'][0,:]
        #NOTE: ASSUME ONLY UNIFORMLY UNPOLARIZED INCIDENT RADIATION!
        stokes_out = np.zeros_like(geom['incident_stokes'])
        photon_direction = self.scatter_event['dir_in']
        photon_position = self.scatter_event['photon_position']
        source_dir = self.scatter_event['dir_out']

        sca_wl = np.array([geom['wavelength'][0], geom['wavelength'][geom['wavelength'].size // 2], geom['wavelength'][-1]])
        #NOTE: ASSUME ONLY 3 DIFFERENT MULLER MATRICES IN THE BAND
        R_fwd = forward_rotation(photon_direction, source_dir)

        if self.scatter_event['type'] == 'reflection':
            b_idx = self.scatter_event['b_idx']
            R, _, _ = refl_funs[b_idx].reflection(photon_direction, source_dir, geom['boundary_normal'][b_idx](photon_position), sca_wl)

            stokes = R_fwd @ R @ stokes_in
            stokes_interp = interp1d(sca_wl,stokes,axis=0,kind='linear')
            #NOTE: WHEN DOING LINEAR INTERPOLATION, THE RESULTING STOKES MAY BE INTERPOLATED
            #IF MORE SMOOTH INTERPOLATION IS DESIRED, THEN THE MULLER NEEDS TO BE INTERPOLATED
            stokes_out[:] = stokes_interp(geom['wavelength'])

        elif self.scatter_event['type'] == 'scattering':
            idx_sca = self.scatter_event['idx_sca']
            scattering_length = self.scatter_event['scattering_length']
            (S, _) = sca_funs[idx_sca % len(sca_funs)][0](photon_direction, source_dir, sca_wl)

            scattering_efficiency = (interp_funs['scattering_cross_section'][idx_sca](photon_position)[:]
                                   * interp_funs['scatterer'][idx_sca](photon_position)
                                   * cm_in_km
                                   * scattering_length)
            stokes = R_fwd @ S @ stokes_in

            if len(stokes.shape) == 1:
                # In this case sca_funs doesn't return a wavelength-dependent
                # scattering mueller matrix
                stokes = stokes.reshape((1,stokes.size)).repeat(sca_wl.size,axis=0)

            stokes_interp = interp1d(sca_wl,stokes,axis=0,kind='linear')
            #NOTE: WHEN DOING LINEAR INTERPOLATION, THE RESULTING STOKES MAY BE INTERPOLATED
            #IF MORE SMOOTH INTERPOLATION IS DESIRED, THEN THE MULLER NEEDS TO BE INTERPOLATED
            stokes_out[:] = (scattering_efficiency * stokes_interp(geom['wavelength']).T).T
        elif self.scatter_event['type'] == 'pass-through':
            # direct and limb measurements
            if geom['source_direct_illumination'](photon_direction,photon_position):
                stokes_out = geom['incident_stokes'].copy()
            else:
                stokes_out = np.zeros_like(geom['incident_stokes'])

        stokes_out = stokes_out.T
        np.multiply(seg0.transmissivity,stokes_out,out=stokes_out)
        if not self.scatter_event['type'] == 'pass-through':
            np.multiply(seg1.transmissivity,stokes_out,out=stokes_out)
        stokes_out = stokes_out.T
        self.transmissivity = stokes_out

class Segment:
    def __init__(self,subseg_ids,segment_ids):
        self.subseg_ids = subseg_ids
        self.segment_ids = segment_ids
        self.computed = False
    def compute_transmissivity(self,geom,subsegments,segments,interp_funs):
        if self.computed:
            return
        n_wl = geom['wavelength'].size
        transmissivity = np.ones(n_wl)
        for idx in self.segment_ids:
            seg = segments[idx]
            if not seg.computed:
                seg.compute_transmissivity(geom,subsegments,segments,interp_funs)
            np.multiply(transmissivity,seg.transmissivity,out=transmissivity)
        for idx in self.subseg_ids:
            subseg = subsegments[idx]
            if not subseg.computed:
                subseg.compute(interp_funs)
            np.multiply(transmissivity,subseg.transmissivity,out=transmissivity)
        self.transmissivity = transmissivity
        self.computed = True
    def release_subsegments(self,subsegments):
        for ssid in self.subseg_ids:
            subsegments[ssid].release()

class Subsegment:
    def __init__(self,r_start,r_end):
        self.computed = False
        self.deleted = False
        self.r_start = r_start
        self.r_end = r_end

    #from numba import jit
    #@jit(nopython=True)
    def compute(self,interp_funs):
        step_length = np.linalg.norm(self.r_end - self.r_start)

        if not super_accurate:
            tau = None
            for idx_ext in range(len(interp_funs['extinction'])):
                if type(tau) == type(None):
                    tau = ( interp_funs['extinction'][idx_ext](self.r_start)
                          + interp_funs['extinction'][idx_ext](self.r_end) )
                else:
                    tau = ( tau
                          + interp_funs['extinction'][idx_ext](self.r_start)
                          + interp_funs['extinction'][idx_ext](self.r_end) )

            self.transmissivity = np.exp(-0.5 * step_length * tau)#.repeat(4,axis=0).T
            # NOTE: If polarization-dependent transmissivity is needed, then
            # we need to uncomment the repeat method call.
        else:
            points = 10
            small_step = step_length / points
            step_norm = (self.r_end - self.r_start) / step_length
            r_curr = self.r_start
            ext_tot = small_step * interp_funs['extinction'](self.r_start)
            for s in range(1,points):
                r_curr = s * small_step * step_norm + self.r_start
                ext_tot = ext_tot + small_step * interp_funs['extinction'](r_curr)
            self.transmissivity = np.exp(-ext_tot)
        self.computed = True
    def release(self):
        del(self.transmissivity)
        self.deleted = True

def trace_photon_optimized(geom,
                 photon_position, photon_direction,
                 photon_in_domain,
                 step_length):

    current_position = photon_position
    next_position = photon_position + step_length * photon_direction
    boundary_info = {}
    if not photon_in_domain(next_position):
        # We need to check if the boundary blocks radiation and to do that,
        # we have to figure out which boundary was crossed.
        b_idx = find_crossed_boundary(geom['boundary'], current_position, next_position)
        # NOTE: boundary_normal on the next line should probably be called with an
        # exact point of crossing, which could be computed. If this turns out
        # to be a problem, let's fix it.

        photon_exiting_position = move_photon_to_boundary(geom['boundary'][b_idx], geom['boundary_normal'][b_idx], photon_position, next_position)

        next_position = photon_exiting_position # this is to ensure that the compute_stokes knows to exit too
        boundary_info['b_idx'] = b_idx

    return (next_position, boundary_info)

def get_main_path_event(refl_funs,b_idx):
    """
    This is to determine which event the main path does. Case-dependent, it will
    be
    'reflection' : in nadir and glint observation modes
    'pass-through' : in limb scattering and direct modes, as well as in ground-
    based direct and open sky observation modes
    """
    if refl_funs[b_idx].refl_type == 0:
        # this is the pass-through boundary
        main_path_event = 'pass-through'
    else:
        main_path_event = 'reflection'
    return main_path_event


import copy
def compute_stokes_optimized(geom, interp_funs, sca_funs, refl_funs):
    """
    The idea here is to first trace all the paths through the atmosphere. This
    is easy in the single-scattering case. After all the paths are traced, then
    transmissivity is computed afterwards.
    """

    subsegments = []
    main_beam_segments = []
    scatter_segments = []
    paths = []

    n_sca = len(interp_funs['scatterer'])
    idx_wl = geom['wavelength_index']

    wl = geom['wavelength'] if geom['wavelength'].shape == () else geom['wavelength'][idx_wl]

    stokes_in = geom['incident_stokes']
    stokes_out = np.zeros_like(stokes_in)

    photon_position = geom['instrument_position']
    photon_direction = geom['view_vector']

    # For the sake of simplicity, let's assume that all the medium is inside the
    # defined domain, which is bounded by the boundaries. That is, it doesn't get
    # scattered or attenuated anywhere else.

    # NOTE: This only works in cases with two boundary functions. That might cause
    # problems later on.
    photon_in_domain = lambda r : geom['boundary'][0](r) * geom['boundary'][1](r) < 0

    distance_travelled = 0.0

    step_length = main_beam_step_length

    next_position = np.nan * np.ones((3,))

    while not photon_in_domain(photon_position):
        # TODO: In this part, we should check also the exact crossing point
        # AND any effect which might come from crossing some certain boundary.
        next_position = photon_position + step_length * photon_direction
        distance_travelled = distance_travelled + step_length
        if distance_travelled > max_distance:
            #The photon missed the domain completely. Who aimed this ray?!
            if geom['source_direct_illumination'](photon_direction, photon_position):
                # if the photon hits the source without passing through the medium
                return stokes_in
            else:
                return stokes_out

        else:
            prev_position = photon_position
            photon_position = next_position

    b_idx = find_crossed_boundary(geom['boundary'], prev_position, photon_position)
    photon_position = move_photon_to_boundary(geom['boundary'][b_idx], geom['boundary_normal'][b_idx], prev_position, photon_position)

    photon_exiting_domain = False

    while not photon_exiting_domain:

        (next_position, boundary_info) = trace_photon_optimized(geom,
                                                     photon_position, photon_direction,
                                                     photon_in_domain,
                                                     step_length)
        subsegments.append(Subsegment(photon_position,next_position))
        subsegid = len(subsegments) - 1
        mb_segid = len(main_beam_segments) - 1
        if len(main_beam_segments) == 0:
            main_beam_segments.append(Segment([subsegid],[]))
        else:
            main_beam_segments.append(Segment([subsegid],[mb_segid]))

        if 'b_idx' in boundary_info:
            # If the photon escapes the domain, we'll stop tracing it
            # The photon has hit a boundary and next_position is outside the domain, while
            # photon_position is still on the inside.
            b_idx = boundary_info['b_idx']
            if single_scattering:
                scattering_length = np.linalg.norm(next_position - photon_position)
            photon_exiting_domain = True
        else:
            if single_scattering:
                scattering_length = step_length
        photon_position = next_position

        if single_scattering:
            # scattering peel-off
            mb_segid = len(main_beam_segments) - 1
            #Note: mb_segid needs to be updated at this point in time due to
            #scattering Paths

            source_dir = geom['source_dir'](photon_position)
            scatter_event = {'type' : 'scattering',
                             'dir_in' : photon_direction,
                             'dir_out' : source_dir,
                             'scattering_length' : scattering_length,
                             'photon_position' : photon_position}
            boundary_info_scatter = {}
            current_position_sca = photon_position
            ss_subsegid_start = len(subsegments)
            while 'b_idx' not in boundary_info_scatter:
                (next_position_sca, boundary_info_scatter) = trace_photon_optimized(geom,
                                                         current_position_sca, source_dir,
                                                         photon_in_domain,
                                                         scattering_step_length)
                subsegments.append(Subsegment(current_position_sca,next_position_sca))
                current_position_sca = next_position_sca
                #NOTE: ASSUMED NO TERMINATOR EFFECT HERE
            ss_subsegid_end = len(subsegments)
            scatter_segments.append(Segment(range(ss_subsegid_start,ss_subsegid_end),[]))
            sca_segid = len(scatter_segments) - 1
            for idx_sca in range(n_sca):
                scatter_event['idx_sca'] = idx_sca
                paths.append(Path([mb_segid,sca_segid],copy.deepcopy(scatter_event)))

    #Create the final path
    source_dir = geom['source_dir'](photon_position)
    main_path_scatter_event = get_main_path_event(refl_funs,boundary_info['b_idx'])
    scatter_event = {'type' : main_path_scatter_event,
                              'dir_in' : photon_direction,
                              'dir_out' : source_dir,
                              'photon_position' : photon_position,
                              'b_idx' : boundary_info['b_idx']}
    boundary_info_ref = {}
    current_position_ref = photon_position
    ref_subsegid_start = len(subsegments)
    while 'b_idx' not in boundary_info_ref:
        #this moves the photon to the boundary
        (next_position_ref, boundary_info_ref) = trace_photon_optimized(geom,
                                                 current_position_ref, source_dir,
                                                 photon_in_domain,
                                                 main_beam_step_length)
        subsegments.append(Subsegment(current_position_ref,next_position_ref))
        current_position_ref = next_position_ref
        #NOTE: ASSUMED NO TERMINATOR EFFECT HERE
    ref_subsegid_end = len(subsegments)
    scatter_segments.append(Segment(range(ref_subsegid_start,ref_subsegid_end),[]))
    mb_segid = len(main_beam_segments) - 1
    sca_segid = len(scatter_segments) - 1
    paths.append(Path([mb_segid,sca_segid],scatter_event))

    #compute all path transmissivities here
    for seg in main_beam_segments:
        seg.compute_transmissivity(geom,subsegments,main_beam_segments,interp_funs)
        seg.release_subsegments(subsegments)
    if False:
        #This is somehow more slow than sequential????
        import pathos
        def process_wrapper(seg):
            seg.compute_transmissivity(geom,subsegments,main_beam_segments,interp_funs)
            seg.release_subsegments(subsegments)
        with pathos.multiprocessing.Pool(5) as p:
            p.map(process_wrapper,scatter_segments)
    else:
        for seg in scatter_segments:
            seg.compute_transmissivity(geom,subsegments,scatter_segments,interp_funs)
            seg.release_subsegments(subsegments)
    for seg in subsegments:
        if not seg.computed:
            print(seg)
    for p in paths:
        p.compute_transmissivity(main_beam_segments,
                    scatter_segments,subsegments,geom,interp_funs,
                    sca_funs,refl_funs)
        np.add(stokes_out,p.transmissivity,out=stokes_out)

    return stokes_out

def run_task_list(core_id, task_list, interp_funs, sca_funs, refl_funs, instrument, source, boundary,n_source, n_los, n_output_wavelength):
    """
    This is ran per single core.
    """

    n_stokes = 4
    output_stokes = np.zeros((n_source,
                              n_los,
                              n_output_wavelength,
                              n_stokes))
    tasks = task_list[task_list[:,3] == core_id]
    task_amt = len(tasks)

    if show_progress:
        pb = tqdm.tqdm(total=task_amt,desc='Thread %d' % core_id,position=core_id)
    for task in tasks:
        source_idx = task[0]
        los_idx = task[1]
        wl_idx = task[2]
        geom = set_up_geometry(source_idx,los_idx,wl_idx,instrument,source,boundary)

        if optimized_mode:
            stokes_f = compute_stokes_optimized
        else:
            stokes_f = compute_stokes
        stokes = stokes_f(geom,
                        interp_funs,
                        sca_funs,
                        refl_funs)

        if single_trace_for_all_wl:
            output_stokes[source_idx,los_idx,:,:] = stokes
        else:
            output_stokes[source_idx,los_idx,wl_idx,:] = stokes
        if show_progress:
            pb.update(1)
    if show_progress:
        pb.close()
    return output_stokes

def get_compute_cores():
    if parallelization:
        if process_count > 0:
            n_compute_cores = min([process_count,multiprocessing.cpu_count()])
        else:
            n_compute_cores = multiprocessing.cpu_count()
    else:
        n_compute_cores = 1
    return n_compute_cores

def task_wrapper(datatup):
    (core_id, filename) = datatup

    (medium, instrument, source, boundary) = read_input_nc(filename)
    n_compute_cores = get_compute_cores()
    n_source = source['type'][:].shape[0]
    n_los = instrument['view_vector'][:].shape[0]
    n_output_wavelength = source['output_wavelength'][:].shape[0]
    task_list = set_up_RT_tasks(n_source,n_los,n_output_wavelength,n_compute_cores)

    wavelengths_in = source['input_wavelength']
    wavelengths_out = source['output_wavelength']

    interp_funs = set_up_interpolation(medium, wavelengths_in, wavelengths_out)

    sca_funs = set_up_scattering(medium)

    refl_funs = set_up_reflection(boundary, phase_folder, brdf_file)

    return run_task_list(core_id, task_list, interp_funs, sca_funs, refl_funs, instrument, source, boundary,n_source, n_los, n_output_wavelength)

def run_simulation(filename):
    """
    This is the main function of raysca.
    """
    read_settings()

    (medium, instrument, source, boundary) = read_input_nc(filename)

    n_source = source['type'][:].shape[0]
    n_los = instrument['view_vector'][:].shape[0]
    n_output_wavelength = source['output_wavelength'][:].shape[0]

    n_compute_cores = get_compute_cores()

    n_stokes = 4
    output_stokes = np.zeros((n_source,
                              n_los,
                              n_output_wavelength,
                              n_stokes))

    if parallelization:
        core_ids = range(n_compute_cores)
        datatups = zip(core_ids,(filename,) * len(core_ids))

        with multiprocessing.Pool(n_compute_cores) as p:
            out_stokes_list = p.map(task_wrapper, datatups)
    else:
        out_stokes_list = [task_wrapper((0,filename))] # 0 is the core_id of the singular core

    for single_core_stokes in out_stokes_list:
        output_stokes = output_stokes + single_core_stokes

    return output_stokes @ IIUV2IQUV if stokes_output_format == 'IQUV' else output_stokes

#The two lines below are for running the debugger through Spyder
#out_stokes = run_simulation('../../solver_input/co2m_aerosol_test/single_spectrum.nc')
#out_stokes = run_simulation('ch4_rayleigh_spherical_low_res.nc')
#out_stokes = run_simulation('ch4_rayleigh_spherical-new.nc')
