#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:23:03 2021

@author: mikkonea
"""

import numpy as np

class Boundary:
    self.abs = []
    self.sca = []
    self.emi = []
    
    def __init__(self,absorbers,scatterers,emitters):
        self.abs = absorbers
        self.sca = scatterers
        self.emi = emitters
    
    def b(self):
        pass
    
    def bn(self):
        pass
        
def create_wl_stack(a,n):
    a_ = a[np.newaxis,:,:]
    return np.repeat(a_,n,axis=0)
    
def reflection_passthrough(dir_in,dir_out,surf_norm,wl,refl_param):
    if not np.all(np.isclose(dir_in,dir_out)):
        dir_out = dir_in
    # The pass-through surface doesn't reflect anything, but attenuates whatever goes
    # through it.
    T = (1.0 - refl_param) * np.eye(4)
    T = create_wl_stack(T,wl.size)
    R = np.zeros((wl.size,4,4))
    return (R, T, dir_out)

def reflection_semispecular(dir_in,dir_out,surf_norm,wl,refl_param):
    # This semispecular reflection is depolarizing. The reflection parameter
    # denotes the maximum radiation reflected by the surface (at the glint point).
    # Warning! This is just a quick hack of a reflection function! Energy may not
    # be conserved!
    # TODO: Check if energy is conserved
    diff_refl = 1.0 # The portion of diffusively reflected radiation
    spec_refl = 1.0 - diff_refl # The portion of specularly reflected radiation

    # The np.max is here just to desperately fix some cases where we would check
    # reflection through the boundary (in which case it should be transmission)
    lambertian_coeff = np.max([np.dot(surf_norm,dir_out) / np.pi,0])
    R_diff = np.zeros((4,4))
    R_diff[0,0] = 0.5
    R_diff[1,1] = 0.5
    R_diff[0,1] = 0.5
    R_diff[1,0] = 0.5
    R_diff = lambertian_coeff * R_diff

    R_spec = np.zeros((4,4))
    spec_dot_width_degrees = 2
    spec_dot_width = np.pi / 180 * spec_dot_width_degrees
    spec_refl_vec = dir_in - 2 * surf_norm * np.dot(surf_norm, dir_in)
    semispec_ref = np.exp(-(np.arccos(np.dot(spec_refl_vec,dir_out)) / spec_dot_width) ** 2)
    R_spec[0,0] = semispec_ref
    R_spec[1,1] = semispec_ref
    R_spec[0,1] = semispec_ref
    R_spec[1,0] = semispec_ref

    T = np.zeros((wl.size,4,4))
    R = refl_param * diff_refl * R_diff + spec_refl * R_spec
    R = create_wl_stack(R,wl.size)
    return (R, T, dir_out)

def reflection_lambertian(dir_in, dir_out, surf_norm,wl,refl_param): 
    lambertian_coeff = refl_param * np.max([np.dot(surf_norm,dir_out) / np.pi,0])
    R_diff = np.zeros((4,4))
    R_diff[0,0] = 0.5
    R_diff[1,1] = 0.5
    R_diff[0,1] = 0.5
    R_diff[1,0] = 0.5
    R = lambertian_coeff * R_diff
    T = np.zeros((wl.size,4,4))
    R = create_wl_stack(R,wl.size)
    return(R, T, dir_out)
    
def _get_projection_onto_plane(u,v):
    #yields some vector perpendicular to vector v if u is parallel to v
    #otherwise returns normalized u flattened to the plane defined by v 
    eps = 1e-12
    u_proj = u - np.dot(u,v) * v
    if np.linalg.norm(u_proj) < eps:
        c = np.cross(v,np.array([1.0,0.0,0.0]))
        if np.linalg.norm(c) < eps:
            return np.cross(v,np.array([0.0,1.0,0.0]))
        else:
            return c
    else:
        return u_proj / np.linalg.norm(u_proj)
    
    
def _get_azi_angle(dir_in, dir_out, surf_norm):
    #figuring out the surface tangential vector into which azimuth = 0
    #NOTE: The minus sign in the dir_in changes the viewing vector to 
    #correct base.
    dir_azi0 = _get_projection_onto_plane(-dir_in,surf_norm)
    dir_azi_out = _get_projection_onto_plane(dir_out,surf_norm)
    
    dir_perp_azi = np.cross(surf_norm,-dir_azi0)
    
    azi_ang = np.arccos(np.dot(dir_azi0,dir_azi_out))
    
    return azi_ang if np.dot(dir_perp_azi,dir_azi_out) >= 0 else 2 * np.pi - azi_ang 

def reflection_brdf(dir_in, dir_out, surf_norm, wl, refl_param):
    R = np.zeros((4,4))
    T = np.zeros((4,4))
    lambertian_coeff = np.max([np.dot(surf_norm,dir_out) / np.pi,0])
    
    zen_in = np.arccos(np.dot(surf_norm,-dir_in)) #NOTE: -dir_in points away from the boundary
    zen_out = np.arccos(np.dot(surf_norm,dir_out))
    if zen_out > np.pi / 2: # in this case the outgoing ray is transmitted, not reflected
        return(create_wl_stack(R,wl.size),create_wl_stack(T,wl.size),dir_out)
    azi_out = _get_azi_angle(dir_in, dir_out, surf_norm)
    
    intfun_zen_in = interp1d(brdf_muller['zen_inc'], brdf_muller['muller'], axis=0)
    temp_muller = intfun_zen_in(zen_in)
    
    intfun_wl = interp1d(brdf_muller['wavelength'],temp_muller,axis=0,fill_value='extrapolate')
    temp_muller = intfun_wl(wl)
    
    intfun_zen_out = interp1d(brdf_muller['zen_ref'],temp_muller,axis=1)
    temp_muller = intfun_zen_out(zen_out)
    
    intfun_azi_out = interp1d(brdf_muller['azi_ref'],temp_muller,axis=1)
    R = intfun_azi_out(azi_out) * lambertian_coeff
    
    T = np.zeros((wl.size,4,4))
    return(R,T,dir_out)

def set_up_reflection(boundary):
    """
    Sets up reflection functions for each of the boundary surfaces. This could
    include setting up look-up-tables from files for example.

    The reflection functions do the sampling (if applicable) of the reflection
    direction and the generation of reflection/transmission Muller matrices.

    Returns a list of all reflection functions used in the computation.
    """
    global brdf_muller
    refl_fun_idxs = boundary['reflection_kernel'][:]
    refl_fun_params = boundary['reflection_kernel_parameter'][:,:]
    refl_funs = []

    for (fun_idx,current_fun_selection) in enumerate(refl_fun_idxs):
        if current_fun_selection == 0:
            # Pass-through boundary
            passthrough_fun = lambda dir_in, dir_out, surf_norm, wl, refl_param=refl_fun_params[:,fun_idx] : reflection_passthrough(dir_in, dir_out, surf_norm, wl, refl_param)
            refl_funs.append(passthrough_fun)
        elif current_fun_selection == 1:
            # Lambertian reflection
            lambert_fun = lambda dir_in, dir_out, surf_norm, wl, refl_param=refl_fun_params[:,fun_idx] : reflection_lambertian(dir_in, dir_out, surf_norm, wl, refl_param)
            refl_funs.append(lambert_fun)
        elif current_fun_selection == 2:
            # Semispecular reflection
            semispec_fun = lambda dir_in, dir_out, surf_norm, wl, refl_param=refl_fun_params[:,fun_idx] : reflection_semispecular(dir_in, dir_out, surf_norm, wl, refl_param)
            refl_funs.append(semispec_fun)
        elif current_fun_selection == 3:
            # BRDF reflection
            import netCDF4 
            with netCDF4.Dataset(phase_folder + brdf_file) as ds:
                # Since the BRDF data may be big, let's store it into a global variable
                brdf_muller = {}
                for key in ds.variables.keys():
                    brdf_muller[key] = ds[key][:].data
            brdf_fun = lambda dir_in, dir_out, surf_norm, wl, refl_param=refl_fun_params[:,fun_idx] : reflection_brdf(dir_in, dir_out, surf_norm, wl, refl_param)
            refl_funs.append(brdf_fun)
        else:
            raise ValueError("No reflection function '%d' defined in raysca!" % current_fun_selection)

    return refl_funs


