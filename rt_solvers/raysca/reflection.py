#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:23:03 2021

@author: mikkonea
"""

import numpy as np
from scipy.interpolate import interp1d
import netCDF4

class BRDF:
    zen_inc = []
    wavelength = []
    zen_ref = []
    azi_ref = []
    muller = []
    
    def __init__(self,zen_inc,wavelength,zen_ref,azi_ref,muller):
        self.zen_inc = zen_inc
        self.wavelength = wavelength
        self.zen_ref = zen_ref
        self.azi_ref = azi_ref
        self.muller = muller
    
    def reflection_brdf(self,dir_in, dir_out, surf_norm, wl):
        R = np.zeros((4,4))
        T = np.zeros((4,4))
        lambertian_coeff = np.max([np.dot(surf_norm,dir_out) / np.pi,0])
        
        zen_in = np.arccos(np.dot(surf_norm,-dir_in)) #NOTE: -dir_in points away from the boundary
        zen_out = np.arccos(np.dot(surf_norm,dir_out))
        if zen_out > np.pi / 2: # in this case the outgoing ray is transmitted, not reflected
            return(create_wl_stack(R,wl.size),create_wl_stack(T,wl.size),dir_out)
        azi_out = _get_azi_angle(dir_in, dir_out, surf_norm)
        
        intfun_zen_in = interp1d(self.zen_inc, self.muller, axis=0)
        temp_muller = intfun_zen_in(zen_in)
        
        intfun_wl = interp1d(self.wavelength,temp_muller,axis=0,fill_value='extrapolate')
        temp_muller = intfun_wl(wl)
        
        intfun_zen_out = interp1d(self.zen_ref,temp_muller,axis=1)
        temp_muller = intfun_zen_out(zen_out)
        
        intfun_azi_out = interp1d(self.azi_ref,temp_muller,axis=1)
        R = intfun_azi_out(azi_out) * lambertian_coeff
        
        T = np.zeros((wl.size,4,4))
        
        return(R, T, dir_out)
        
    
class SpatVarLamb:
    def __init__(self,wavelength,reflectances):
        self.wavelength = wavelength
        self.reflectances = reflectances
        self.current_pixel = 0 
        # NOTE: in this case the refl_fun interface doesn't take position into
        # account. This is a problem with the varying surface reflectance. A
        # quick fix is to run the code sequentially and let the function to 
        # store a 'state'. This is very ugly and will cause problems later on,
        # but necessary now.
    def reflection_lambertian_spatial(self,dir_in, dir_out, surf_norm, wl):
        #shape_size = 15
        shape_size = self.reflectances.shape[1]
        j = self.current_pixel % shape_size
        i = self.current_pixel // shape_size
        # This might seem wrong, but this is in line with the current ARSCA
        # camera_fov function.
        self.current_pixel += 1
        alb_interp_fun = interp1d(self.wavelength,self.reflectances[i,j,:],fill_value='extrapolate')
        albedo = alb_interp_fun(wl)
        lambertian_coeff = albedo * np.max([np.dot(surf_norm,dir_out) / np.pi,0])
        R_diff = np.zeros((4,4))
        R_diff[0,0] = 0.5
        R_diff[1,1] = 0.5
        R_diff[0,1] = 0.5
        R_diff[1,0] = 0.5
        R = create_wl_stack(R_diff,wl.size)
        R = lambertian_coeff[:, np.newaxis, np.newaxis] * R
        T = np.zeros((wl.size,4,4))
        return (R,T,dir_out)
    
class Reflection:
    refl_type = 0
    refl_param = []
    input_file = ''
    
    def __init__(self,refl_type,refl_param,input_file=''):
        self.refl_type = refl_type
        self.refl_param = refl_param
        self.input_file = input_file
        self.set_up_reffun()
        
    def set_up_reffun(self):
        if self.refl_type == 0:
            # Pass-through boundary
            self.reflection = self.reflection_passthrough
        elif self.refl_type == 1:
            # Lambertian reflection
            self.reflection = self.reflection_lambertian
        elif self.refl_type == 2:
            # Lambertian with spatial variance
            if self.input_file == '':
                print("No input file for spatially varying Lambertian specified but it is required!")
            with netCDF4.Dataset(self.input_file) as ds:
                b = {}
                for key in ds.variables.keys():
                    b[key] = ds[key][:].data
            self.spatvarlamb = SpatVarLamb(b['wavelengths'], b['reflectance'])
            self.reflection = self.spatvarlamb.reflection_lambertian_spatial
        elif self.refl_type == 3:
            # BRDF reflection
            if self.input_file == '':
                print("No input file specified but BRDF requested!")
            with netCDF4.Dataset(self.input_file) as ds:
                b = {}
                for key in ds.variables.keys():
                    b[key] = ds[key][:].data
            self.brdf = BRDF(b['zen_inc'],b['wavelength'],b['zen_ref'],b['azi_ref'],b['muller'])
            self.reflection = self.brdf.reflection_brdf
        else:
            raise ValueError("No reflection function '%d' defined in raysca!" % self.refl_type)
        
    def reflection_passthrough(self,dir_in,dir_out,surf_norm,wl):
        if not np.all(np.isclose(dir_in,dir_out)):
            dir_out = dir_in
        # The pass-through surface doesn't reflect anything, but attenuates whatever goes
        # through it.
        T = (1.0 - self.refl_param) * np.eye(4)
        T = create_wl_stack(T,wl.size)
        R = np.zeros((wl.size,4,4))
        return (R, T, dir_out)
    
    def reflection_lambertian(self,dir_in, dir_out, surf_norm,wl): 
        lambertian_coeff = self.refl_param * np.max([np.dot(surf_norm,dir_out) / np.pi,0])
        R_diff = np.zeros((4,4))
        R_diff[0,0] = 0.5
        R_diff[1,1] = 0.5
        R_diff[0,1] = 0.5
        R_diff[1,0] = 0.5
        R = lambertian_coeff * R_diff
        T = np.zeros((wl.size,4,4))
        R = create_wl_stack(R,wl.size)
        return(R, T, dir_out)
        
def create_wl_stack(a,n):
    a_ = a[np.newaxis,:,:]
    return np.repeat(a_,n,axis=0)
    
if False:
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


def set_up_reflection(boundary,phase_folder,brdf_file):
    """
    Sets up reflection functions for each of the boundary surfaces. This could
    include setting up look-up-tables from files for example.

    The reflection functions do the sampling (if applicable) of the reflection
    direction and the generation of reflection/transmission Muller matrices.

    Returns a list of all reflection functions used in the computation.
    """
    refl_fun_idxs = boundary['reflection_kernel'][:]
    refl_fun_params = boundary['reflection_kernel_parameter'][:,:]
    reflections = []

    input_file = phase_folder + brdf_file

    for idx,choice in enumerate(refl_fun_idxs):
        reflections.append(Reflection(choice,refl_fun_params[:,idx],input_file))

    return reflections


