#Data transformation module for ARSCA

#handle creation of nc-idf compatible input data from various different sources
#handle simpler tasks, such as wn2wl, etc...

# TODO: Add SSA-extinction to abs-sca basis

# TODO: Add raysca source-LoS cartesian product handling

import numpy as np

# --------------------------
# Geometry generation stuffs
# --------------------------

def limb_geometries(sat_altitude,tan_altitude,solar_zenith,solar_azimuth):
    """
    Generates limb viewing geometries for satellite at altitude sat_altitude.
    Returns len(tan_altitude) * len(solar_zenith) geometries.
    """
    #TODO: This hardcoded 6371 might cause problems some time
    tan_surf_point = np.array([0.0,0.0,6371.0])
    view_dir = np.array([1.0,0.0,0.0])
    sat_polar_vec = np.array([0.0,1.0,0.0])
    tan_surf_norm = norm_vec(tan_surf_point)
    n_geometry = tan_altitude.size * solar_zenith.size
    n_coordinate = 3
    sat_pos = np.zeros((n_geometry,n_coordinate))
    sat_view_vec = np.zeros((n_geometry,n_coordinate))
    sat_polar = np.zeros((n_geometry,n_coordinate))
    solar_dir = np.zeros((n_geometry,n_coordinate))
    solar_polar = np.zeros((n_geometry,n_coordinate))
    for i in range(tan_altitude.size):
        for j in range(solar_zenith.size):
            geom_idx = i * solar_zenith.size + j
            sat_pos[geom_idx,:] = sat_position(sat_altitude,tan_altitude[i],tan_surf_point,view_dir)
            sat_view_vec[geom_idx,:] = view_dir
            sat_polar[geom_idx,:] = sat_polar_vec
            solar_dir[geom_idx,:] = solar_direction(tan_surf_point,view_dir,solar_zenith[j],solar_azimuth[j])
            solar_polar[geom_idx,:] = perp_norm_vec(solar_dir[geom_idx,:],sat_polar_vec) #NOTE: This is actually probably not right.
    return (sat_pos,sat_view_vec,sat_polar,solar_dir,solar_polar)

def sat_position(sat_alt,tan_alt,tan_surf_point,view_dir):
    """
    The satellite position is in a point where the line with a point
    tan_surf_point + norm_vec(tan_surf_point) * tan_alt
    along a vector
    view_dir
    intersects a sphere with radius
    np.norm(tan_surf_point) + sat_alt
    """
    R_pos = np.linalg.norm(tan_surf_point) + sat_alt
    p_tan = tan_surf_point + norm_vec(tan_surf_point) * tan_alt
    a = np.dot(view_dir,view_dir)
    b = 2 * np.dot(view_dir,p_tan)
    c = np.dot(p_tan,p_tan) - R_pos**2
    s = (- b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    return p_tan + s * view_dir

def direct_geometries(sat_alt,solar_zeniths):
    surf_vec = np.array([0.0,0.0,sat_alt])
    sat_pos_base = surf_vec
    view_dir_base = np.array([1.0,0.0,0.0])
    solar_azimuth = 0.0
    n_geometry = len(solar_zeniths)
    n_coordinate = 3
    sat_pos = np.zeros((n_geometry,n_coordinate))
    sat_view_vec = np.zeros((n_geometry,n_coordinate))
    incident_dir = np.zeros((n_geometry,n_coordinate))
    
    for zen_idx,zen_angle in enumerate(solar_zeniths):
        incident_dir[zen_idx,:] = -solar_direction(surf_vec,view_dir_base,zen_angle,solar_azimuth)
        sat_pos[zen_idx,:] = sat_pos_base
        sat_view_vec[zen_idx,:] = -incident_dir[zen_idx,:]
    
    return (sat_pos, sat_view_vec, incident_dir)
# -------------------------
# Solar stuffs
# -------------------------

def solar_direction(surf_vec,view_dir,solar_zenith,solar_azimuth):
    """
    Input: The point on the surface (surf_vec) from where the zenith and azimuth
    vectors are derived.

    surf_norm points directly to the zenith = 0. view_dir points directly
    to azimuth = 0.
    
    Zenith and azimuth angles are in degrees.
    """
    surf_norm = norm_vec(surf_vec)
    d2r = np.pi / 180.0 #degrees to radians
    #first rotate around (surf_norm X view_dir) with the zenith angle
    #and then rotate around (surf_norm) with the azimuth angle
    zen_rot_ax = perp_norm_vec(surf_norm,view_dir)
    solar_dir = surf_norm #we start with this
    solar_dir = arb_rotation(solar_dir,d2r * solar_zenith,zen_rot_ax)
    solar_dir = arb_rotation(solar_dir,d2r * solar_azimuth,surf_norm)
    return solar_dir

# TODO: a function (solar_dir, view_vec, surf_vec) -> (solar_zenith, solar_azimuth)

# -------------------------
# Spectrum stuff
# -------------------------

def wn2wl(wn):
    """
    Transforms a wavenumber array to a wavelength array.

    Wavelengths are in nanometres and wavenumbers are 1/cm.
    """
    conversion_fun = lambda x : 1e7 / x
    if type(wn) == type(np.array([])):
        return conversion_fun(wn[::-1])
    elif type(wn) == type([]):
        return list(map(conversion_fun,wn[::-1]))
    else:
        #this is supposed to be a single number
        return conversion_fun(wn)

def wl2wn(wl):
    """
    Transforms a wavelength array to a wavenumber array.

    Since the function is symmetrical, this function is just a wrapper.
    """
    return wn2wl(wl)

def raysca_stokes2scalar(spec,retain_shape=True):
    """
    RaySca outputs spectra in shape (n_source, n_los, n_wavelength, n_stokes)
    
    This function will flatten the stokes into a single value. 
    
    If retain_shape is True, then output will be of shape 
    (n_source, n_los, n_wavelength). Else, the output shape will be 
    (n_source * n_los * n_wavelength,). 

    """
    
    if retain_shape:
        return spec.sum(axis=3)
    else:
        return spec.sum(axis=3).ravel()


def column_average(alt_grid, gas_prof, atmos_file, day_idx):
    """
    Computes the column average of a gas
    
    Note that this function needs to be modified to compute 'the proper' column
    average as this one just uses some h2o, gravity and pressure profiles. They 
    would need to be supplied.
    
    Potential source: https://amt.copernicus.org/articles/3/947/2010/amt-3-947-2010.pdf
    """
    import netCDF4
    from scipy.interpolate import interp1d
    avo = 6.02214129e23
    m_h2o = 18.02/1e3/avo
    m_dryair = 28.964/1e3/avo
    with netCDF4.Dataset(atmos_file) as ds:
        center_alts = 0.5 * (alt_grid[1:] + alt_grid[:-1])
        gas_cent_intpfun = interp1d(alt_grid,gas_prof)
        gas_cent = gas_cent_intpfun(center_alts)
        grav_cent_intpfun = interp1d(alt_grid,ds['prior_gravity'][day_idx,:])
        grav_cent = grav_cent_intpfun(center_alts)
        h2o_cent_intpfun = interp1d(alt_grid,ds['prior_h2o'][day_idx,:])
        h2o_cent = h2o_cent_intpfun(center_alts)
        #pres_intpfun = interp1d(alt_grid,ds['prior_Pressure'])
        #h2o_cent = h2o_cent_intpfun(center_alts)
        try:
            dP = abs(ds['prior_Pressure'][day_idx,:-1] - ds['prior_Pressure'][day_idx,1:])
        except IndexError:
            dP = abs(ds['prior_pressure'][day_idx,:-1] - ds['prior_pressure'][day_idx,1:])
    h2o_dry = h2o_cent / (1 - h2o_cent)
    h = dP / (m_dryair * grav_cent * (1 + h2o_dry * (m_h2o / m_dryair)))
    h = h / np.sum(h)
    xgas = np.dot(h,gas_cent / (1 - h2o_cent))
    return xgas

def dummy_column_average(alt_grid, gas_prof):
    return column_average(alt_grid, gas_prof, './files_used_in_examples/so20090516_20181030.public.nc', 800)
    
# -------------------------
# Some math stuffs
# -------------------------

def perp_norm_vec(v1,v2):
    cross = np.cross(v1,v2)
    tol = 1e-8
    if np.linalg.norm(cross) < tol:

        cross_cand1 = np.cross(v1,np.array([0.0, 0.0, 1.0]))
        cross_cand2 = np.cross(v1,np.array([0.0, 1.0, 0.0]))
        if np.linalg.norm(cross_cand1) < tol:
            cross = cross_cand2
        else:
            cross = cross_cand1
    return norm_vec(cross)

def norm_vec(v):
    return v / np.linalg.norm(v)

def arb_rotation(v,ang,ax):
    """
    Rotates the vector v around the vector ax ang much.
    """
    ax = norm_vec(ax) #normalized for good measure

    slow_and_enormous_rotation_matrix = False
    if slow_and_enormous_rotation_matrix:
        #this is taken from Worlfram MathWorld page
        #http://mathworld.wolfram.com/RodriguesRotationFormula.html
        cosa = np.cos(ang)
        mcosa = 1 - cosa
        sina = np.sin(ang)
        ax0 = ax[0]
        ax1 = ax[1]
        ax2 = ax[2]
        M = np.array([
        [cosa + ax0**2*mcosa, ax0*ax1*mcosa - ax2*sina, ax0*ax2*mcosa + ax1*sina],
        [ax1*ax0*mcosa + ax2*sina, cosa + ax1**2*mcosa, ax1*ax2*mcosa - ax0*sina],
        [ax2*ax0*mcosa - ax1*sina, ax2*ax1*mcosa + ax0*sina, cosa + ax2**2*mcosa]
        ])
        return M@v
    else:
        #this is a simplified version from the Rodrigues Rotation Wikipedia page:
        return (v * np.cos(ang) + np.cross(ax,v) * np.sin(ang)
         + ax * np.dot(ax,v) * (1 - np.cos(ang)))
