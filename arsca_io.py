#The IO-module for ARSCA

#handle various input file reading
#handle nc-idf creation

import os
import time
import uuid
import numpy as np
import netCDF4

from arsca_utilities import unique_filename, if_folder_not_exist_then_create

input_folder = "./solver_input/"

case = 'foo'
configuration = 'bar'

float_acc = 'f8'
int_acc = 'i4'

def read_atmos(path,gases,atmos_choice=0):
    """
    Reads atmosphere data of format so20140319.mav, or whatever that is
    """

    with open(path, "r", encoding="utf-8") as file_handle:
        atmos_data = {'Height' : [], 'Temp' : [],
                          'Pres' : [], 'Density' : []}
        for gas in gases:
            atmos_data[gas] = []

        #in files with multiple atmoses, let's pick the selected one
        nexts = 0
        for line in file_handle:
            split_line = line.split()
            if split_line[0] == 'Next':
                if nexts == atmos_choice:
                    break
                else:
                    nexts = nexts + 1
        else:
            print('ERROR: atmos_choice invalid.')
            return

        datalines = []
        data_reached = False
        #let's find the line containing the column names
        for line in file_handle:
            split_line = line.split()
            if split_line[0] == 'Height':
                col_names = split_line
                data_reached = True
            elif data_reached:
                if split_line[0] == 'Next':
                    break
                datalines.append(split_line)

        columns = {}
        #then we'll find indexes of the columns of interest
        for i in range(len(col_names)):
            if col_names[i] in atmos_data:
                columns[i] = col_names[i]

        #read the data to the atmos_data
        for dataline in datalines:
            for col_idx in columns:
                atmos_data[columns[col_idx]].append(float(dataline[col_idx]))

    for key in atmos_data:
        atmos_data[key] = np.array(atmos_data[key])

    # Some files may have weird negative altitudes. Let's remove those
    mask = atmos_data['Height'] >= 0
    for key in atmos_data:
        atmos_data[key] = atmos_data[key][mask]

    return atmos_data


def read_nc_atmos(nc_file,atmos_ind,formatting):
    atmosdict = {}
    #co2 from ppm to fraction
    #ch4 from ppb to fraction
    #pressure from hPa to atm
    
    if formatting == 'tccon':
        coeff_dict = {'co2' : 1e-6, 'ch4' : 1e-9, 'Pressure' : 0.000987}
    elif formatting == 'ggg':
        coeff_dict = {}
    with netCDF4.Dataset(nc_file) as ds:
        for k in ds.variables.keys():
            for coeff in coeff_dict:
                if coeff in k:
                    c = coeff_dict[coeff]
                    break
            else:
                c = 1
            #6: skips the 'prior_' prefix of the keys
            try:
                if 'prior' in k:
                    atmosdict[k[6:]] = ds[k][atmos_ind,:].data * c
                else:
                    atmosdict[k] = ds[k][atmos_ind,:].data * c
            except ValueError:
                atmosdict[k] = ds[k][:].data * c
    return atmosdict

def read_afgl_atmos(filename):
    atmosdict = {}
    atmosarr = np.genfromtxt(filename,skip_header=2)
    atmosarr = atmosarr[::-1,:] # the file starts with the highest altitude, so
    #we'll flip it
    atmosdict['Height'] = atmosarr[:,0]
    atmosdict['Pres'] = atmosarr[:,1] / 1013.25 # convert millibars to atmos
    atmosdict['Temp'] = atmosarr[:,2]
    atmosdict['Density'] = atmosarr[:,3]
    atmosdict['o3'] = atmosarr[:,4] / atmosdict['Density']
    atmosdict['o2'] = atmosarr[:,5] / atmosdict['Density']
    atmosdict['h2o'] = atmosarr[:,6] / atmosdict['Density']
    atmosdict['co2'] = atmosarr[:,7] / atmosdict['Density']
    atmosdict['no2'] = atmosarr[:,8] / atmosdict['Density']
    return atmosdict
#
#The Intermediate NetCDF File Format

#TODO: Phase function NetCDF Files
#TODO: BRDF NetCDF files

"""
Unit conventions:
    Distances, altitudes and positions are all in kilometers.
    Temperatures are in Kelvin.
    Pressures are in
    Wavelengths are all in nanometers.
    Densities are in 1 / cm^3.
    Cross-sections are in cm^2.
    Stokes vectors are IQUV.
    Source angular radius is in degrees.
    Muller matrices are normalized so that TODO
    Theta angles are zero toward the propagation direction and pi opposite to
        the propagation direction. Theta is in [0,pi].
    Phi angles are zero toward the first polarization direction, pi/2 toward the
        second polarization direction. Phi is in [0,2*pi]

"""

#TODO: n_parameter dimension to other stuff too

#TODO: Add checks for input parameters
#TODO: Fill all missing values into the dicts

#TODO: 
def create_simulator_input(medium,instrument,source,
                            boundary,overwrite=True):
    """
    This function defines and creates the NetCDF4 file which acts as the input
    for various simulator kernels, such as ARSCA or Siro.

    The dimensions in the generated NetCDF4 file are as follows:
        n_medium_position : With how many basis functions are used in describing the
        (atmospheric) medium.

        n_input_wavelength : How many wavelengths are used in defining
        cross-sections and refractive indices.

        n_output_wavelength : How many wavelengths are used in the
        simulations. The distinction with input wavelengths is to enable
        different wavelength discretizations for simulations using the same
        cross-sections.

        n_absorber : How many absorbing substances are there in the atmosphere.

        n_scatterer : How many scattering substances are there in the
        atmosphere.

        n_parameter : How many parameters a scattering or reflecting kernel
        uses. This dimension is unlimited.

        n_emitter : How many emitting substances are there in the atmosphere.

        n_los : How many line-of-sight geometries are simulated with
        these parameters. For example, in limb measurements you could have
        several tangent altitudes with the same atmosphere.

        n_source : How many sources are used in the simulations.

        n_coordinate : How many dimensions are in all position or direction
        vectors. This is 3.

        n_stokes : How many Stokes parameters are used in polarized
        transmittances. This is 4.

    The structure of the file is as follows:
    root_group
        medium (Group)
            position (n_medium_position, n_coordinate) ('f8')
            absorber (n_medium_position, n_absorber) ('f8')
            scatterer (n_medium_position, n_scatterer) ('f8')
            emitter (n_medium_position, n_emitter) ('f8')
            absorbing_cross_section (n_medium_position, n_input_wavelength, n_absorber) ('f8')
            scattering_cross_section (n_medium_position, n_input_wavelength, n_scatterer) ('f8')
            medium_emissivity (n_medium_position, n_input_wavelength, n_stokes, n_emitter) ('f8')
            refractive_index (n_medium_position, n_input_wavelength) ('f8')
            scatterer_kernel (n_scatterer) ('i4')
            scatterer_kernel_parameter (n_parameter, n_scatterer) ('f8')
            interpolation_function (n_medium_positions) ('i4')
            interpolation_parameter (n_medium_positions) ('i4')
        instrument (Group)
            position (n_los, n_coordinate) ('f8')
            view_vector (n_los, n_coordinate) ('f8')
        source (Group)
            source_type (n_source) ('i4')
            incident_direction (n_source, n_coordinate) ('f8')
            source_position (n_source, n_coordinate) ('f8')
            input_wavelength (n_input_wavelength) ('f8')
            output_wavelength (n_output_wavelength) ('f8')
            incident_stokes (n_source, n_input_wavelength, n_stokes) ('f8')
            source_angular_radius (n_source) ('f8')
        boundary (Group)
            shape (n_boundary) ('i4')
            parameter (n_boundary) ('f8')
            reflection_kernel (n_boundary) ('i4')
            reflection_kernel_parameter (n_parameter, n_boundary) ('f8')

    In-depth description is found in the variable attributes defined below.
    """


    input_case_folder = input_folder + case + '/'
    if_folder_not_exist_then_create(input_case_folder)
    input_file_path = unique_filename(input_case_folder,configuration,'.nc',overwrite)
    
    # TODO: Finish this function definition and enable the following line
    #(medium, instrument, source, boundary) = validate_and_fill_in_input(medium, instrument, source, boundary)
    
    with netCDF4.Dataset(input_file_path, 'w', format="NETCDF4") as ds:
        #create the groups and description
        ds.description = "An intermediate radiative transfer simulation data input file."
        ds.history = "Created %s." % time.ctime(time.time())
        grp_medium = ds.createGroup('/medium')
        grp_medium.description = """
            Contains the information about the medium, such as the
            absorption, scattering and emission of the medium."""
        grp_instrument = ds.createGroup('/instrument')
        grp_instrument.description = """
            Contains the information about the measurement instrument (observer).
            For RT simulations, the necessary information is the position of the instrument
            and its viewing direction."""
        grp_source = ds.createGroup('/source')
        grp_source.description = """
            Contains the information about the radiation sources. A source may have
            several type"""
        grp_boundary = ds.createGroup('/boundary')
        grp_boundary.description = """
            Contains the information about the bondary of the computation domain.
            Only two types of bondaries are currently implemented: a sphere and a
            plane. Their reflection functions are also described here."""
    
        #creating the dimensions
        #TODO: Assert that all the input dimensions match
        #medium['position'].shape[0] > medium['absorbing_cross_section'].shape[0]
        n_medium_position = medium['position'].shape[0]
        #n_medium_ medium['absorbing_cross_section'].shape[0]
        n_absorber = medium['absorber'].shape[1]
        n_scatterer = medium['scatterer'].shape[1]
        n_emitter = medium['emitter'].shape[1]
        n_input_wavelength = source['input_wavelength'].size
        n_output_wavelength = source['output_wavelength'].size
        n_source = source['type'].size
        n_boundary = boundary['shape'].size
        n_los = instrument['view_vector'].shape[0]
    
        ds.createDimension('n_medium_position', n_medium_position)
        ds.createDimension('n_absorber', n_absorber)
        ds.createDimension('n_scatterer', n_scatterer)
        ds.createDimension('n_parameter', None)
        ds.createDimension('n_emitter', n_emitter)
        ds.createDimension('n_input_wavelength',n_input_wavelength)
        ds.createDimension('n_output_wavelength',n_output_wavelength)
        ds.createDimension('n_los', n_los)
        ds.createDimension('n_source',n_source)
        ds.createDimension('n_boundary', n_boundary)
        ds.createDimension('n_coordinate', 3)
        ds.createDimension('n_stokes', 4)
    
        #MEDIUM VARIABLE DEFINITIONS
        var_position = grp_medium.createVariable('position', float_acc, ('n_medium_position','n_coordinate'))
        var_position.description = """
            The medium is described using basis functions. The basis functions are
            defined in more detail in the interpolation_function description and
            ARSCA documentation."""
        var_absorber = grp_medium.createVariable('absorber', float_acc, ('n_medium_position', 'n_absorber'))
        var_absorber.description = """
            The absorbing part of the medium at each of the basis functions."""
        grp_medium.createVariable('scatterer', float_acc, ('n_medium_position', 'n_scatterer'))
        grp_medium.createVariable('emitter', float_acc, ('n_medium_position', 'n_emitter'))
        grp_medium.createVariable('refractive_index', float_acc, ('n_medium_position','n_input_wavelength'))
        grp_medium.createVariable('absorbing_cross_section', float_acc, ('n_medium_position', 'n_input_wavelength', 'n_absorber'))
        grp_medium.createVariable('scattering_cross_section', float_acc, ('n_medium_position', 'n_input_wavelength', 'n_scatterer'))
        grp_medium.createVariable('medium_emissivity', float_acc, ('n_medium_position', 'n_input_wavelength', 'n_stokes', 'n_emitter'))
        var_interp_fun = grp_medium.createVariable('interpolation_function', int_acc, ('n_medium_position',))
        var_interp_fun.description = """
            The interpolation function used to interpolate between points outside
            the positions of medium basis functions. The function selections are as
            follows:
                0: Interpolation is done linearly along the first coordinate of the
                points. If the interpolation point is outside the coordinate distribution,
                it is set to 0. Interpolation parameter has no effect.
    
                1: Interpolation is done linearly along the length of the position
                vectors. This stands for the usual spherical shell discretization
                of the medium. If the interpolation point is outside the coordinate distribution,
                it is set to 0. Interpolation parameter has no effect.
    
                2: Interpolation is radially from the medium position.
                There's no dependence between different medium positions. The basis
                function is a 3-dimensional Gaussian hump. The interpolation parameter
                is the standard deviation of the Gaussian."""
        grp_medium.createVariable('interpolation_parameter', float_acc, ('n_medium_position',))
    
        grp_medium.createVariable('scatterer_kernel', int_acc, ('n_scatterer',))
        grp_medium.createVariable('scatterer_kernel_parameter', float_acc, ('n_parameter','n_scatterer'))
    
        #INSTRUMENT VARIABLE DEFINITIONS
        grp_instrument.createVariable('position', float_acc, ('n_los', 'n_coordinate'))
        grp_instrument.createVariable('view_vector', float_acc, ('n_los', 'n_coordinate'))
    
        var_shape = grp_boundary.createVariable('shape', int_acc, ('n_boundary'))
        var_shape.description = """The shape of the boundary. The possible choices are:
            0: Flat surface. Parameter is the z-coordinate of the infinite plane.
            1: Spherical surface. Parameter is the radius of an origin-centered sphere."""
        grp_boundary.createVariable('parameter', float_acc, ('n_boundary'))
        var_ref_kernel = grp_boundary.createVariable('reflection_kernel', int_acc, ('n_boundary'))
        var_ref_kernel.description = """The reflection kernel used for this boundary.
        The kernel parameter is usually the albedo of the surface, if nothing is mentioned.
            0: Pass-through boundary. Usually the top of the atmosphere is this.
            10: Constant Lambertian. The most usual case.
            11: Spatially varying Lambertian in flat geometry. The variations should be defined in another file.
            20: BRDF. The reflectance function definitions are given in another file.
        """
        grp_boundary.createVariable('reflection_kernel_parameter', float_acc, ('n_parameter','n_boundary'))
    
        var_source_type = grp_source.createVariable('type', int_acc, ('n_source',))
        var_source_type.description = """
            0: far-field source. source/position disregarded.
            1: point source. incident_direction disregarded.
            2: cone source. source/parameter stands for the angular width of the cone."""
        grp_source.createVariable('incident_direction', float_acc, ('n_source', 'n_coordinate'))
        grp_source.createVariable('position', float_acc, ('n_source', 'n_coordinate'))
        grp_source.createVariable('input_wavelength', float_acc, ('n_input_wavelength',))
        grp_source.createVariable('output_wavelength', float_acc, ('n_output_wavelength',))
        grp_source.createVariable('incident_stokes', float_acc, ('n_source','n_output_wavelength', 'n_stokes'))
        grp_source.createVariable('parameter', float_acc, ('n_parameter','n_source'))
        var_source_ang_rad = grp_source.createVariable('source_angular_radius', float_acc, ('n_source',))
        var_source_ang_rad.description = """If the viewing beam points to the source with an angle
            below this value, it is considered looking at the source."""
    
        #populate the data
        grp_medium['position'][:,:] = medium['position']
        grp_medium['absorber'][:,:] = medium['absorber']
        grp_medium['scatterer'][:,:] = medium['scatterer']
        grp_medium['emitter'][:,:] = medium['emitter']
        grp_medium['refractive_index'][:,:] = medium['refractive_index']
        grp_medium['absorbing_cross_section'][:,:,:] = medium['absorbing_cross_section']
        grp_medium['scattering_cross_section'][:,:,:] = medium['scattering_cross_section']
        grp_medium['medium_emissivity'][:,:,:,:] = medium['medium_emissivity']
        grp_medium['interpolation_function'][:] = medium['interpolation_function']
        grp_medium['interpolation_parameter'][:] = medium['interpolation_parameter']
        grp_medium['scatterer_kernel'][:] = medium['scatterer_kernel']
        max_sca_ker_size = medium['scatterer_kernel_parameter'].shape[0]
        grp_medium['scatterer_kernel_parameter'][:max_sca_ker_size,:] = medium['scatterer_kernel_parameter']
    
        grp_instrument['position'][:,:] = instrument['position']
        grp_instrument['view_vector'][:,:] = instrument['view_vector']
    
        grp_boundary['shape'][:] = boundary['shape']
        grp_boundary['parameter'][:] = boundary['parameter']
        grp_boundary['reflection_kernel'][:] = boundary['reflection_kernel']
        max_ref_ker_size = boundary['reflection_kernel_parameter'].shape[0]
        grp_boundary['reflection_kernel_parameter'][:max_ref_ker_size,:] = boundary['reflection_kernel_parameter']
    
        grp_source['type'][:] = source['type']
        grp_source['incident_direction'][:,:] = source['incident_direction']
        grp_source['position'][:,:] = source['position']
        grp_source['input_wavelength'][:] = source['input_wavelength']
        grp_source['output_wavelength'][:] = source['output_wavelength']
        grp_source['incident_stokes'][:,:,:] = source['incident_stokes']
        max_src_par_size = source['parameter'].shape[0]
        grp_source['parameter'][:max_src_par_size,:] = source['parameter']
        grp_source['source_angular_radius'][:] = source['source_angular_radius']

    return input_file_path

def validate_and_fill_in_input(medium, instrument, source, boundary):
    """
    This function checks if the all necessary values are in there and if they're
    properly defined
    """
    # The two following dicts contain the required fields and other fields which
    # depend on them. The fields defined in required_fields are required in each
    # of the computations and the fields in optional_fields may be included in
    # the definition. If some input fields aren't defined, then an error is 
    # raised if they're required, but if they are optinal, they're just filled
    # in.
    required_fields = {'medium' : ['position','interpolation_function',
                                   'interpolation_parameter'],
                       'instrument' : ['position','view_vector'],
                       'boundary' : [],
                       'source' : ['input_wavelength', 'output_wavelength']}
    optional_fields = {'medium' : [['absorber', 'absorbing_cross_section'],
                                  ['scatterer', 'scattering_cross_section',
                                   'scatterer_kernel', 'scatterer_kernel_parameter'],
                                  ['emitter', 'medium_emissivity'],
                                  ['refractive_index']],
                      'instrument' : [[]],
                      'boundary' : [['shape', 'parameter', 'reflection_kernel', 
                                    'reflection_kernel_parameter']],
                      'source' : [['type', 'incident_direction', 'position', 
                                   'input_wavelength', 'output_wavelength', 
                                   'parameter', 'source_angular_radius']]}
    
    missing_fields_required = []
    missing_fields_optional = []
                  
    for key in required_fields.keys():
        for field in required_fields[key]:    
            if not eval(field + ' in ' + key):
                missing_fields_required.append(key + '[' + field +']')
    
    for key in optional_fields.keys():
        for field_group in optional_fields[key]:
            field_in_group = []
            for field in field_group:
                field_in_group.append(eval(field + ' in ' + key))
            if not all(field_in_group) and any(field_in_group):
                # in this some of the variables are defined
                missing_fields_optional.append(field_group)
                
    #Get the dimensions
    try:
        n_medium_position = medium['position'].shape[0]
    except KeyError:
        missing_fields_required.append('position')

    try:    
        n_absorber = medium['absorber'].shape[1]
    except KeyError:
        missing_fields_optional.append('absorber')

    try:    
        n_scatterer = medium['scatterer'].shape[1]
    except KeyError:
        missing_fields_optional.append('scatterer')

    try:    
        n_emitter = medium['emitter'].shape[1]
    except KeyError:
        missing_fields_optional.append('emitter')

    try:    
        n_input_wavelength = source['input_wavelength'].size
    except KeyError:
        missing_fields_required.append('input_wavelength')

    try:
        n_output_wavelength = source['output_wavelength'].size
    except KeyError:
        missing_fields_required.append('output_wavelength')
        
    try:
        n_source = source['type'].size
    except KeyError:
        missing_fields_required.append('source type')
        
    try: 
        n_boundary = boundary['shape'].size
    except KeyError:
        missing_fields_required.append('shape')
        
    try:
        n_los = instrument['view_vector'].shape[0]
    except KeyError:
        missing_fields_required.append('view_vector')
        
    if missing_fields_required:
        print("The following fields are missing:")
        print(missing_fields_required)
        
    if missing_fields_optional:
        print("Only some of the following optional fields are defined:")
        print(missing_fields_optional)
    
    return (medium, instrument, source, boundary)

def create_muller_matrices():
    """
    Creates the NetCDF4 file for Müller matrices for scatterers or reflectors.

    The dimensions in the generated NetCDF4 file are as follows:

        n_stokes : How many Stokes parameters are used in polarized
        transmittances. This is 4.

        n_wavelength : How many wavelengths are included in the scattering
        matrices. If you intend to have just one scattering kernel for the whole
        wavelength band, then set this to 1.

        n_theta : How many angles are in the theta discretization vector.

        n_phi : How many angles are in the phi discretization vector. If you
        intend to have no phi dependence, then set this to 1.
    root_group

    """
