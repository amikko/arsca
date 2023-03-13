#Radiative transfer simulator wrapper for ARSCA
#wrapper for launching the simulators and reading the spectral data afterwards
#also for the configuration of the simulators themselves

import subprocess
import time
import os
import numpy as np
import sys
import netCDF4

path_to_folder_containing_this_file = os.path.dirname(os.path.realpath(__file__))

case = 'foo'
configuration = 'bar'

siro_folder = path_to_folder_containing_this_file + "/rt_solvers/siro/"
arsca_solver_folder = path_to_folder_containing_this_file + "/rt_solvers/arsca-solver/"
raysca_folder = path_to_folder_containing_this_file + "/rt_solvers/raysca/"

relative_root_folder = "../../"

def run(solver,input_file,extra_params=[]):
    """
    Runs solver simulator using input_file.
    """
    solvers = ['siro','arsca-solver','raysca','libradtran']
    if solver not in solvers:
        raise KeyError("Chosen solver '%s' not recognized! \nAvailable solvers: %s." % (solver,solvers))
    if solver == 'siro':
        outfilename = configuration
        launch_string = "./siro {} {}".format(relative_root_folder + input_file,outfilename)
        os.chdir(siro_folder)
        siro_proc = subprocess.Popen(launch_string,shell=True)
        siro_proc.communicate()
        if siro_proc.returncode != 0:
            raise RuntimeError("Siro crashed!")
        os.chdir(relative_root_folder)
        return load_siro_result_file('radiance' + outfilename)
    elif solver == 'arsca-solver':
        sys.path.append(arsca_solver_folder)
        import arsca_solver
        import sounding
        os.chdir(arsca_solver_folder)
        sndg, xsecs = sounding.sounding_from_nc(relative_root_folder + input_file)
        flags = {'interface' : True}
        spec = arsca_solver.rt_simu(sndg, flags, xsecs) #spectrum gets unpacked here
        os.chdir(relative_root_folder)
        return spec
    elif solver == 'raysca':
        #sys.path.append(raysca_folder)
        #os.chdir(raysca_folder)
        #print(os.getcwd())
        #import raysca
        #Below: this would be more ideal way?
        import rt_solvers.raysca.raysca as raysca
        raysca.set_paths(raysca.__file__)
        #stokes = raysca.run_simulation(relative_root_folder + input_file)
        stokes = raysca.run_simulation(input_file)
        #os.chdir(relative_root_folder)
        return stokes
    elif solver == 'libradtran':
        run_libradtran(input_file,extra_params)
        return read_libradtran_outfile(extra_params['rte_solver'])

def load_siro_result_file(dfilename_in):
    #This function reads the result file and returns the simulated transmittance
    #spectrum. This function hails from the ALTIUS era, so nomenclature might be
    #a bit old.
    ddirname = siro_folder + "output/"
    if not os.path.isdir(ddirname):
        print("The specified case folder doesn't exist. Exiting...")
        return 0

    siro_settings_filename = siro_folder + "siro_settings.nml"
    with open(siro_settings_filename) as sirosetfile:
        settings_full = sirosetfile.read()
        settings = settings_full.split('\n')
        for setting in settings:
            if 'usepolar' in setting:
                polarization_enabled = 'true' in setting
                break

    if polarization_enabled:
        suffixes = ['.psiroI','.psiroQ','.psiroU','.psiroV']
    else:
        suffixes = ['.siro']

    unformatted_radiances = []
    for suffix in suffixes:
        """
        files = ([f for f in os.listdir(ddirname)
        if os.path.isfile(os.path.join(ddirname, f))])
        for f in files:
            if f[-len(suffix):] == suffix:
                dfilename = f
        try:
            dfilepath = ddirname + dfilename
        except UnboundLocalError:
            print("Data not found. Most likely an erronous run of Siro.")
            return 0
            """
        dfilepath = ddirname + dfilename_in + suffix
        radiance_data = np.genfromtxt(dfilepath)
        param_data = np.genfromtxt(dfilepath[:-len(suffix)] + ".param")
        noph = param_data[1]
        geoms = np.unique(radiance_data[:,0])
        wls = np.unique(radiance_data[:,1])
        wlAmt = len(wls)
        altAmt = len(geoms)
        ordAmt = radiance_data.shape[1] - 2
        result_data = np.zeros((wlAmt,altAmt,ordAmt))

        rownum = 0
        for idxAlt in range(altAmt):
            for idxWl in range(wlAmt):
                result_data[idxWl,idxAlt,:] = radiance_data[rownum,2:]
                rownum = rownum + 1
        #now siro is tuned so that it outputs the real radiances, not photon counts
        #result_data = result_data / noph
        unformatted_radiances.append(result_data)

    if polarization_enabled:
        stokesAmt = 4
        radshape = unformatted_radiances[0].shape
        radiance_data = np.zeros(radshape + (stokesAmt,))
        for i in range(stokesAmt):
            radiance_data[:,:,:,i] = unformatted_radiances[i]
    else:
        radiance_data = unformatted_radiances[0]
    return radiance_data

def change_raysca_settings(setting,value):
    """
    Changes raysca settings in a settings file in usual arsca installation
    Both setting and value need to be in string format.
    """
    settings_path = path_to_folder_containing_this_file + "/rt_solvers/raysca/raysca_settings.yaml"
    row_idx = "" # This needs to be preset so it exists in this scope
    # and preferably it needs to be non-index so an error gets raised if
    # such setting isn't find (just an extra countermeasure for the for-else)
    with open(settings_path,'r') as f:
        lines = f.readlines()
        for l_idx, l in enumerate(lines):
            ls = l.strip()
            lsl = ls.split(':')
            if lsl[0].strip() == setting:
                row_idx = l_idx
                break
        else:
            raise ValueError('Could not find setting "%s" in raysca_settings.yaml!' % setting)
    with open(settings_path,'w') as f:
        lines[row_idx] = "%s : %s\n" % (setting,value)
        f.writelines(lines)

def create_siro_settings(custom_settings={}):
    """
    Creates the settings .nml file for Siro.

    The parameters and their default values for the Siro simulation kernel are:
        noph : The number of photons simulated for a single wavelength and
        and a single geometry.
            Default: 50000

        usepolar : If set to true, the full Stokes vector will be simulated. If
        not, then only the unpolarized component is used.
            Default: True

        userefrac : If set to true, refraction will be used in optical path
        calculations.
            Default: False

        atmos_layers : Into how many layers is the atmosphere discretized
        internally. This can be different than the atmospheric layer count in
        the input.
            Default: 10001

        maxnoord : Maximum number of scattering orders distincted in the output
        file. If set to zero, then all the transmittances are summed together.
        Otherwise, each of the orders of scattering up to maxnoord are saved in
        separate output files.
            Default: 3

        maxnolay : The line of sight discretization for scattering point
        calculations.
            Default: 501

        minweight : The minimum weight of a photon packet. If a photon weight
        falls under this value, it is eliminated. Photon weight starts from 1.0
        and it is reduced through absorption and scattering.
            Default: 1.0e-6

        req : The radius of the Earth in kilometers.
            Default: 6371.0

        ratm : The radius of the top of the atmosphere in kilometers.
            Default: 6471.0

        step : The integration step length used in absorption and scattering
        calculations along a path through the atmoshere in kilometers.
            Default: 2.0
            
        BRF_FILENAME : The path to the netcdf4 file containing the surface BRDF
        model for the simulator.
            Default: "input/brdf/new_snow.nc4"
            
        AER_FILENAME : The path to the aerosol file used.
            Default: "input/miefiles/aer_ALTIUS.dat"
    """
    settings = {'noph' : 40000,
                'usepolar' : True,
                'userefrac' : False,
                'atmos_layers' : 1001,
                'maxnoord' : 6,
                'maxnolay' : 501,
                'minweight' : 1e-6,
                'req' : 6371.0,
                'ratm' : 6471.0,
                'step' : 2.0,
                'brdf_reflection' : False,
                'BRF_FILENAME' : "'input/brdf/new_snow.nc4'",
                'AER_FILENAME' : "'input/miefiles/aer_ALTIUS.dat'"}

    #replace the default settings with new values if desired
    default_settings_keys = list(settings.keys())
    custom_settings_keys = list(custom_settings.keys())
    while len(custom_settings_keys) > 0:
        key = custom_settings_keys[-1]
        if key in default_settings_keys:
            settings[key] = custom_settings[key]
            custom_settings_keys.remove(key)
        else:
            raise KeyError("Invalid settings key '%s'" % key)

    #transform the booleans into fortran format
    fortran_bool = lambda bool : '.true.' if bool else '.false.'
    for key in settings.keys():
        if type(settings[key]) == bool:
            settings[key] = fortran_bool(settings[key])

    #format and join the settings into a long string
    type_formats = {int : '%d',
                    float : '%1.5e',
                    str : '%s'}

    max_key_len = max(map(len,settings.keys()))

    settings_rows = []
    for key in settings.keys():
        value = settings[key]
        spacing = (max_key_len - len(key)) * " "
        formatting = type_formats[type(value)]
        settings_rows.append(("\t%s%s = %s" % (key, spacing, formatting)) % value)

    settings_file_rows = (["",
    " $param"] +
    settings_rows +
    [" $end\n"])

    nml_string = "\n".join(settings_file_rows)

    #finally write all that into a file
    settings_file_name = path_to_folder_containing_this_file + '/rt_solvers/siro/siro_settings.nml'
    with open(settings_file_name,'w') as nml_file:
        nml_file.writelines(nml_string)


lrt_basepath = path_to_folder_containing_this_file + "/rt_solvers/libradtran/"
lrt_temp_foldername = "arsca/"
lrt_rel_to_bin = "../" + lrt_temp_foldername
lrt_output_fname = 'arsca_output.out'

def find_lrt_folder(): 
    lrt_list = os.listdir(lrt_basepath)
    for name in lrt_list:
        if os.path.isdir(lrt_basepath + name) and 'libRadtran' in name:
            lrt_path = lrt_basepath + name + '/'
            return lrt_path
    raise ValueError("Couldn't find libRadtran installation!")

def run_libradtran(nc_file,libradtran_settings):
    lrt_path = find_lrt_folder()
    lrt_input_folder = lrt_path + lrt_temp_foldername
    definition_dict = {}
    if not os.path.isdir(lrt_input_folder):
        os.mkdir(lrt_input_folder)
    with netCDF4.Dataset(nc_file) as ds:
        # Save the wavelengths
        wl_file = 'wavelength.dat'
        wl = ds['source']['input_wavelength'][:].data
        np.savetxt(lrt_input_folder + wl_file,wl,header='ARSCA-generated wavelength input file') # TODO: Check if this format is okay!
        definition_dict['wavelength_grid_file'] = lrt_rel_to_bin + wl_file
        
        #Set up the atmosphere
        atmos_file = 'basic_atmos.dat'
        #ds['medium']['scatterer'][:,0]
        #TODO: create this one, libradtran manuel p. 81
        #Kuulema pitää säätää ideaalikaasuyhtälön avulla paineet sopiviksi, wtf
        definition_dict['atmosphere_file'] = lrt_rel_to_bin + atmos_file
        
        # Save the cross-sections and profiles
        med_pos = ds['medium']['position'][:].data #shape: n_medium_pos,3
        alts = np.linalg.norm(med_pos,axis=1)
        alts = alts - np.min(alts)
        definition_dict['mol_abs_param'] = 'crs'
        for i_gas,gas in enumerate(libradtran_settings['gases']):
            # The gas profile
            gas_filename = gas + '.dat'
            gas_arr = np.zeros((alts.size,2))
            gas_arr[:,0] = alts
            gas_arr[:,1] = ds['medium']['absorber'][:,i_gas].data
            gas_arr = gas_arr[::-1,:] # for some reason the spectra are upside down
            np.savetxt(lrt_input_folder + gas_filename,gas_arr)
            definition_dict['mol_file %s' % gas] = lrt_rel_to_bin + gas_filename
            
            # The gas absorption cross-section
            # TODO: Come up a way how to handle the varying cross-sections in each altitude
            
            crs_filename = gas + '_crs.dat'
            crs_arr = np.zeros((wl.size,2))
            crs_arr[:,0] = wl
            crs_arr[:,1] = ds['medium']['absorbing_cross_section'][0,:,i_gas].data
            np.savetxt(lrt_input_folder + crs_filename,crs_arr)
            definition_dict['crs_file %s' % gas] = lrt_rel_to_bin + crs_filename
        
        
        #Set up the spectrum file
        spec_file = 'solar_spectrum.dat'
        spec_arr = np.zeros((wl.size,2))
        spec_arr[:,0] = wl
        spec_arr[:,1] = ds['source']['incident_stokes'][0,:,0]
        np.savetxt(lrt_input_folder + spec_file, spec_arr)
        definition_dict['source solar'] = lrt_rel_to_bin + spec_file
        
        #Set the albedo
        albedo = ds['boundary']['reflection_kernel_parameter'][0][0]
        definition_dict['albedo'] = "%f" % albedo
        
        #Set the solar zenith angle
        inc_dir = ds['source']['incident_direction'][0,:].data
        surf_norm = med_pos[1,:] - med_pos[0,:]
        surf_norm = surf_norm / np.linalg.norm(surf_norm)
        sza = 180.0 / np.pi * np.arccos(np.dot(surf_norm, -inc_dir))
        definition_dict['sza'] = "%f" % sza
        
        #Set the viewing geometry
        umu = np.dot(surf_norm,-ds['instrument']['view_vector'][0,:].data)
        print("Only nadir available!")
        definition_dict['umu'] = "%f" % umu
        #definition_dict['zout'] = "TOA"
        
        #definition_dict['phi'] = "%f" % 0
        #definition_dict['phi0'] = "%f" % 0
        
        
    # Create the problem definition file
    input_fname = 'arsca_input.inp'
    skip_keys = ['gases']
    with open(lrt_input_folder + input_fname,'w') as fhandle:
        for key in definition_dict:
            fhandle.write("%s %s\n" % (key, definition_dict[key]))
        for key in libradtran_settings:
            if key not in skip_keys:
                fhandle.write("%s %s\n" % (key, libradtran_settings[key]))
        
    # Change to bin folder and run the simulation
    original_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(lrt_path + 'bin/')
    os.system('./uvspec < %s > %s' % (lrt_rel_to_bin + input_fname, lrt_rel_to_bin + lrt_output_fname))
    os.chdir(original_path)
    
def read_libradtran_outfile(solver):
    lrt_path = find_lrt_folder()
    if solver == 'disort':
        lrt_out_ = np.genfromtxt(lrt_path + lrt_temp_foldername + lrt_output_fname,usecols=(0,1))
        lrt_out = lrt_out_[1::2,1]
        return lrt_out
    elif solver == 'mystic':
        lrt_out = np.genfromtxt(lrt_path + lrt_temp_foldername + lrt_output_fname)
        return lrt_out[:,1] # TODO: Check this!
    #elif solver == 'polradtran':
    #    return lrt_out[:,1::2]
    
libradtran_settings = {
        'gases' : ['O2'],
        'rte_solver' : 'disort',
        #'rte_solver' : 'polradtran',
        #'polradtran nstokes' : 4
        }
if __name__ == "__main__":
    run_libradtran('./solver_input/snow_cover_study_o2a/siro_brdf_o2a.nc',libradtran_settings)
