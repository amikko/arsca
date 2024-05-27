#ARSCA: Atmospheric Radiation Simulation Computing Application

import arsca_io as io
import arsca_datatf as tf
import arsca_simulators as simu
import arsca_scattering as scatter
import arsca_reflectance as reflect
import arsca_cross_sections as xsec
import arsca_instrument as inst
import arsca_utilities as util
import arsca_emission as emi

"""
Cases and Configurations

A case is a broader set of configurations. Case is usually limited to a single
instrument. Configuration consists of some particular geometry, atmosphere and
solar angles.

These are just some words with which to organize the cross-section databases and
intermediate input files.

For example, case could be called "ALTIUS" and a configuration could be called
"only-rayleigh".

"""

case = "foo"
configuration = "bar"

def set_case(casename):
    case = casename
    io.case = casename
    simu.case = casename
    xsec.case = casename

def set_configuration(confname):
    configuration = confname
    io.configuration = confname
    simu.configuration = confname
    xsec.configuration = confname
