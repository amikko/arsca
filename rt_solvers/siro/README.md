# Backward Monte Carlo RT model Siro

First version published by Liisa Oikarinen (1965-2002) at Finnish Meteorological Institute in 1999 [1].

Further developed by Simo Tukiainen and Erkki Kyrölä.

This version has been modified by Antti Mikkonen and correspondence should be directed to antti.mikkonen@fmi.fi.

## Compilation

Siro should be compilable with a simple `make` call.

Siro requires `netcdf` (TODO: check version number) and it should be installed beforehand.

Installation can be easily done on Ubuntu with the command `sudo apt install libnetcdff-dev`

## Usage
In principle, Siro is ran with the following command: `./siro <input file name> <output file name>`

In practice, however, Siro should be ran through ARSCA with which you are able to create the input file and read the output file through Python.

Siro reads four types of files:
* ARSCA-made netCDF input files
* siro_settings.nml configuration file
* ./input/brdf/<brdf_file>
* ./input/miefiles/<scatterer_file>

The first two files are created and modified using ARSCA. The last two files need to be provided in the corresponding folders. Examples of the file formats are provided.

## Known issues
* Polarization of the BRDF surface is disabled due to unverified.
* Some aerosols may cause abnormally high shot noise. The cause is most likely inadequate handling of singularities in the polaris subroutine.
* Refraction is hard-coded to according to some pressure profile and it is not properly functioning for arbitrary profiles currently.

## References
[1]: https://doi.org/10.1029/1999JD900969
