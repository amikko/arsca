
OBJECTS =   parameters.o \
	    global_variables.o \
	    basic_functions.o \
	    atmosphere.o \
	    refra.o \
	    polarisation.o \
            routines.o \
            singlescattering.o \
            netcdf_io.o \
	    siro.o\

COMPILE = gfortran

COPTIONS = -c -O3 -fbounds-check -mtune=native -W -Wall -pedantic -Wconversion -I/usr/local/include -L/usr/local/lib -lnetcdff -L/lib -lnetcdf -lnetcdf -ldl -lm -fopenmp

EXECUTABLE = siro

${EXECUTABLE} : ${OBJECTS}
	${COMPILE} -o ${EXECUTABLE} ${OBJECTS} -I/usr/local/include -L/usr/local/lib -lnetcdff -fopenmp

parameters.o : parameters.f90
	${COMPILE} ${COPTIONS} parameters.f90 
global_variables.o : global_variables.f90
	${COMPILE} ${COPTIONS} global_variables.f90
basic_functions.o : basic_functions.f90
	${COMPILE} ${COPTIONS} basic_functions.f90 
atmosphere.o : atmosphere.f90
	${COMPILE} ${COPTIONS} atmosphere.f90
refra.o : refra.f90
	${COMPILE} ${COPTIONS} refra.f90
polarisation.o : polarisation.f90
	${COMPILE} ${COPTIONS} polarisation.f90
routines.o : routines.f90
	${COMPILE} ${COPTIONS} routines.f90
singlescattering.o : singlescattering.f90
	${COMPILE} ${COPTIONS} singlescattering.f90
netcdf_io.o : netcdf_io.f90
	${COMPILE} ${COPTIONS} netcdf_io.f90
siro.o : siro.f90 
	${COMPILE} ${COPTIONS} siro.f90

clean : 
	\rm -f *.o *.L *.l *.m "#"*"#" *~ "#."*"#" .*~ *.mod siro







