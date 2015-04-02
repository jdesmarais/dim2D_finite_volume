ifeq ($(cluster), newcluster)

#intel libraries: parallel
ifeq ($(AUGEANSTABLES_PARALLEL),true)

	INCLUDE_IFORT  = -I /cm/shared/apps/openmpi/intel/64/1.6.5/include
	INCLUDE_IFORT += -I /cm/shared/apps/netcdf/ifort/64/4.2/include

	LIBS_IFORT     = -lcurl
	LIBS_IFORT    += -L /cm/shared/apps/netcdf/ifort/64/4.2/lib -lnetcdff -lnetcdf
	LIBS_IFORT    += -L /cm/shared/apps/hdf5/1.8.13/intel64/lib -lhdf5_hl -lhdf5 -lz
	LIBS_IFORT    += -L /cm/shared/apps/openmpi/intel/64/1.6.5/lib

#intel libraries: sequential
else

	INCLUDE_IFORT = -I /cm/shared/apps/netcdf/ifort/64/4.2/include
	LIBS_IFORT    = -L /cm/shared/apps/netcdf/ifort/64/4.2/lib -lnetcdff -lnetcdf

endif

#gfortran libraries
	INCLUDE_GFORTRAN = -I /cm/shared/apps/netcdf/gfortran/64/4.2-gcc481/include
	LIBS_GFORTRAN    = -L /cm/shared/apps/netcdf/gfortran/64/4.2-gcc481/lib -lnetcdff -lnetcdf

	INCLUDE_GCC_PAR	 = -I /cm/shared/apps/openmpi/gcc/64/1.6.5/include
	LIBS_GCC_PAR	 = -L /cm/shared/apps/openmpi/gcc/64/1.6.5/lib64
endif