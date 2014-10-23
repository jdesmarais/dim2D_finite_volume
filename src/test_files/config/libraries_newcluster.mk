ifeq ($(cluster), newcluster)
#intel libraries
	INCLUDE_IFORT = -I /cm/shared/apps/netcdf/open64/64/4.3.0/include
	LIBS_IFORT    = -L /cm/shared/apps/netcdf/open64/64/4.3.0/lib

#GFORTRAN libraries
	INCLUDE_GFORTRAN = -I /cm/shared/apps/netcdf/gcc/64/4.3.0/include
	LIBS_GFORTRAN    = -L /cm/shared/apps/netcdf/gcc/64/4.3.0/lib

	INCLUDE_GCC_PAR	 = -I /cm/shared/apps/openmpi/gcc/64/1.6.5/include
	LIBS_GCC_PAR	 = -L /cm/shared/apps/openmpi/gcc/64/1.6.5/lib64
endif