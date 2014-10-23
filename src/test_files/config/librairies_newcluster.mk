ifeq ($(cluster), newcluster)
#intel libraries
	INCLUDE_IFORT           = -I $(NETCDFINCLUDE)

	LIBS_IFORT		= -L /cm/shared/apps/netcdf/open64/64/4.3.0/lib


#GFORTRAN libraries
	INCLUDE_GFORTRAN = -I /cm/shared/apps/netcdf/gcc/64/4.2.1.1-par/include
	INCLUDE_GFORTRAN+= -I /cm/shared/apps/hdf5/1.8.10-par/include
	INCLUDE_GFORTRAN+= -I /cm/shared/apps/openmpi/gcc/64/1.6.1/include

	LIBS_GFORTRAN = -lcurl
	LIBS_GFORTRAN+= -L /cm/shared/apps/netcdf/gcc/64/4.2.1.1-par/lib -lnetcdff -lnetcdf
	LIBS_GFORTRAN+= -L /cm/shared/apps/hdf5/1.8.10-par/lib -lhdf5_hl -lhdf5 -lmfhdf -ldf
	LIBS_GFORTRAN+= -L /cm/shared/apps/openmpi/gcc/64/1.6.1/lib64

	INCLUDE_GCC_PAR		= -I /cm/shared/apps/openmpi/gcc/64/1.6.1/include
	LIBS_GCC_PAR		= -L /cm/shared/apps/openmpi/gcc/64/1.6.1/lib64
endif