#------------------------------------------------------------------------
#cluster
#------------------------------------------------------------------------
cluster = $(LERNAEANHYDRA_CLUSTER)

#------------------------------------------------------------------------
#user main options
#------------------------------------------------------------------------
profile       = $(LERNAEANHYDRA_PROFILE)        #ddt/totalview debugging
analyse       = $(LERNAEANHYDRA_ANALYSE)        #gprof
trace         = $(LERNAEANHYDRA_TRACE)          #vampir tracing
report        = $(LERNAEANHYDRA_OPTREPORT)      #optimization report
profileGuided = $(LERNAEANHYDRA_PROFILE_GUIDED) #profile guided

#------------------------------------------------------------------------
#source files directories
#------------------------------------------------------------------------
src_dir	= /home/jdesmara/Code/lernaeanhydra_opt/src
include $(src_dir)/test_files/config/makefile_folders

#------------------------------------------------------------------------
#compiler and options
#------------------------------------------------------------------------
FC=$(LERNAEANHYDRA_COMPILER)
FFLAGS	= 
LDFLAGS =
INCLUDE =
LIBS	=

#ifort compiler options------------------------------------------------
FFLAGS_OP_IFORT		 = -O2 -132 -ipo-c -opt-report3 -opt-report-phase=ipo_inl
FFLAGS_PROFILE_IFORT	 = -O0 -g -132 -warn all -check all -r8 -traceback	#profiling options 
FFLAGS_OPT_REPORT_IFORT  = -opt-report3 -opt-report-phase hpo_vectorization -opt-report-file report_intel.txt

LDFLAGS_OP_IFORT     	 = -O2 -132 -ipo -opt-report3 -opt-report-phase=ipo_inl #-warn all -check all #GOOD options
LDFLAGS_PROFILE_IFORT	 = -O0 -g -i-dynamic -132 -r8			#profiling options
LDFLAGS_OPT_REPORT_IFORT = -opt-report3 -opt-report-phase hpo_vectorization -opt-report-file report_intel.txt

#gfortran compiler options---------------------------------------------
FFLAGS_OP_GFORTRAN     	= -O3 -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #optimized options
FFLAGS_PROFILE_GFORTRAN	= -g  -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #profiling options
FFLAGS_ANALYSE_GFORTRAN	= -pg -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #analyse options


LDFLAGS_OP_GFORTRAN	= -O3 -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #optimized options
LDFLAGS_PROFILE_GFORTRAN= -g  -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #profiling options
LDFLAGS_ANALYSE_GFORTRAN= -pg -ffixed-line-length-132 -Wall -Wextra -pedantic -fbacktrace #analyse options

#path for the mpi, hdf5, netcdf libraries depending on the cluster used----------------------
ifeq ($(cluster),bolt)
#intel librairies
	INCLUDE_IFORT = -I /cm/shared/apps/openmpi/intel/64/current/include
	INCLUDE_IFORT+= -I /cm/shared/apps/netcdf/gcc/64/4.1.1/include
	LIBS_IFORT    = -L /cm/shared/apps/netcdf/gcc/64/4.1.1/lib -netcdf

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


ifeq ($(cluster), cartesius)
#intel libraries
	INCLUDE_IFORT 		= -I $(SURFSARA_IMPI_INCLUDE)
	INCLUDE_IFORT          += -I $(SURFSARA_NETCDF_INCLUDE)

	LIBS_IFORT		 = -lcurl
	LIBS_IFORT		+= -L $(SURFSARA_NETCDF_LIB) -lnetcdff -lnetcdf
	LIBS_IFORT		+= -L $(SURFSARA_HDF5_LIB) -lhdf5_hl -lhdf5 -lz
	LIBS_IFORT		+= -L $(SURFSARA_IMPI_LIB)

#LIBS_IFORT		= -L $(SURFSARA_NETCDF_LIB) -lnetcdff -lnetcdf
#LIBS_IFORT            += -L $(SURFSARA_IMPI_LIB)

ifeq ($(trace),true)
	LIBS_IFORT             += -L $(VT_LIB_DIR) -lVT $(VT_ADD_LIBS)
endif

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


#options handeling depending on the user main options: debug, trace, analyse--
#IFORT options
ifeq ($(FC),ifort)
#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_IFORT)
	LDFLAGS+=$(LDFLAGS_PROFILE_IFORT)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_IFORT)
	LDFLAGS+=$(LDFLAGS_ANALYSE_IFORT)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_IFORT)
	LDFLAGS+=$(LDFLAGS_OP_IFORT)

#OPTIMIZATION REPORT
ifeq ($(strip $(report)),true)
	FFLAGS+=-opt-report3
	FFLAGS+=-opt-report-phase hpo_vectorization
	FFLAGS+=-opt-report-file intel_report.txt 
endif


endif
endif
INCLUDE +=$(INCLUDE_IFORT)
LIBS	+=$(LIBS_IFORT)
endif

#GFORTRAN options
ifeq ($(FC),gfortran)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_PROFILE_GFORTRAN)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_ANALYSE_GFORTRAN)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_OP_GFORTRAN)

endif
endif

INCLUDE +=$(INCLUDE_GFORTRAN)
LIBS	+=$(LIBS_GFORTRAN)
endif

#mpiifort options
ifeq ($(FC),mpiifort)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_IFORT)
	LDFLAGS+=$(LDFLAGS_PROFILE_IFORT)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_IFORT)
	LDFLAGS+=$(LDFLAGS_ANALYSE_IFORT)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_IFORT)
	LDFLAGS+=$(LDFLAGS_OP_IFORT)

#OPTIMIZATION REPORT
ifeq ($(strip $(report)),true)
	FFLAGS+=$(FFLAGS_OPT_REPORT_IFORT)
	LDFLAGS+=$(LDFLAGS_OPT_REPORT_IFORT)
endif

endif
endif

INCLUDE +=$(INCLUDE_IFORT)
LIBS	+=$(LIBS_IFORT)
endif


#mpif90 options
ifeq ($(FC),mpif90)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_PROFILE_GFORTRAN)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_ANALYSE_GFORTRAN)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_OP_GFORTRAN)

endif
endif

INCLUDE +=$(INCLUDE_GFORTRAN)
INCLUDE +=$(INCLUDE_GCC_PAR)
LIBS	+=$(LIBS_GFORTRAN)
LIBS	+=$(LIBS_GCC_PAR)
endif

#------------------------------------------------------------------------
#source code files
#------------------------------------------------------------------------
SRC	= $(wildcard *.f)	#fortran source files
SRCS	= $(SRC)		#fortran program files
OBJS	= $(SRC:.f=.o)		#objects (.o and .mod)
CMD	= $(SRCS:.f=)		#executables

#------------------------------------------------------------------------
#general rules
#------------------------------------------------------------------------
help:
	@(echo 'make all      : create executable for all test files')
	@(echo 'make clean    : remove .o files')
	@(echo 'make cleanall : remove .o and executable files')
	@(echo 'make test     : display makefile options')

	@(echo 'to change the cluster configuration, change the env. var. LERNAEANHYDRA_CLUSTER')
	@(echo 'to change the compiler, change the env. var. LERNAEANHYDRA_COMPILER')
	@(echo 'to activate profiling, change the env. var. LERNAEANHYDRA_PROFILE')
	@(echo 'to activate analyzing, change the env. var. LERNAEANHYDRA_ANALYSE')
	@(echo 'to activate tracing, change the env. var. LERNAEANHYDRA_TRACE')
	@(echo 'to activate report, change the env. var. LERNAEANHYDRA_OPTREPORT')
	@(echo '')

doc:
	./config/generate_doc.sh
	@(echo 'documentation generated')

test:
	@(echo '')

	@(echo 'source code folder: ' $(src_dir))
	@(echo 'cluster:            ' $(cluster))
	@(echo '')

	@(echo 'profile: ' $(profile))
	@(echo 'analyse: ' $(analyse))
	@(echo 'trace:   ' $(trace))
	@(echo 'report:  ' $(report))
	@(echo '')

	@(echo 'compiler: ' $(FC))
	@(echo 'FFLAGS:   ' $(FFLAGS))
	@(echo 'LDFLAGS:  ' $(LDFLAGS))
	@(echo 'INCLUDE:  ' $(INCLUDE))
	@(echo 'LIBS:     ' $(LIBS))
	@(echo '')

%.o:	%.f
	@(basename $@)
	@($(FC) $(FFLAGS) -c $< $(INCLUDE))

%:	%.o
	@(basename $@)
	@($(FC) $(LDFLAGS) -o $@ $^ $(LIBS))

all:	$(CMD)					#main commmand

clean:	
	rm -f *.o *.mod core *~			#cleaning

cleanall:
	rm -f $(subst .exe,_mono.exe,$(CMD))\
	$(CMD) *.o *.mod core *~       		#complete cleaning


#------------------------------------------------------------------------
#dependencies
#------------------------------------------------------------------------

#parameter_dir						
$(param_dir)/parameters_kind.o:
$(param_dir)/parameters_constant.o:$(param_dir)/parameters_kind.o

#fields
$(field_dir)/surrogate_class.o:
$(field_dir)/field_class.o:	$(field_dir)/surrogate_class.o\
				$(param_dir)/parameters_kind.o

#space discretization operators
$(sd_dir)/interface_primary.o:	$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o

$(sd_dir)/sd_operators_class.o:	$(field_dir)/field_class.o\
				$(sd_dir)/interface_primary.o\
				$(param_dir)/parameters_kind.o

$(sd_dir)/cg_operators_class.o:	$(field_dir)/field_class.o\
				$(sd_dir)/interface_primary.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/sd_operators_class.o

#physical models
$(phy_eq_dir)/phy_model_eq_class.o:$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/sd_operators_class.o

#dim2d_equations
$(dim2d_dir)/dim2d_parameters.o:$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_prim_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_eq_class.o:	$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/sd_operators_class.o

#tests
test_field.o:			$(param_dir)/parameters_kind.o\
				$(field_dir)/field_class.o

test_field:			parameters_kind.o\
				field_class.o\
				surrogate_class.o

test_cg_operators.o:		$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/cg_operators_class.o

test_cg_operators:		dim2d_prim_module.o\
				field_class.o\
				interface_primary.o\
				parameters_kind.o\
				cg_operators_class.o\
				surrogate_class.o

test_dim2d_prim.o:		$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o

test_dim2d_prim:		dim2d_parameters.o\
				dim2d_prim_module.o\
				field_class.o\
				parameters_kind.o\
				surrogate_class.o