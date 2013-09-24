#------------------------------------------------------------------------
#cluster
#------------------------------------------------------------------------
cluster = $(LERNAEANHYDRA_CLUSTER)


#------------------------------------------------------------------------
#user main options for profiling, optimizing ...
#------------------------------------------------------------------------
profile       = $(LERNAEANHYDRA_PROFILE)        #ddt/totalview debugging
analyse       = $(LERNAEANHYDRA_ANALYSE)        #gprof
trace         = $(LERNAEANHYDRA_TRACE)          #vampir tracing
report        = $(LERNAEANHYDRA_OPTREPORT)      #optimization report
profileGuided = $(LERNAEANHYDRA_PROFILE_GUIDED) #profile guided
scalasca      = false				#scalasca instrumentation


#------------------------------------------------------------------------
#user main options
#-----------------------------------------------------------------------
src_dir	   = $(lernaeanhydra)/src
config_dir = $(src_dir)/test_files/config
dep_dir	   = $(src_dir)/test_files/dep
bc_choice  = wall_xy_choice

ifeq ($(strip $(bc_choice)), periodic_xy_choice)
	sbc_dir=$(pbc_dir)
	sim_dim2d_dep=$(dep_dir)/sim_dim2d_periodic.dep
	sim_dim2d_par_dep=$(dep_dir)/sim_dim2d_par_periodic.dep
endif
ifeq ($(strip $(bc_choice)), reflection_xy_choice)
	sbc_dir=$(rbc_dir)
	sim_dim2d_dep=$(dep_dir)/sim_dim2d_reflection.dep
	sim_dim2d_par_dep=$(dep_dir)/sim_dim2d_par_reflection.dep
endif
ifeq ($(strip $(bc_choice)), wall_xy_choice)
	sbc_dir=$(wbc_dir)
	sim_dim2d_dep=$(dep_dir)/sim_dim2d_wall.dep
	sim_dim2d_par_dep=$(dep_dir)/sim_dim2d_par_wall.dep
endif
#------------------------------------------------------------------------
#source files directories
#------------------------------------------------------------------------
#define the folder paths for the source code compilation
#------------------------------------------------------------------------
include $(config_dir)/makefile_folders


#------------------------------------------------------------------------
#compiler and options
#------------------------------------------------------------------------
PREP=
FC=$(LERNAEANHYDRA_COMPILER)
FFLAGS	= 
LDFLAGS =
INCLUDE =
LIBS	=

#ifort compiler options------------------------------------------------
#FFLAGS_OP_IFORT		 = -O2 -p -g -ipo-c -132 -inline-factor=300 -warn all
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=ipo_inl -inline-factor=500
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=ipo_align -inline-factor=500 #check if the data are aligned
FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -inline-factor=500 #-warn noalignment #check if there are data that are unaligned
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -vec-report3 -inline-factor=300
#FFLAGS_OP_IFORT		 = -O2 -132 -ipo-c #-D	ALIGNED #-inline-factor=400
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -vec-report3 -guide -inline-factor=400
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=hlo_unroll -inline-factor=400

FFLAGS_PROFILE_IFORT	 = -O0 -g -132 -warn all -check all -r8 -traceback	#profiling options 
FFLAGS_OPT_REPORT_IFORT  = -opt-report3 -opt-report-phase hpo_vectorization -opt-report-file report_intel.txt

#LDFLAGS_OP_IFORT     	 = -O2 -p -g -ipo -132 -inline-factor=300 -warn all
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=ipo_inl -inline-factor=500 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=ipo_align -inline-factor=500 #check if the data are aligned
LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -inline-factor=500 #-warn noalignment #check if there are data that are unaligned
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -vec-report3 -inline-factor=300 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -132 -ipo #-DALIGNED #-inline-factor=400 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -vec-report3 -guide -inline-factor=400 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=hlo_unroll -inline-factor=400 #-warn all -check all #GOOD options
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
include $(config_dir)/bolt_libraries
include $(config_dir)/cartesius_libraries

#options handeling depending on the user main options: debug, trace, analyse--
#scalasca instrumentation
ifeq ($(strip $(scalasca)),true)
	PREP=skin
endif


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
	@(echo 'make clear    : remove .err .out and .nc')
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
	@(echo 'opening documentation...')
	@(firefox ./doc/html/files.html &)

test:
	@(echo '')

	@(echo 'source code folder: ' $(src_dir))
	@(echo 'cluster:            ' $(cluster))
	@(echo '')

	@(echo 'profile:   ' $(profile))
	@(echo 'analyse:   ' $(analyse))
	@(echo 'trace:     ' $(trace))
	@(echo 'report:    ' $(report))
	@(echo 'scalasca:  ' $(scalasca))
	@(echo '')

	@(echo 'PREP:    ' $(PREP))
	@(echo 'FC:      ' $(FC))
	@(echo 'FFLAGS:  ' $(FFLAGS))
	@(echo 'LDFLAGS: ' $(LDFLAGS))
	@(echo 'INCLUDE: ' $(INCLUDE))
	@(echo 'LIBS:    ' $(LIBS))
	@(echo '')

tests:	test_field\
	test_cg_operators\
	test_dim2d_prim\
	test_dim2d_fluxes\
	test_dim2d_eq\
	test_fv_operators\
	test_bc_operators\
	test_rk3tvd\
	test_nf90_operators\
	test_dim2d_ic

tests_par: test_field_par\
	test_mpi_mg_bc

sims:	sim_dim2d

%.o:	%.f
	@(basename $@)
	@($(PREP) $(FC) $(FFLAGS) -c $< $(INCLUDE))

%:	%.o
	@(basename $@)
	@($(PREP) $(FC) $(LDFLAGS) -o $@ $^ $(LIBS))

all:	$(CMD)					#main commmand

clear:
	rm *.err *.out *.nc

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
$(param_dir)/parameters_input.o:$(param_dir)/parameters_constant.o


#fields
$(field_dir)/surrogate_class.o:
$(field_dir)/field_class.o:	$(field_dir)/surrogate_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(field_dir)/field_par_class.o:	$(field_dir)/field_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
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
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/cg_operators_class.o

$(simpletest_dir)/simpletest_eq_class.o:$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(phy_eq_dir)/phy_model_eq_class.o\
				$(sd_dir)/cg_operators_class.o

$(dim2d_dir)/dim2d_parameters.o:$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_prim_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_fluxes_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/cg_operators_class.o

$(dim2d_dir)/dim2d_state_eq_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_eq_class.o:	$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_fluxes_module.o\
				$(dim2d_ic)/dim2d_steadystate_module.o\
				$(dim2d_ic)/dim2d_dropretract_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(phy_eq_dir)/phy_model_eq_class.o

$(dim2d_ic)/dim2d_steadystate_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_state_eq_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_dropretract_module.o:$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_state_eq_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o


#time discretization methods
$(td_dir)/td_operators_class.o:	$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(phy_eq_dir)/phy_model_eq_class.o

$(fv_dir)/fv_operators_class.o:	$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(phy_eq_dir)/phy_model_eq_class.o\
				$(td_dir)/td_operators_class.o

#boundary conditions
$(bc_dir)/bc_abstract_class.o:	$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(bc_dir)/bc_abstract_par_class.o:$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o

#periodic boundary conditions
$(pbc_dir)/bc_operators_class.o:$(bc_dir)/bc_abstract_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(pbc_dir)/periodic_xy_module.o:$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(pbc_dir)/periodic_xy_par_module.o:$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(mpi_dir)/mpi_tag_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(pbc_dir)/bc_operators_par_class.o:$(bc_dir)/bc_abstract_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(pbc_dir)/periodic_xy_par_module.o

#reflection boundary conditions
$(rbc_dir)/reflection_xy_module.o:$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(rbc_dir)/reflection_xy_par_module.o:$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(mpi_dir)/mpi_tag_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rbc_dir)/reflection_xy_module.o

$(rbc_dir)/bc_operators_class.o:$(bc_dir)/bc_abstract_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rbc_dir)/reflection_xy_module.o

$(rbc_dir)/bc_operators_par_class.o:$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(bc_dir)/bc_abstract_par_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(rbc_dir)/reflection_xy_par_module.o

#wall boundary conditions
$(wbc_dir)/wall_prim_module.o:	$(dim2d_dir)/dim2d_parameters.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_kind.o

$(wbc_dir)/wall_xy_module.o:	$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wbc_dir)/wall_prim_module.o

$(wbc_dir)/bc_operators_class.o:$(bc_dir)/bc_abstract_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wbc_dir)/wall_xy_module.o

#time integration methods
$(ti_dir)/td_integrator_class.o:$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_kind.o

$(ti_dir)/td_integrator_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_kind.o

$(rk3tvd_dir)/rk3tvd_steps_class.o:$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(rk3tvd_dir)/rk3tvd_class.o:	$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(ti_dir)/td_integrator_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_class.o

$(rk3tvd_dir)/rk3tvd_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(ti_dir)/td_integrator_par_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_class.o

#io operators
$(io_dir)/io_operators_module.o:

$(io_dir)/io_operators_class.o:	$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(phy_eq_dir)/phy_model_eq_class.o
				$(nf90_dir)/nf90_operators_wr_class.o

$(nf90_dir)/nf90_operators_module.o:$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(dim2d_dir)/dim2d_parameters.o

$(nf90_dir)/nf90_operators_wr_class.o:$(field_dir)/field_class.o\
				$(io_dir)/io_operators_module.o\
				$(nf90_dir)/nf90_operators_module.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/dim2d_eq_class.o


$(nf90_dir)/nf90_operators_wr_par_class.o:$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(io_dir)/io_operators_module.o\
				$(mpi_dir)/mpi_process_class.o\
				$(nf90_dir)/nf90_operators_module.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\

#mpi processes
$(mpi_dir)/mpi_tag_module.o:

$(mpi_dir)/mpi_process_class.o:	$(param_dir)/parameters_input.o


#mpi messenger bc
$(mpi_bc_dir)/mpi_mg_neighbours.o:$(param_dir)/parameters_constant.o

$(mpi_bc_dir)/mpi_mg_derived_types.o:$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o

$(mpi_bc_dir)/mpi_mg_bc_class.o:$(field_dir)/field_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(mpi_bc_dir)/mpi_mg_neighbours.o\
				$(mpi_bc_dir)/mpi_mg_derived_types.o\
				$(param_dir)/parameters_input.o

$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o:$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o

$(mpi_bc_dir)/mpi_mg_construct.o:$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o

$(mpi_bc_dir)/mpi_mg_bc_ext_class.o:$(field_dir)/field_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(mpi_bc_dir)/mpi_mg_construct.o\
				$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o\
				$(param_dir)/parameters_input.o



#simulations
sim_dim2d.o:			$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(nf90_dir)/nf90_operators_wr_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_class.o

sim_dim2d_par.o:		$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(nf90_dir)/nf90_operators_wr_par_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_par_class.o

include $(sim_dim2d_dep)
include $(sim_dim2d_par_dep)


#tests
test_field.o:			$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(field_dir)/field_class.o

test_field:			field_class.o\
				parameters_input.o\
				parameters_kind.o\
				surrogate_class.o

test_cg_operators.o:		$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/cg_operators_class.o

test_cg_operators:		dim2d_prim_module.o\
				field_class.o\
				interface_primary.o\
				parameters_input.o\
				parameters_kind.o\
				cg_operators_class.o\
				surrogate_class.o

test_dim2d_prim.o:		$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_prim:		dim2d_parameters.o\
				dim2d_prim_module.o\
				field_class.o\
				parameters_input.o\
				parameters_kind.o\
				surrogate_class.o

test_dim2d_fluxes.o:		$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_fluxes_module.o\
				$(field_dir)/field_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_fluxes:		cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				field_class.o\
				parameters_kind.o\
				surrogate_class.o

test_dim2d_eq.o:		$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_eq:			cg_operators_class.o\
				dim2d_eq_class.o\
				dim2d_fluxes_module.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_dropretract_module.o\
				dim2d_steadystate_module.o\
				dim2d_state_eq_module.o\
				field_class.o\
				parameters_input.o\
				parameters_kind.o\
				surrogate_class.o

test_fv_operators.o:		$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/simpletest_eq_class.o

test_fv_operators:		cg_operators_class.o\
				field_class.o\
				fv_operators_class.o\
				parameters_input.o\
				parameters_kind.o\
				simpletest_eq_class.o\
				surrogate_class.o\
				td_operators_class.o

test_bc_operators.o:		$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_operators:		bc_operators_class.o\
				cg_operators_class.o\
				field_class.o\
				parameters_constant.o\
				parameters_input.o\
				parameters_kind.o\
				periodic_xy_module.o\
				reflection_xy_module.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				surrogate_class.o

test_bc_periodic.o:		$(pbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_periodic:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o

test_bc_reflection.o:		$(rbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_reflection:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				reflection_xy_module.o\
				bc_operators_class.o

test_wall_xy_module.o:		$(wbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wbc_dir)/wall_xy_module.o

test_wall_xy_module:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				wall_prim_module.o\
				wall_xy_module.o\
				bc_operators_class.o


test_bc_wall.o:			$(wbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_wall:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				wall_prim_module.o\
				wall_xy_module.o\
				bc_operators_class.o

test_rk3tvd.o:			$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/simpletest_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_class.o

test_rk3tvd:			bc_operators_class.o\
				cg_operators_class.o\
				field_class.o\
				fv_operators_class.o\
				parameters_constant.o\
				parameters_input.o\
				parameters_kind.o\
				periodic_xy_module.o\
				phy_model_eq_class.o\
				rk3tvd_class.o\
				simpletest_eq_class.o\
				surrogate_class.o\
				td_integrator_class.o\
				td_operators_class.o

test_rk3tvd_dim2d.o:		$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_class.o

test_rk3tvd_dim2d:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				periodic_xy_module.o\
				bc_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				td_operators_class.o\
				fv_operators_class.o\
				td_integrator_class.o\
				rk3tvd_class.o


test_nf90_operators.o:		$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(nf90_dir)/nf90_operators_wr_class.o

test_nf90_operators:		surrogate_class.o\
				parameters_kind.o\
				parameters_input.o\
				field_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_state_eq_module.o\
				dim2d_dropretract_module.o\
				dim2d_steadystate_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				parameters_constant.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				nf90_operators_module.o\
				nf90_operators_wr_class.o\
				io_operators_module.o

test_dim2d_ic.o:		$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(nf90_dir)/nf90_operators_wr_class.o

test_dim2d_ic:			surrogate_class.o\
				parameters_input.o\
				parameters_kind.o\
				field_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				dim2d_state_eq_module.o\
				parameters_constant.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				nf90_operators_module.o\
				nf90_operators_wr_class.o\
				io_operators_module.o



test_field_par.o:		$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_kind.o

test_field_par:			field_par_class.o\
				field_class.o\
				surrogate_class.o\
				mpi_process_class.o\
				parameters_input.o\
				parameters_kind.o

test_mpi_mg_bc.o:		$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_mpi_mg_bc:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_process_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o

test_mpi_mg_bc_ext.o:		$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o

test_mpi_mg_bc_ext:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_process_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o

test_reflection_xy_par_module.o:$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rbc_dir)/reflection_xy_par_module.o	

test_reflection_xy_par_module:	parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				reflection_xy_module.o\
				reflection_xy_par_module.o

test_reflection_xy_par.o:	$(rbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_reflection_xy_par:		parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				reflection_xy_module.o\
				reflection_xy_par_module.o\
				bc_abstract_par_class.o\
				bc_operators_par_class.o

test_mpi_mg_construct.o:	$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_bc_dir)/mpi_mg_construct.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o


test_mpi_mg_construct:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_mg_construct.o

test_periodic_xy_par_module.o:	$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pbc_dir)/periodic_xy_par_module.o

test_periodic_xy_par_module:	parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o

test_periodic_xy_par.o:		$(pbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_periodic_xy_par:		parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o\
				bc_abstract_par_class.o\
				bc_operators_par_class.o

test_rk3tvd_dim2d_par.o:	$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_par_class.o

test_rk3tvd_dim2d_par:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				mpi_process_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_construct.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_bc_ext_class.o\
				bc_abstract_par_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o\
				bc_operators_par_class.o\
				td_operators_class.o\
				fv_operators_class.o\
				td_integrator_par_class.o\
				rk3tvd_par_class.o

test_nf90_operators_par.o:	$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(nf90_dir)/nf90_operators_wr_par_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_nf90_operators_par:	surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropretract_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				mpi_process_class.o\
				field_par_class.o\
				io_operators_module.o\
				nf90_operators_module.o\
				nf90_operators_wr_par_class.o




