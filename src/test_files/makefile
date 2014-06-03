#-----------------------------------------------------------------------
#cluster
#-----------------------------------------------------------------------
cluster = $(ERYMANTHIANBOAR_CLUSTER)


#-----------------------------------------------------------------------
#user main options for profiling, optimizing ...
#-----------------------------------------------------------------------
profile       = $(ERYMANTHIANBOAR_PROFILE)      #ddt/totalview debugging
analyse       = $(ERYMANTHIANBOAR_ANALYSE)      #gprof
trace         = $(ERYMANTHIANBOAR_TRACE)        #vampir tracing
report        = $(ERYMANTHIANBOAR_OPTREPORT)    #optimization report
profileGuided = $(ERYMANTHIANBOAR_PROFILE_GUIDED)#profile guided
scalasca      = false			       #scalasca instrumentation


#-----------------------------------------------------------------------
#user main options
#-----------------------------------------------------------------------
src_dir	   = $(erymanthianboar)/src
config_dir = $(src_dir)/test_files/config
dep_dir	   = $(src_dir)/test_files/dep
pm_choice  = simpletest_choice
bc_choice  = reflection_xy_choice

ifeq ($(strip $(pm_choice)), simpletest_choice)
	spm_dir=$(simpletest_dir)
endif
ifeq ($(strip $(pm_choice)), dim2d_choice)
	spm_dir=$(dim2d_dir)
endif

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
ifeq ($(strip $(bc_choice)), wall_x_reflection_y_choice)
	sbc_dir=$(wrbc_dir)
	sim_dim2d_dep=$(dep_dir)/sim_dim2d_wall_reflection.dep
	sim_dim2d_par_dep=$(dep_dir)/sim_dim2d_par_wall_reflection.dep
endif


#-----------------------------------------------------------------------
#compiler and options
#-----------------------------------------------------------------------
PREP	=
FC	= $(ERYMANTHIANBOAR_COMPILER)
FFLAGS	= 
LDFLAGS =
INCLUDE =
LIBS	=

#compiler options------------------------------------------------
include $(config_dir)/options_ifort.mk
include $(config_dir)/options_gfortran.mk

#path for the mpi, hdf5, netcdf libraries depending on the cluster used-
include $(config_dir)/libraries_bolt.mk
include $(config_dir)/libraries_cartesius.mk

#options handeling depending on the user main options: debug, trace,
#and the type of compiler used (ifort, gfortran, mpiifort...)
include $(config_dir)/options_choice.mk


#-----------------------------------------------------------------------
#source code files of the current dir
#-----------------------------------------------------------------------
SRC	= $(wildcard *.f)	#fortran source files
SRCS	= $(SRC)		#fortran program files
OBJS	= $(SRC:.f=.o)		#objects (.o and .mod)
CMD	= $(SRCS:.f=)		#executables


#-----------------------------------------------------------------------
#source files directories
#-----------------------------------------------------------------------
#define the folder paths for the source code compilation
#-----------------------------------------------------------------------
include $(config_dir)/makefile_folders.mk


#----------------------------------------------------------------------
#general rules
#----------------------------------------------------------------------
help:
	@(echo 'make all      : create executable for all test files')
	@(echo 'make clean    : remove .o files')
	@(echo 'make cleanall : remove .o and executable files')
	@(echo 'make clear    : remove .err .out and .nc')
	@(echo 'make test     : display makefile options')

	@(echo 'to change the cluster configuration, change the env. var. ERYMANTHIANBOAR_CLUSTER')
	@(echo 'to change the compiler, change the env. var. ERYMANTHIANBOAR_COMPILER')
	@(echo 'to activate profiling, change the env. var. ERYMANTHIANBOAR_PROFILE')
	@(echo 'to activate analyzing, change the env. var. ERYMANTHIANBOAR_ANALYSE')
	@(echo 'to activate tracing, change the env. var. ERYMANTHIANBOAR_TRACE')
	@(echo 'to activate report, change the env. var. ERYMANTHIANBOAR_OPTREPORT')
	@(echo 'all these variables can also be set in the corresponding module file')
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

	@(echo $(test_dir))

tests:	test_field\
	test_cg_operators\
	test_dim2d_prim\
	test_dim2d_fluxes\
	test_dim2d_eq\
	test_fv_operators\
	test_bc_periodic\
	test_bc_reflection\
	test_wall_xy_module\
	test_bc_wall\
	test_rk3tvd\
	test_nf90_operators\
	test_dim2d_ic\
	test_rk3tvd_dim2d

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


#-----------------------------------------------------------------------
#dependencies
#-----------------------------------------------------------------------

#buffer layer dep
$(bf_layer_dir)/parameters_bf_layer.o:\
	$(param_dir)/parameters_input.o

$(bf_layer_dir)/bf_layer_errors_module.o:

$(bf_layer_dir)/bf_layer_ini_grdptID_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_allocate_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_reallocate_module.o:\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_merge_module.o:\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_exchange_module.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_class.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/bf_layer_reallocate_module.o\
	$(bf_layer_dir)/bf_layer_merge_module.o\
	$(bf_layer_dir)/bf_layer_exchange_module.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_sublayer_class.o:\
	$(bf_layer_dir)/bf_layer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(nbf_layer_dir)/nbf_element_class.o:\
	$(bf_layer_dir)/bf_layer_class.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(nbf_layer_dir)/nbf_list_class.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(param_dir)/parameters_kind.o

$(nbf_layer_dir)/nbf_interface_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(bf_layer_dir)/parameters_bf_layer.o

$(bf_layer_dir)/bf_mainlayer_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

#$(bf_layer_dir)/bf_layer_path_abstract_class.o:\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(param_dir)/parameters_kind.o

#$(bf_layer_dir)/bf_layer_corner_check_module.o:\
#	$(bf_layer_dir)/bf_layer_path_abstract_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/bf_mainlayer_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o

#$(bf_layer_dir)/bf_sublayers_merge_module.o:\
#	$(bf_layer_dir)/bf_mainlayer_abstract_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/parameters_bf_layer.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o

#$(bf_layer_dir)/bf_layer_update_allocation_module.o:\
#	$(bf_layer_dir)/bf_layer_path_abstract_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o
#


$(bf_layer_dir)/bf_mainlayer_pointer_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

#$(bf_layer_dir)/interface_abstract_class.o:\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/bf_mainlayer_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_pointer_class.o\
	$(nbf_layer_dir)/nbf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_path_class.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/dbf_element_class.o:\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/dbf_list_class.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_i_list_class.o:\
	$(bf_layer_dir)/bf_interface_class.o\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(dbf_layer_dir)/dbf_list_class.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_nbc_template_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o

$(bf_layer_dir)/bf_interface_icr_class.o:\
	$(dbf_layer_dir)/bf_detector_i_list_class.o\
	$(bf_layer_dir)/bf_layer_path_class.o\
	$(bf_layer_dir)/bf_nbc_template_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

#tests for buffer layer
test_bf_layer_module.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	$(bf_layer_dir)/bf_layer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_prog.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	$(bf_layer_dir)/bf_layer_class.o\
	test_bf_layer_module.o

test_bf_sublayer_prog.o:\
	$(bf_layer_dir)/bf_sublayer_class.o

test_bf_mainlayer_module.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o

test_bf_mainlayer_prog.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_bf_mainlayer_module.o

test_nbf_element_prog.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_nbf_list_module.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

test_nbf_list_prog.o:\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_nbf_interface_prog.o:\
	$(nbf_layer_dir)/nbf_interface_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_bf_layer_path_prog.o:\
	$(bf_layer_dir)/bf_layer_path_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

#test_interface_abstract_prog.o:\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/bf_mainlayer_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o\
#	test_bf_mainlayer_module.o

test_bf_interface_prog.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_dbf_element_prog.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(param_dir)/parameters_kind.o

test_dbf_list_prog.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(dbf_layer_dir)/dbf_list_class.o\
	$(param_dir)/parameters_kind.o

test_bf_detector_i_list_prog.o:\
	$(dbf_layer_dir)/bf_detector_i_list_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_icr_prog.o:\
	$(dbf_layer_dir)/bf_detector_i_list_class.o\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

#test_bf_corner_prog.o:\
#	$(bf_layer_dir)/bf_corner_module.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o

test_cases_interface_module.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

#test_bf_layer_adapt_corner_prog.o:\
#	$(bf_layer_dir)/bf_layer_adapt_corner_module.o\
#	$(bf_layer_dir)/bf_interface_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o

test_bf_layer_update_grdpts_id_prog.o:\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_cases_interface_module.o

#test_cases_path_module.o:\
#	$(bf_layer_dir)/bf_layer_path_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o

#interface_print_module.o:\
#	$(bf_layer_dir)/bf_mainlayer_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o

#sublayers_ini_module.o:\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o

#test_bf_layer_corner_check_prog.o:\
#	$(bf_layer_dir)/bf_layer_corner_check_module.o\
#	$(bf_layer_dir)/bf_layer_path_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	test_cases_path_module.o\
#	interface_print_module.o

#test_bf_sublayers_merge_prog.o:\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/bf_sublayers_merge_module.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	test_cases_path_module.o\
#	interface_print_module.o\
#	sublayers_ini_module.o
#
#test_bf_layer_update_allocation_prog.o:\
#	$(bf_layer_dir)/bf_mainlayer_class.o\
#	$(bf_layer_dir)/bf_layer_update_allocation_module.o\
#	$(bf_layer_dir)/bf_layer_path_class.o\
#	$(bf_layer_dir)/bf_sublayer_class.o\
#	$(bf_layer_dir)/interface_abstract_class.o\
#	$(param_dir)/parameters_constant.o\
#	$(param_dir)/parameters_input.o\
#	$(param_dir)/parameters_kind.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	test_cases_path_module.o\
#	sublayers_ini_module.o\
#	interface_print_module.o

test_bf_layer_reallocate_prog.o:\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_cases_interface_module.o

test_bf_layer_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_bf_layer.o\
	parameters_input.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	test_bf_layer_module.o


test_bf_sublayer_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o


test_bf_mainlayer_prog:\
	parameters_bf_layer.o\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_layer_ini_grdptID_module.o\
	bf_layer_merge_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_bf_mainlayer_module.o

test_nbf_element_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	nbf_element_class.o\
	test_bf_layer_module.o

test_nbf_list_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_nbf_interface_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o




#test_interface_abstract_prog:\
#	parameters_kind.o\
#	bf_layer_abstract_class.o\
#	parameters_bf_layer.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_layer_exchange_module.o\
#	bf_layer_compute_module.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_layer_allocate_module.o\
#	bf_layer_print_module.o\
#	bf_layer_class.o\
#	bf_sublayer_class.o\
#	bf_mainlayer_abstract_class.o\
#	bf_sublayers_merge_module.o\
#	bf_mainlayer_class.o\
#	interface_abstract_class.o\
#	bf_layer_update_grdpts_module.o\
#	test_bf_layer_module.o\
#	test_bf_mainlayer_module.o

test_bf_interface_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	test_bf_layer_module.o

test_bf_layer_path_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	bf_layer_path_class.o\
	test_bf_layer_module.o

test_dbf_element_prog:\
	dbf_element_class.o\
	parameters_kind.o

test_dbf_list_prog:\
	dbf_element_class.o\
	dbf_list_class.o\
	parameters_kind.o

test_bf_detector_i_list_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_i_list_class.o\
	test_bf_layer_module.o

test_bf_interface_icr_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_i_list_class.o\
	bf_layer_path_class.o\
	bf_nbc_template_module.o\
	bf_interface_icr_class.o\
	test_bf_layer_module.o


#test_bf_corner_prog:\
#	parameters_kind.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_corner_module.o\
#	parameters_bf_layer.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_layer_merge_module.o\
#	bf_layer_class.o\
#	test_bf_layer_module.o

#test_bf_layer_adapt_corner_prog:\
#	parameters_bf_layer.o\
#	parameters_kind.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_layer_merge_module.o\
#	bf_layer_adapt_corner_module.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_corner_module.o\
#	bf_layer_class.o\
#	bf_sublayer_class.o\
#	bf_mainlayer_class.o\
#	bf_mainlayer_pointer_class.o\
#	bf_interface_class.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o


#test_bf_layer_corner_check_prog:\
#	parameters_kind.o\
#	bf_layer_abstract_class.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_layer_exchange_module.o\
#	bf_layer_compute_module.o\
#	parameters_bf_layer.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_layer_allocate_module.o\
#	bf_layer_print_module.o\
#	bf_layer_class.o\
#	bf_sublayer_class.o\
#	bf_layer_path_abstract_class.o\
#	bf_mainlayer_abstract_class.o\
#	bf_sublayers_merge_module.o\
#	bf_mainlayer_class.o\
#	interface_abstract_class.o\
#	bf_layer_corner_check_module.o\
#	bf_layer_path_class.o\
#	bf_layer_update_grdpts_module.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	test_cases_path_module.o\
#	interface_print_module.o
#
#test_bf_sublayers_merge_prog:\
#	parameters_kind.o\
#	bf_layer_abstract_class.o\
#	parameters_bf_layer.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_layer_exchange_module.o\
#	bf_layer_compute_module.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_layer_allocate_module.o\
#	bf_layer_print_module.o\
#	bf_layer_class.o\
#	bf_sublayer_class.o\
#	bf_sublayers_merge_module.o\
#	bf_mainlayer_class.o\
#	interface_abstract_class.o\
#	bf_layer_update_grdpts_module.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	bf_layer_path_abstract_class.o\
#	bf_layer_corner_check_module.o\
#	bf_layer_path_class.o\
#	test_cases_path_module.o\
#	interface_print_module.o\
#	sublayers_ini_module.o
#
#test_bf_layer_update_allocation_prog:\
#	parameters_kind.o\
#	bf_layer_abstract_class.o\
#	parameters_bf_layer.o\
#	parameters_constant.o\
#	parameters_input.o\
#	bf_layer_exchange_module.o\
#	bf_layer_compute_module.o\
#	bf_layer_ini_grdptID_module.o\
#	bf_layer_allocate_module.o\
#	bf_layer_print_module.o\
#	bf_layer_class.o\
#	bf_sublayer_class.o\
#	bf_layer_path_abstract_class.o\
#	bf_mainlayer_abstract_class.o\
#	bf_sublayers_merge_module.o\
#	bf_mainlayer_class.o\
#	interface_abstract_class.o\
#	bf_layer_update_allocation_module.o\
#	bf_layer_corner_check_module.o\
#	bf_layer_path_class.o\
#	bf_layer_update_grdpts_module.o\
#	test_bf_layer_module.o\
#	test_cases_interface_module.o\
#	test_cases_path_module.o\
#	sublayers_ini_module.o\
#	interface_print_module.o

test_bf_layer_reallocate_prog:\
	parameters_bf_layer.o\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	bf_interface_class.o\
	test_bf_layer_module.o\
	test_cases_interface_module.o	

test_bf_layer_update_grdpts_id_prog:\
	parameters_bf_layer.o\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_layer_ini_grdptID_module.o\
	bf_layer_merge_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	bf_interface_class.o\
	test_bf_layer_module.o\
	test_cases_interface_module.o


#-----------------------------------------------------------------------
#parameter_dir						
$(param_dir)/parameters_kind.o:

$(param_dir)/parameters_constant.o:\
			$(param_dir)/parameters_kind.o

$(param_dir)/parameters_input.o:\
			$(param_dir)/parameters_constant.o

#fields
$(field_dir)/surrogate_class.o:
$(field_dir)/field_class.o:\
			$(field_dir)/surrogate_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(field_dir)/field_par_class.o:\
			$(field_dir)/field_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o


#space discretization operators
$(sd_dir)/interface_primary.o:\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_kind.o

$(sd_dir)/sd_operators_class.o:\
			$(field_dir)/field_class.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o

$(sd_dir)/cg_operators_class.o:\
			$(field_dir)/field_class.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/sd_operators_class.o

#physical models
$(phy_eq_dir)/phy_model_eq_class.o:\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/cg_operators_class.o

#simple test equations
$(simpletest_dir)/dim2d_eq_class.o:\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/phy_model_eq_class.o\
			$(sd_dir)/cg_operators_class.o

#diffuse interface model equations
$(dim2d_dir)/dim2d_parameters.o:\
			$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_prim_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_fluxes_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/cg_operators_class.o

$(dim2d_dir)/dim2d_state_eq_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(param_dir)/parameters_kind.o

$(dim2d_dir)/dim2d_eq_class.o:\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_ic)/dim2d_bubble_ascending_module.o\
			$(dim2d_ic)/dim2d_drop_collision_module.o\
			$(dim2d_ic)/dim2d_drop_retraction_module.o\
			$(dim2d_dir)/dim2d_fluxes_module.o\
			$(dim2d_ic)/dim2d_homogeneous_module.o\
			$(dim2d_ic)/dim2d_phase_separation_module.o\
			$(dim2d_ic)/dim2d_steadystate_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/phy_model_eq_class.o

$(dim2d_ic)/dim2d_steadystate_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_dropbubble_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_vortex_module.o:\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_bubble_ascending_module.o:\
			$(dim2d_ic)/dim2d_dropbubble_module.o\
			$(dim2d_ic)/dim2d_vortex_module.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_drop_collision_module.o:\
			$(dim2d_ic)/dim2d_dropbubble_module.o\
			$(dim2d_ic)/dim2d_vortex_module.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_drop_retraction_module.o:\
			$(dim2d_ic)/dim2d_dropbubble_module.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_phase_separation.o:\
			$(dim2d_ic)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(dim2d_ic)/dim2d_homogeneous_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

#boundary conditions
$(bc_dir)/bc_abstract_class.o:\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(bc_dir)/bc_abstract_par_class.o:\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o

#periodic boundary conditions
$(pbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_abstract_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(pbc_dir)/periodic_xy_module.o:\
			$(sd_dir)/cg_operators_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(pbc_dir)/periodic_xy_par_module.o:\
			$(field_dir)/field_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(pbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_abstract_par_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pbc_dir)/periodic_xy_par_module.o

#reflection boundary conditions
$(rbc_dir)/reflection_xy_module.o:\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(rbc_dir)/reflection_xy_par_module.o:\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o

$(rbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_abstract_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o

$(rbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_abstract_par_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(spm_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_par_module.o

#wall boundary conditions
$(wbc_dir)/wall_prim_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_kind.o

$(wbc_dir)/wall_xy_module.o:\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_prim_module.o

$(wbc_dir)/wall_xy_par_module.o:\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

$(wbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_abstract_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

$(wbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_abstract_par_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_par_module.o

#wall_x reflection_y boundary conditions
$(wrbc_dir)/wall_x_reflection_y_par_module.o:\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(wbc_dir)/wall_xy_module.o

$(wrbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_abstract_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(wbc_dir)/wall_xy_module.o

$(wrbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_abstract_par_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(dim2d_dir)/dim2d_eq_class.o\
			$(field_dir)/field_par_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wrbc_dir)/wall_x_reflection_y_par_module.o\
			$(wbc_dir)/wall_xy_par_module.o

#time discretization methods
$(td_dir)/td_operators_class.o:\
			$(sbc_dir)/bc_operators_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(spm_dir)/dim2d_eq_class.o

$(td_dir)/td_operators_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(spm_dir)/dim2d_eq_class.o

$(fv_dir)/fv_operators_class.o:	$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(td_dir)/td_operators_class.o

$(fv_dir)/fv_operators_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(td_dir)/td_operators_par_class.o


#time integration methods
$(ti_dir)/td_integrator_class.o:$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_kind.o

$(ti_dir)/td_integrator_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_par_class.o\
				$(param_dir)/parameters_kind.o

$(rk3tvd_dir)/rk3tvd_steps_class.o:$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(rk3tvd_dir)/rk3tvd_class.o:	$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(ti_dir)/td_integrator_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_class.o

$(rk3tvd_dir)/rk3tvd_par_class.o:$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(spm_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_par_class.o\
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

$(mpi_dir)/mpi_requests_module.o:$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(mpi_dir)/mpi_tag_module.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

$(mpi_dir)/mpi_process_class.o:	$(param_dir)/parameters_input.o


#mpi messenger bc
$(mpi_bc_dir)/mpi_mg_neighbours.o:$(param_dir)/parameters_constant.o

$(mpi_bc_dir)/mpi_mg_derived_types.o:$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
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
				$(fv_dir)/fv_operators_par_class.o\
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

test_dim2d_eq:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_state_eq_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_steadystate_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o

test_fv_operators.o:		$(pbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/dim2d_eq_class.o

test_fv_operators:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o\
				td_operators_class.o\
				fv_operators_class.o

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
				dim2d_dropbubble_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_vortex_module.o\
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
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_state_eq_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_steadystate_module.o\
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
				dim2d_dropbubble_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_vortex_module.o\
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
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				wall_prim_module.o\
				wall_xy_module.o\
				bc_operators_class.o

test_rk3tvd.o:			$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/dim2d_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_class.o

test_rk3tvd:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o\
				td_operators_class.o\
				fv_operators_class.o\
				td_integrator_class.o\
				rk3tvd_steps_class.o\
				rk3tvd_class.o

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
				dim2d_dropletbubble_module.o\
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
				dim2d_dropbubble_module.o\
				dim2d_steadystate_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
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
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
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
				dim2d_dropletbubble_module.o\
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
				dim2d_dropletbubble_module.o\
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
				dim2d_dropletbubble_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
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
				$(fv_dir)/fv_operators_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
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
				dim2d_dropletbubble_module.o\
				dim2d_homogeneous_liquid.o\
				dim2d_phase_separation_module.o\
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
				td_operators_par_class.o\
				fv_operators_par_class.o\
				td_integrator_par_class.o\
				rk3tvd_steps_class.o\
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
				dim2d_dropletbubble_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				mpi_process_class.o\
				field_par_class.o\
				io_operators_module.o\
				nf90_operators_module.o\
				nf90_operators_wr_par_class.o




