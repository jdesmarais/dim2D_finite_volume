#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

#test_mpi_process
test_mpi_process.o:	$(mpi_dir)/mpi_process_class.o

test_mpi_process:	mpi_process_class.o


#test_mpi_tag
test_mpi_tag.o:		$(mpi_dir)/mpi_tag_module.o
test_mpi_tag:		mpi_tag_module.o


#test_mpi_derived_types
test_mpi_derived_types.o:$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_derived_types.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_derived_types:	$(param_dep)\
			check_data_module.o\
			mpi_mg_derived_types.o\
			mpi_process_class.o

#test_mpi_neighbours
test_mpi_mg_neighbours.o:$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_neighbours.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_input.o

test_mpi_mg_neighbours:	$(param_dep)\
			check_data_module.o\
			mpi_process_class.o\
			mpi_mg_neighbours.o


#test_mpi_mg_bc
test_mpi_mg_bc.o:	$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_mg_bc:		$(param_dep)\
			check_data_module.o\
			$(mpi_mg_bc_dep)\
			mpi_process_class.o


#test_mpi_requests
test_mpi_requests.o:	$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_requests:	$(param_dep)\
			check_data_module.o\
			$(mpi_mg_bc_dep)\
			mpi_process_class.o\
			mpi_tag_module.o\
			mpi_requests_module.o


#test_mpi_mg_ini_bc_proc
test_mpi_mg_ini_bc_proc.o:$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o

test_mpi_mg_ini_bc_proc:$(param_dep)\
			mpi_mg_ini_bc_proc.o\
			mpi_process_class.o


#test_mpi_mg_construct
test_mpi_mg_construct.o:$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_bc_dir)/mpi_mg_construct.o\
			$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_mg_construct:	$(param_dep)\
			$(mpi_mg_bc_dep)\
			mpi_mg_construct.o\
			mpi_mg_ini_bc_proc.o\
			mpi_process_class.o\
			mpi_requests_module.o\
			mpi_tag_module.o\
			check_data_module.o


#test_mpi_mg_construct
test_mpi_mg_bc_ext.o:	$(test_dir)/tools/check_data_module.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_mg_bc_ext:	$(param_dep)\
			$(mpi_mg_bc_ext_dep)\
			mpi_process_class.o\
			mpi_requests_module.o\
			mpi_tag_module.o\
			check_data_module.o

