#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

parameters_cst.o:		

compare_type_module.o:			$(param_dir)/parameters_kind.o

compute_error_module.o:			compare_type_module.o\
					parameters_cst.o\
					$(param_dir)/parameters_kind.o

nf90_operators_error_module.o:		compare_type_module.o\
					parameters_cst.o\
					$(param_dir)/parameters_kind.o

nf90_operators_error_max_module.o:	nf90_operators_error_module.o\
					parameters_cst.o\
					$(param_dir)/parameters_kind.o

error_netcdf_file_module.o:		$(param_dir)/parameters_kind.o

cmd_operators_error_class.o:

cmd_operators_error_max_class.o:	

compute_error_file.o:			cmd_operators_error_class.o\
					compare_type_module.o\
					compute_error_module.o\
					nf90_operators_error_module.o\
					parameters_cst.o\
					$(param_dir)/parameters_kind.o

compute_error_file:			cmd_operators_error_class.o\
					compare_type_module.o\
					compute_error_module.o\
					nf90_operators_error_module.o\
					parameters_cst.o\
					parameters_kind.o


compute_error_max_file.o:		cmd_operators_error_max_class.o\
					nf90_operators_error_module.o\
					nf90_operators_error_max_module.o\
					$(param_dir)/parameters_kind.o

compute_error_max_file:			parameters_kind.o\
					compare_type_module.o\
					cmd_operators_error_max_class.o\
					nf90_operators_error_module.o\
					nf90_operators_error_max_module.o

test_nf90_operators_error_max.o:	nf90_operators_error_module.o\
					nf90_operators_error_max_module.o\
					$(param_dir)/parameters_kind.o

test_nf90_operators_error_max:		parameters_kind.o\
					compare_type_module.o\
					nf90_operators_error_module.o\
					nf90_operators_error_max_module.o

test_cmd_operators_error_max.o:		cmd_operators_error_max_class.o

test_cmd_operators_error_max:		cmd_operators_error_max_class.o