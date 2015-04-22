#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

parameters_cst.o:		

nf90_operators_merge_module.o:	parameters_cst.o\
				$(param_dir)/parameters_kind.o

cmd_operators_merge_class.o:


merge_netcdf_file.o:		cmd_operators_merge_class.o\
				nf90_operators_merge_module.o

merge_netcdf_file:		cmd_operators_merge_class.o\
				nf90_operators_merge_module.o\
				parameters_cst.o	


merge_netcdf_files.o:		cmd_operators_merge_class.o\
				nf90_operators_merge_module.o

merge_netcdf_files:		cmd_operators_merge_class.o\
				nf90_operators_merge_module.o\
				parameters_cst.o


merge_netcdf_files_par.o:	cmd_operators_merge_class.o\
				nf90_operators_merge_module.o

merge_netcdf_files_par:		cmd_operators_merge_class.o\
				nf90_operators_merge_module.o\
				parameters_cst.o

# test files
#-----------------------------------------------------------------------
test_nf90_operators_merge.o:	nf90_operators_merge_module.o\
				parameters_cst.o\
				$(param_dir)/parameters_kind.o

test_nf90_operators_merge:	parameters_cst.o\
				parameters_kind.o\
				nf90_operators_merge_module.o


