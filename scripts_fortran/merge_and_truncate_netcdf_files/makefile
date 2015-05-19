#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

parameters_cst.o:		

nf90_operators_mt_module.o:	file_merge_param_class.o\
				parameters_cst.o\
				$(param_dir)/parameters_kind.o

file_merge_param_class.o:	$(param_dir)/parameters_kind.o

cmd_operators_mt_class.o:	$(param_dir)/parameters_kind.o

progress_mg_module.o:

merge_and_truncate_netcdf_files.o:\
				cmd_operators_mt_class.o\
				file_merge_param_class.o\
				nf90_operators_mt_module.o\
				$(param_dir)/parameters_kind.o\
				progress_mg_module.o

merge_and_truncate_netcdf_files:\
				cmd_operators_mt_class.o\
				file_merge_param_class.o\
				nf90_operators_mt_module.o\
				parameters_kind.o\
				progress_mg_module.o

# test files
#-----------------------------------------------------------------------
test_nf90_operators_mt.o:	$(test_dir)/tools/check_data_module.o\
				nf90_operators_mt_module.o\
				$(param_dir)/parameters_kind.o

test_nf90_operators_mt:		check_data_module.o\
				nf90_operators_mt_module.o\
				parameters_kind.o

test_file_merge_param.o:	$(test_dir)/tools/check_data_module.o\
				file_merge_param_class.o\
				$(param_dir)/parameters_kind.o

test_file_merge_param:		check_data_module.o\
				file_merge_param_class.o\
				parameters_kind.o



