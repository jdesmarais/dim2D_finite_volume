#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

check_data_module.o:		$(param_dir)/parameters_kind.o

parameters_cst.o:

nf90_operators_truncate_module.o:\
				parameters_cst.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

# test files
#-----------------------------------------------------------------------
test_nf90_operators_truncate.o:	$(test_dir)/tools/check_data_module.o\
				nf90_operators_truncate_module.o\
				parameters_cst.o\
				$(param_dir)/parameters_kind.o

test_nf90_operators_truncate:	parameters_cst.o\
				parameters_kind.o\
				nf90_operators_truncate_module.o\
				check_data_module.o


