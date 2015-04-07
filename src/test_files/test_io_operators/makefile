#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_io_operators.o:	$(io_dir)/io_operators_module.o\
			$(test_dir)/tools/check_data_module.o

test_io_operators:	parameters_kind.o\
			io_operators_module.o\
			check_data_module.o
