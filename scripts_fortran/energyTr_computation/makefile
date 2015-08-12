#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

error_dir = ../error_computation
test_dir  = $(augeanstables)/src/test_files
tools_dir = $(test_dir)/tools


# module files
#-----------------------------------------------------------------------
check_data_module.o:				$(param_dir)/parameters_kind.o

cmd_operators_extract_class.o:			$(param_dir)/parameters_kind.o


$(error_dir)/parameters_cst.o:

$(error_dir)/compare_type_module.o:		$(param_dir)/parameters_kind.o

$(error_dir)/nf90_operators_error_module.o:	$(error_dir)/compare_type_module.o\
						$(error_dir)/parameters_cst.o\
						$(param_dir)/parameters_kind.o

energyTr_computation_module.o:			$(param_dir)/parameters_kind.o\
						$(error_dir)/parameters_cst.o


# main files
#-----------------------------------------------------------------------
extract_energyTr.o:		cmd_operators_extract_class.o\
				energyTr_computation_module.o\
				$(error_dir)/nf90_operators_error_module.o\
				$(error_dir)/parameters_cst.o\
				parameters_kind.o

extract_energyTr:		cmd_operators_extract_class.o\
				energyTr_computation_module.o\
				nf90_operators_error_module.o\
				parameters_kind.o\
				compare_type_module.o\
				parameters_cst.o


# test files
#-----------------------------------------------------------------------
test_energyTr_computation.o:	$(tools_dir)/check_data_module.o\
				energyTr_computation_module.o\
				$(param_dir)/parameters_kind.o

test_energyTr_computation:	parameters_kind.o\
				check_data_module.o\
				parameters_cst.o\
				energyTr_computation_module.o


