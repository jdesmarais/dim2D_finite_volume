#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_gaussian_perturbation.o:	$(test_dir)/tools/check_data_module.o\
				$(field_dir)/field_class.o\
				$(dim2d_ic)/gaussian_perturbation/gaussian_perturbation_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_gaussian_perturbation:	check_data_module.o\
				gaussian_perturbation_module.o\
				parameters_constant.o\
				parameters_input.o\
				parameters_kind.o\
				$(field_dep)\
				$(sim_dep)
