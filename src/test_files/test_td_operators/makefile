#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_td_operators.o:	$(pbc_dir)/bc_operators_class.o\
			$(test_dir)/tools/check_data_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(simpletest_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(td_cdir)/td_operators_class.o

test_td_operators:	$(param_dep)\
			$(mt_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			check_data_module.o
