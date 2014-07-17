#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_bf_compute_prog.o:	$(bf_layer_dir)/parameters_bf_layer.o\
			$(cbf_layer_dir)/bf_compute_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_compute_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			bf_compute_class.o


test_bf_layer_prog.o:	$(bf_layer_dir)/parameters_bf_layer.o\
			$(bf_layer_dir)/bf_layer_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_layer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)


test_bf_mainlayer_prog.o:$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_mainlayer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)\
			bf_sublayer_class.o\
			bf_mainlayer_class.o
