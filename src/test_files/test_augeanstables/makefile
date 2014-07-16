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
			bf_compute_class.o