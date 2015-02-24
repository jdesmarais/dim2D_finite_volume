#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

test_bf_layer_extract.o:	$(sbf_layer_dir)/bf_layer_extract_module.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_layer_extract:		bf_layer_extract_module.o\
				parameters_input.o\
				parameters_kind.o


test_bf_layer_basic.o:		$(pbf_layer_dir)/bf_layer_basic_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_layer_basic:		bf_layer_extract_module.o\
				bf_layer_basic_class.o\
				parameters_input.o\
				parameters_kind.o

test_bf_layer_nf90_operators_prog.o:\
				$(iobf_layer_dir)/bf_layer_nf90_operators_module.o\
				$(dim2d_dir)/pmodel_eq_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(param_dir)/parameters_input.o

test_bf_layer_nf90_operators_prog:\
				$(param_dep)\
				bf_layer_errors_module.o\
				parameters_bf_layer.o\
				bf_layer_nf90_operators_module.o\
				$(sd_dep)\
				sd_operators_fd_module.o\
				$(dim2d_dep)

test_bf_layer_print.o:		$(pbf_layer_dir)/bf_layer_print_class.o\
				$(dim2d_dir)/pmodel_eq_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(param_dir)/parameters_input.o

test_bf_layer_print:		$(param_dep)\
				bf_layer_errors_module.o\
				parameters_bf_layer.o\
				$(sd_dep)\
				sd_operators_fd_module.o\
				$(dim2d_dep)\
				bf_layer_extract_module.o\
				bf_layer_basic_class.o\
				bf_layer_nf90_operators_module.o\
				bf_layer_print_class.o

test_bf_layer_exchange.o:	$(sbf_layer_dir)/bf_layer_exchange_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_layer_exchange:		bf_layer_errors_module.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				parameters_bf_layer.o\
				bf_layer_exchange_module.o\
				check_data_module.o
