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

test_bf_layer_basic:		parameters_input.o\
				parameters_kind.o\
				$(bf_layer_basic_dep)

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
				$(bf_layer_print_dep)

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

test_bf_layer_sync.o:		$(pbf_layer_dir)/bf_layer_sync_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_layer_sync:		$(param_dep)\
				bf_layer_errors_module.o\
				parameters_bf_layer.o\
				$(bf_layer_sync_dep)\
				check_data_module.o

test_bf_newgrdpt_prim.o:	$(cbf_layer_dir)/bf_newgrdpt_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o\
				$(mt_dir)/sd_operators_fd_module.o\
				$(dim2d_dir)/dim2d_parameters.o

test_bf_newgrdpt_prim:		check_data_module.o\
				$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				n_coords_module.o\
				bf_newgrdpt_class.o\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				bf_layer_newgrdpt_procedure_module.o

test_bf_newgrdpt_procedure.o:	$(cbf_layer_dir)/bf_newgrdpt_class.o\
				$(cbf_layer_dir)/bf_newgrdpt_procedure_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o\
				$(mt_dir)/sd_operators_fd_module.o\
				$(dim2d_dir)/dim2d_parameters.o



test_bf_newgrdpt_procedure:	check_data_module.o\
				$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				n_coords_module.o\
				bf_newgrdpt_class.o\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				bf_layer_newgrdpt_procedure_module.o\
				bf_newgrdpt_procedure_module.o

test_bf_newgrdpt_verification.o:$(cbf_layer_dir)/bf_newgrdpt_verification_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o


test_bf_newgrdpt_verification:	parameters_bf_layer.o\
				bf_newgrdpt_verification_module.o\
				check_data_module.o

test_bf_layer_bc_procedure.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_procedure:\
			parameters_kind.o\
			parameters_bf_layer.o\
			bf_layer_bc_procedure_module.o