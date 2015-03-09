#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

#bf_layer_sync
test_bf_layer_extract.o:	$(sbf_layer_dir)/bf_layer_extract_module.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

test_bf_layer_extract:		bf_layer_extract_module.o\
				check_data_module.o\
				parameters_bf_layer.o\
				parameters_input.o\
				parameters_kind.o

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

#bf_layer_parents
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


test_bf_layer_time.o:		$(hobc_dir)/bc_operators_class.o\
				$(pbf_layer_dir)/bf_layer_time_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o\
				$(mt_dir)/sd_operators_class.o\
				$(fv_dir)/td_operators_class.o

test_bf_layer_time:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				rk3tvd_steps_module.o\
				$(bf_layer_time_dep)\
				check_data_module.o


#buffer_layer/bf_layer_newgrdpt
test_bf_newgrdpt_prim.o:	$(ngbf_layer_dir)/bf_newgrdpt_class.o\
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
				bf_layer_errors_module.o

test_bf_newgrdpt_procedure.o:	$(ngbf_layer_dir)/bf_newgrdpt_class.o\
				$(ngbf_layer_dir)/bf_newgrdpt_procedure_module.o\
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
				bf_newgrdpt_procedure_module.o

test_bf_newgrdpt_verification.o:$(ngbf_layer_dir)/bf_newgrdpt_verification_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o

test_bf_newgrdpt_verification:	parameters_bf_layer.o\
				bf_newgrdpt_verification_module.o\
				check_data_module.o


#buffer_layer/bf_layer_bc_sections
test_bf_layer_bc_procedure.o:	$(bbf_layer_dir)/bf_layer_bc_procedure_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_procedure:	parameters_kind.o\
				parameters_bf_layer.o\
				bf_layer_bc_procedure_module.o


test_bf_layer_bc_sections_overlap.o:\
				$(bbf_layer_dir)/bf_layer_bc_sections_overlap_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

test_bf_layer_bc_sections_overlap:\
				bf_layer_bc_sections_overlap_module.o\
				bf_layer_errors_module.o\
				parameters_bf_layer.o\
				check_data_module.o


test_bf_interior_bc_sections.o:	$(bbf_layer_dir)/bf_interior_bc_sections_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_interior_bc_sections:	$(param_dep)\
				parameters_bf_layer.o\
				bf_interior_bc_sections_module.o


test_bf_layer_bc_sections.o:	$(bbf_layer_dir)/bf_layer_bc_sections_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_sections:	$(param_dep)\
				bf_layer_bc_procedure_module.o\
				bf_layer_errors_module.o\
				bf_layer_bc_sections_overlap_module.o\
				bf_layer_bc_sections_class.o\
				check_data_module.o


test_bf_layer_bc_checks.o:	$(bbf_layer_dir)/bf_layer_bc_checks_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_layer_bc_checks:	parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				bf_layer_bc_checks_module.o


test_bf_layer_bc_anticorner.o:	$(bbf_layer_dir)/bf_layer_bc_anticorner_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_bf_layer_bc_anticorner:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				parameters_bf_layer.o\
				n_coords_module.o\
				check_data_module.o\
				bf_layer_extract_module.o\
				bf_layer_bc_anticorner_module.o


#buffer_layer/bf_layer_computations
test_bf_compute.o:		$(hobc_dir)/bc_operators_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(cbf_layer_dir)/bf_compute_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o\
				$(mt_dir)/sd_operators_class.o\
				$(fv_dir)/td_operators_class.o

test_bf_compute:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				rk3tvd_steps_module.o\
				bf_compute_class.o\
				check_data_module.o


test_bf_remove.o:		$(cbf_layer_dir)/bf_remove_module.o\
				$(sbf_layer_dir)/bf_layer_extract_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wave2d_dir)/pmodel_eq_class.o

test_bf_remove:			$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				parameters_bf_layer.o\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				check_data_module.o\
				bf_remove_module.o

#main objects
test_bf_sublayer.o:		$(bf_layer_dir)/bf_sublayer_class.o

test_bf_sublayer:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(bf_layer_dep)\
				$(fv_dep)\
				bf_sublayer_class.o


test_bf_mainlayer.o:		$(hobc_dir)/bc_operators_class.o\
				$(bf_layer_dir)/bf_mainlayer_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o\
				$(mt_dir)/sd_operators_class.o\
				$(fv_dir)/td_operators_class.o

test_bf_mainlayer:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(bf_layer_dep)\
				$(fv_dep)\
				$(bf_mainlayer_dep)\
				check_data_module.o\
				rk3tvd_steps_module.o


#buffer mainlayer parents
test_bf_mainlayer_basic.o:	$(mpbf_layer_dir)/bf_mainlayer_basic_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

test_bf_mainlayer_basic:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_sublayer_dep)\
				bf_mainlayer_basic_class.o\
				check_data_module.o

test_bf_mainlayer_print.o:	$(mpbf_layer_dir)/bf_mainlayer_print_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_mainlayer_print:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_sublayer_dep)\
				$(bf_mainlayer_print_dep)

test_bf_mainlayer_sync.o:	$(mpbf_layer_dir)/bf_mainlayer_sync_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_mainlayer_sync:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_sublayer_dep)\
				$(bf_mainlayer_sync_dep)\
				check_data_module.o


#mainlayer_interface parents
test_mainlayer_interface_basic.o:\
				$(mbf_layer_dir)/mainlayer_interface_basic_class.o\
				$(bf_layer_dir)/bf_layer_errors_module.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o

test_mainlayer_interface_basic: $(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_basic_dep)\
				check_data_module.o

test_mainlayer_interface_sync.o:$(bf_layer_dir)/bf_layer_errors_module.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(mbf_layer_dir)/mainlayer_interface_sync_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o

test_mainlayer_interface_sync:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_sync_dep)\
				check_data_module.o

test_mainlayer_interface_dyn.o:	$(test_dir)/tools/check_data_module.o\
				$(mbf_layer_dir)/mainlayer_interface_dyn_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o

test_mainlayer_interface_dyn:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_dyn_dep)\
				check_data_module.o