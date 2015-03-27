#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

#bf_layer_sync
test_bf_layer_extract.o:	$(sbf_layer_dir)/bf_layer_extract_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
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


test_bf_layer_newgrdpt.o:	$(pbf_layer_dir)/bf_layer_newgrdpt_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_bf_layer_newgrdpt:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_layer_newgrdpt_dep)\
				check_data_module.o


test_bf_layer_grdpts_id_update.o:$(pbf_layer_dir)/bf_layer_grdpts_id_update_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_grdpts_id_update:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_layer_grdpts_id_update_dep)\
				check_data_module.o


test_bf_layer_icr.o:		$(pbf_layer_dir)/bf_layer_icr_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o

test_bf_layer_icr:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_layer_icr_dep)\
				check_data_module.o


test_bf_layer.o:		$(bf_layer_dir)/bf_layer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o

test_bf_layer:			$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_layer_dep)


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

test_bf_newgrdpt_dispatch.o:	$(ngbf_layer_dir)/bf_newgrdpt_dispatch_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_bf_newgrdpt_dispatch:	bf_layer_errors_module.o\
				check_data_module.o\
				parameters_bf_layer.o\
				parameters_kind.o\
				parameters_constant.o\
				bf_layer_extract_module.o\
				bf_newgrdpt_class.o\
				bf_newgrdpt_procedure_module.o\
				bf_newgrdpt_verification_module.o\
				bf_newgrdpt_dispatch_module.o\
				$(mt_ext_dep)\
				$(dim2d_dep)


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
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_interior_bc_sections:	$(param_dep)\
				parameters_bf_layer.o\
				bf_interior_bc_sections_module.o\
				check_data_module.o


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

test_bf_mainlayer_bc_sections.o:$(bbf_layer_dir)/bf_mainlayer_bc_sections_module.o\
				$(bf_layer_dir)/bf_mainlayer_pointer_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_mainlayer_bc_sections:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_mainlayer_pointer_dep)\
				bf_mainlayer_bc_sections_module.o\
				check_data_module.o

#buffer_layer/bf_layer_computations
test_bf_compute_basic.o:	$(cbf_layer_dir)/bf_compute_basic_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_compute_basic:		$(bf_compute_basic_dep)\
				check_data_module.o\
				parameters_bf_layer.o\
				parameters_input.o\
				parameters_kind.o

test_bf_compute_time.o:		$(hobc_dir)/bc_operators_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(cbf_layer_dir)/bf_compute_time_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o\
				$(mt_dir)/sd_operators_class.o\
				$(fv_dir)/td_operators_class.o

test_bf_compute_time:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				rk3tvd_steps_module.o\
				$(bf_compute_time_dep)\
				check_data_module.o

test_bf_compute_newgrdpt.o:	$(cbf_layer_dir)/bf_compute_newgrdpt_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_bf_compute_newgrdpt:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_compute_newgrdpt_dep)\
				check_data_module.o




#main objects
test_bf_sublayer.o:		$(bf_layer_dir)/bf_sublayer_class.o

test_bf_sublayer:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(bf_layer_dep)\
				$(fv_dep)\
				bf_sublayer_class.o


#bf_mainlayer parents
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

test_bf_mainlayer_time.o:	$(hobc_dir)/bc_operators_class.o\
				$(mpbf_layer_dir)/bf_mainlayer_time_class.o\
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

test_bf_mainlayer_time:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(bf_layer_dep)\
				$(fv_dep)\
				$(bf_mainlayer_time_dep)\
				check_data_module.o\
				rk3tvd_steps_module.o

test_bf_mainlayer.o:		$(bf_layer_dir)/bf_mainlayer_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_mainlayer:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(bf_layer_dep)\
				$(fv_dep)\
				$(bf_mainlayer_dep)


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

test_mainlayer_interface_dyn.o:	$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(mbf_layer_dir)/mainlayer_interface_dyn_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_mainlayer_interface_dyn:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_dyn_dep)\
				check_data_module.o

test_mainlayer_interface_newgrdpt.o:\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(mbf_layer_dir)/mainlayer_interface_newgrdpt_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o


test_mainlayer_interface_newgrdpt:\
				$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_newgrdpt_dep)\
				check_data_module.o


test_mainlayer_interface_grdpts_id_update.o:\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(mbf_layer_dir)/mainlayer_interface_grdpts_id_update_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_mainlayer_interface_grdpts_id_update:\
				$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(mainlayer_interface_grdpts_id_update_dep)\
				check_data_module.o

#bf_interface parents
test_bf_interface_basic.o:	$(ipbf_layer_dir)/bf_interface_basic_class.o\
				$(bf_layer_dir)/bf_mainlayer_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_interface_basic:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_basic_dep)

test_bf_interface_print.o:	$(ipbf_layer_dir)/bf_interface_print_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_interface_print:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_print_dep)

test_bf_interface_sync.o:	$(ipbf_layer_dir)/bf_interface_sync_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

test_bf_interface_sync:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_sync_dep)\
				check_data_module.o

test_bf_interface_dyn.o:	$(ipbf_layer_dir)/bf_interface_dyn_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(test_dir)/tools/check_data_module.o

test_bf_interface_dyn:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_dyn_dep)\
				check_data_module.o

test_bf_interface_time.o:	$(bc_cdir)/bc_operators_class.o\
				$(ipbf_layer_dir)/bf_interface_time_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(sd_cdir)/sd_operators_class.o\
				$(td_cdir)/td_operators_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o

test_bf_interface_time:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_time_dep)\
				check_data_module.o\
				rk3tvd_steps_module.o

test_bf_interface_coords.o:	$(ipbf_layer_dir)/bf_interface_coords_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bf_interface_coords:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_coords_dep)

test_bf_interface_icr.o:	$(ipbf_layer_dir)/bf_interface_icr_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(ibf_layer_dir)/icr_interface_class.o\
				$(ibf_layer_dir)/icr_path_chain_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o

test_bf_interface_icr:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_icr_dep)\
				check_data_module.o


test_bf_interface_dcr.o:	$(ipbf_layer_dir)/bf_interface_dcr_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wave2d_dir)/pmodel_eq_class.o

test_bf_interface_dcr:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_dcr_dep)\
				check_data_module.o

#test bc_crenels
test_bf_bc_interior_pt_crenel.o:$(gbf_layer_dir)/bf_bc_interior_pt_crenel_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o

test_bf_bc_interior_pt_crenel:	$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				bf_newgrdpt_class.o\
				bf_newgrdpt_procedure_module.o\
				bf_newgrdpt_verification_module.o\
				bf_newgrdpt_dispatch_module.o\
				bf_bc_interior_pt_crenel_module.o\
				parameters_bf_layer.o


test_bf_bc_pt_crenel.o:		$(gbf_layer_dir)/bf_bc_pt_crenel_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o

test_bf_bc_pt_crenel:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				bf_newgrdpt_class.o\
				bf_newgrdpt_procedure_module.o\
				bf_newgrdpt_verification_module.o\
				bf_newgrdpt_dispatch_module.o\
				bf_bc_pt_crenel_module.o\
				parameters_bf_layer.o\
				check_data_module.o

#test bf_layer_increase
test_bf_increase_coords.o:	$(ibf_layer_dir)/bf_increase_coords_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o

test_bf_increase_coords:	$(param_dep)\
				parameters_bf_layer.o\
				bf_increase_coords_module.o


test_bf_sorting.o:		$(ibf_layer_dir)/bf_sorting_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(param_dir)/parameters_kind.o

test_bf_sorting:		parameters_kind.o\
				check_data_module.o\
				bf_sorting_module.o


test_icr_path.o:		$(ipbf_layer_dir)/bf_interface_coords_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(ibf_layer_dir)/icr_path_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_icr_path:			$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_coords_dep)\
				$(icr_path_dep)\
				check_data_module.o


test_icr_path_chain.o:		$(ibf_layer_dir)/icr_path_chain_class.o

test_icr_path_chain:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_coords_dep)\
				$(icr_path_chain_dep)


test_icr_path_list.o:		$(ibf_layer_dir)/icr_path_list_class.o

test_icr_path_list:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_coords_dep)\
				$(icr_path_list_dep)


test_icr_interface.o:		$(ipbf_layer_dir)/bf_interface_coords_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(test_dir)/tools/check_data_module.o\
				$(ibf_layer_dir)/icr_interface_class.o\
				$(ibf_layer_dir)/icr_path_chain_class.o\
				$(ibf_layer_dir)/icr_path_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o

test_icr_interface:		$(param_dep)\
				$(mt_dep)\
				$(dim2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_coords_dep)\
				$(icr_interface_dep)\
				check_data_module.o


#test bf_layer_decrease
test_bf_decrease.o:		$(dbf_layer_dir)/bf_decrease_module.o\
				$(sbf_layer_dir)/bf_layer_extract_module.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wave2d_dir)/pmodel_eq_class.o

test_bf_decrease:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				parameters_bf_layer.o\
				bf_layer_errors_module.o\
				bf_layer_extract_module.o\
				check_data_module.o\
				bf_decrease_module.o


test_dcr_chain.o:		$(dbf_layer_dir)/dcr_chain_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o

test_dcr_chain:			$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_sublayer_dep)\
				$(dcr_chain_dep)


test_dcr_list.o:		$(dbf_layer_dir)/dcr_list_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o

test_dcr_list:			$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_sublayer_dep)\
				$(dcr_list_dep)


test_dcr_interface.o:		$(dbf_layer_dir)/dcr_interface_class.o\
				$(ipbf_layer_dir)/bf_interface_icr_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wave2d_dir)/pmodel_eq_class.o

test_dcr_interface:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(bf_interface_icr_dep)\
				$(dcr_interface_dep)
