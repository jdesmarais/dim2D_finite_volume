#buffer layer dep
$(bf_layer_dir)/parameters_bf_layer.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_errors_module.o:

$(bf_layer_dir)/bf_layer_allocate_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_reallocate_module.o:\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_merge_module.o:\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_exchange_module.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(iobf_layer_dir)/bf_layer_nf90_operators_module.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(cbf_layer_dir)/bf_remove_module.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o

$(cbf_layer_dir)/bf_suspicious_bc_interior_pt_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(cbf_layer_dir)/bf_layer_bc_procedure_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o

$(cbf_layer_dir)/bf_layer_bc_sections_class.o:\
	$(cbf_layer_dir)/bf_layer_bc_procedure_module.o

$(cbf_layer_dir)/bf_interior_bc_sections_module.o:\
	$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(cbf_layer_dir)/bf_layer_newgrdpt_procedure_module.o:\
	$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(cbf_layer_dir)/bf_newgrdpt_class.o:\
	$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
	$(cbf_layer_dir)/bf_layer_newgrdpt_procedure_module.o\
	$(sd_dir)/interface_primary.o\
	$(sd_dir)/n_coords_module.o\
	$(obc_dir)/openbc_operators_module.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(sd_cdir)/sd_operators_fd_module.o\
	$(sd_dir)/sd_operators_fd_n_module.o

$(cbf_layer_dir)/bf_compute_class.o:\
	$(bc_cdir)/bc_operators_class.o\
	$(cbf_layer_dir)/bf_layer_newgrdpt_procedure_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(cbf_layer_dir)/bf_newgrdpt_class.o\
	$(ti_dir)/interface_integration_step.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(sd_cdir)/sd_operators_class.o\
	$(td_cdir)/td_operators_class.o


$(bf_layer_dir)/bf_layer_class.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(bc_cdir)/bc_operators_class.o\
	$(cbf_layer_dir)/bf_compute_class.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/bf_layer_reallocate_module.o\
	$(bf_layer_dir)/bf_layer_merge_module.o\
	$(iobf_layer_dir)/bf_layer_nf90_operators_module.o\
	$(bf_layer_dir)/bf_layer_exchange_module.o\
	$(cbf_layer_dir)/bf_remove_module.o\
	$(cbf_layer_dir)/bf_suspicious_bc_interior_pt_module.o\
	$(ti_dir)/interface_integration_step.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(sd_cdir)/sd_operators_class.o\
	$(td_cdir)/td_operators_class.o

$(bf_layer_dir)/bf_sublayer_class.o:\
	$(bf_layer_dir)/bf_layer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(sbf_layer_dir)/bf_sublayer_pointer_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o

$(sbf_layer_dir)/sbf_list_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(sbf_layer_dir)/bf_sublayer_pointer_class.o

$(nbf_layer_dir)/nbf_element_class.o:\
	$(bf_layer_dir)/bf_layer_class.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(nbf_layer_dir)/nbf_list_class.o:\
	$(cbf_layer_dir)/bf_interior_bc_sections_module.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(sbf_layer_dir)/sbf_list_class.o

$(nbf_layer_dir)/nbf_interface_class.o:\
	$(cbf_layer_dir)/bf_interior_bc_sections_module.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(sbf_layer_dir)/sbf_list_class.o

$(nbf_layer_dir)/nbf_interface_newgrdpt_class.o:\
	$(cbf_layer_dir)/bf_layer_newgrdpt_procedure_module.o\
	$(cbf_layer_dir)/bf_suspicious_bc_interior_pt_module.o\
	$(cbf_layer_dir)/bf_newgrdpt_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o	

$(bf_layer_dir)/bf_mainlayer_class.o:\
	$(cbf_layer_dir)/bf_interior_bc_sections_module.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(ti_dir)/interface_integration_step.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_mainlayer_pointer_class.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(ti_dir)/interface_integration_step.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_restart_module.o:\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_class.o:\
	$(cbf_layer_dir)/bf_interior_bc_sections_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_pointer_class.o\
	$(iobf_layer_dir)/bf_layer_nf90_operators_module.o\
	$(bf_layer_dir)/bf_restart_module.o\
	$(ti_dir)/interface_integration_step.o\
	$(nbf_layer_dir)/nbf_interface_newgrdpt_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(sbf_layer_dir)/sbf_list_class.o

$(bf_layer_dir)/bf_path_icr_class.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o

$(dbf_layer_dir)/dbf_element_class.o:\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/dbf_list_class.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_module.o:\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_icr_list_class.o:\
	$(dbf_layer_dir)/bf_detector_module.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_param_class.o:\
	$(dbf_layer_dir)/bf_detector_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_class.o:\
	$(dbf_layer_dir)/bf_detector_module.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_NS_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_EW_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_N_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_NS_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_S_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_NS_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_E_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_EW_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_dcr_list_W_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_list_EW_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_nbc_template_module.o:\
	$(bf_layer_dir)/parameters_bf_layer.o

$(bf_layer_dir)/bf_interface_icr_class.o:\
	$(dbf_layer_dir)/bf_detector_dcr_param_class.o\
	$(dbf_layer_dir)/bf_detector_icr_list_class.o\
	$(bf_layer_dir)/bf_path_icr_class.o\
	$(bf_layer_dir)/bf_nbc_template_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_dcr_class.o:\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(pm_cdir)/pmodel_eq_class.o\
	$(sbf_layer_dir)/sbf_list_class.o
