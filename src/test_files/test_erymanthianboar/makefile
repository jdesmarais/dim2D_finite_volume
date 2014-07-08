#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


#-----------------------------------------------------------------------
#makefile folders
#-----------------------------------------------------------------------
bf_test_case_dir = $(test_dir)/test_erymanthianboar/test_cases


#-----------------------------------------------------------------------
#dependencies for the test files objects
#-----------------------------------------------------------------------
test_bf_layer_module.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	$(bf_layer_dir)/bf_layer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_prog.o:\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	$(bf_layer_dir)/bf_layer_class.o\
	test_bf_layer_module.o

test_bf_sublayer_prog.o:\
	$(bf_layer_dir)/bf_sublayer_class.o

test_bf_mainlayer_module.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o

test_bf_mainlayer_prog.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_bf_mainlayer_module.o

test_nbf_element_prog.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_nbf_list_module.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

test_nbf_list_prog.o:\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_nbf_interface_prog.o:\
	$(nbf_layer_dir)/nbf_interface_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_bf_path_icr_prog.o:\
	$(bf_layer_dir)/bf_path_icr_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_prog.o:\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	$(sbf_layer_dir)/sbf_list_class.o

test_dbf_element_prog.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(param_dir)/parameters_kind.o

test_dbf_list_prog.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(dbf_layer_dir)/dbf_list_class.o\
	$(param_dir)/parameters_kind.o

test_bf_detector_icr_list_prog.o:\
	$(dbf_layer_dir)/bf_detector_icr_list_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_icr_prog.o:\
	$(dbf_layer_dir)/bf_detector_icr_list_class.o\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_cases_interface_update_module.o:\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_icr_update_prog.o:\
	$(dbf_layer_dir)/bf_detector_icr_list_class.o\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_layer_remove_prog.o:\
	$(bf_layer_dir)/bf_layer_remove_module.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_dcr_update_prog.o:\
	$(bf_layer_dir)/bf_interface_dcr_class.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

$(bf_test_case_dir)/bubble_class.o:\
	$(param_dir)/parameters_kind.o

$(bf_test_case_dir)/rising_bubble_class.o:\
	$(param_dir)/parameters_kind.o

$(bf_test_case_dir)/growing_bubble_class.o:\
	$(param_dir)/parameters_kind.o

$(bf_test_case_dir)/test_case_class.o:\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_test_case_dir)/bubble_class.o\
	$(bf_test_case_dir)/rising_bubble_class.o\
	$(bf_test_case_dir)/growing_bubble_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o

test_bf_interface_dcr_update_prog2.o:\
	$(bf_layer_dir)/bf_interface_dcr_class.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	test_bf_layer_module.o\
	$(bf_test_case_dir)/test_case_class.o


#-----------------------------------------------------------------------
#dependencies for the test files executables
#-----------------------------------------------------------------------
test_bf_layer_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	test_bf_layer_module.o

test_bf_sublayer_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o

test_bf_mainlayer_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_bf_mainlayer_module.o

test_nbf_element_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	nbf_element_class.o\
	test_bf_layer_module.o

test_nbf_list_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_nbf_interface_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_bf_interface_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	test_bf_layer_module.o

test_bf_path_icr_prog:\
	bf_layer_errors_module.o\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	bf_path_icr_class.o\
	test_bf_layer_module.o

test_dbf_element_prog:\
	dbf_element_class.o\
	parameters_kind.o

test_dbf_list_prog:\
	dbf_element_class.o\
	dbf_list_class.o\
	parameters_kind.o

test_bf_detector_icr_list_prog:\
	parameters_kind.o\
	bf_detector_module.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_icr_list_class.o\
	test_bf_layer_module.o

test_bf_interface_icr_prog:\
	parameters_kind.o\
	bf_detector_module.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_icr_list_class.o\
	bf_detector_dcr_list_class.o\
	bf_detector_dcr_list_NS_class.o\
	bf_detector_dcr_list_N_class.o\
	bf_detector_dcr_list_S_class.o\
	bf_detector_dcr_list_EW_class.o\
	bf_detector_dcr_list_E_class.o\
	bf_detector_dcr_list_W_class.o\
	bf_path_icr_class.o\
	bf_nbc_template_module.o\
	bf_interface_icr_class.o\
	test_bf_layer_module.o

test_bf_interface_icr_update_prog:\
	parameters_kind.o\
	bf_detector_module.o\
	parameters_constant.o\
	parameters_input.o\
	parameters_bf_layer.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_activation_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_icr_list_class.o\
	bf_detector_dcr_list_class.o\
	bf_detector_dcr_list_NS_class.o\
	bf_detector_dcr_list_N_class.o\
	bf_detector_dcr_list_S_class.o\
	bf_detector_dcr_list_EW_class.o\
	bf_detector_dcr_list_E_class.o\
	bf_detector_dcr_list_W_class.o\
	bf_path_icr_class.o\
	bf_nbc_template_module.o\
	bf_interface_icr_class.o\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_layer_remove_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_activation_module.o\
	bf_layer_errors_module.o\
	parameters_bf_layer.o\
	bf_layer_remove_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	test_bf_layer_module.o

test_bf_interface_dcr_update_prog:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_activation_module.o\
	bf_detector_module.o\
	bf_detector_dcr_list_class.o\
	bf_detector_dcr_list_NS_class.o\
	parameters_bf_layer.o\
	bf_detector_dcr_list_N_class.o\
	bf_detector_dcr_list_S_class.o\
	bf_detector_dcr_list_EW_class.o\
	bf_detector_dcr_list_E_class.o\
	bf_detector_dcr_list_W_class.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_icr_list_class.o\
	bf_path_icr_class.o\
	bf_nbc_template_module.o\
	bf_interface_icr_class.o\
	bf_interface_dcr_class.o\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_interface_dcr_update_prog2:\
	parameters_kind.o\
	parameters_constant.o\
	parameters_input.o\
	bf_activation_module.o\
	bf_detector_module.o\
	bf_detector_dcr_list_class.o\
	bf_detector_dcr_list_NS_class.o\
	parameters_bf_layer.o\
	bf_detector_dcr_list_N_class.o\
	bf_detector_dcr_list_S_class.o\
	bf_detector_dcr_list_EW_class.o\
	bf_detector_dcr_list_E_class.o\
	bf_detector_dcr_list_W_class.o\
	bf_layer_errors_module.o\
	bf_layer_allocate_module.o\
	bf_layer_reallocate_module.o\
	bf_layer_merge_module.o\
	bf_layer_exchange_module.o\
	bf_layer_remove_module.o\
	bf_layer_class.o\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	bf_mainlayer_pointer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	nbf_interface_class.o\
	bf_interface_class.o\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_icr_list_class.o\
	bf_path_icr_class.o\
	bf_nbc_template_module.o\
	bf_interface_icr_class.o\
	bf_interface_dcr_class.o\
	test_bf_layer_module.o\
	bubble_class.o\
	rising_bubble_class.o\
	growing_bubble_class.o\
	test_case_class.o


#-----------------------------------------------------------------------


#tests
test_field.o:			$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(field_dir)/field_class.o

test_field:			field_class.o\
				parameters_input.o\
				parameters_kind.o\
				surrogate_class.o

test_cg_operators.o:		$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(sd_dir)/cg_operators_class.o

test_cg_operators:		dim2d_prim_module.o\
				field_class.o\
				interface_primary.o\
				parameters_input.o\
				parameters_kind.o\
				cg_operators_class.o\
				surrogate_class.o

test_dim2d_prim.o:		$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_prim_module.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_prim:		dim2d_parameters.o\
				dim2d_prim_module.o\
				field_class.o\
				parameters_input.o\
				parameters_kind.o\
				surrogate_class.o

test_dim2d_fluxes.o:		$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_fluxes_module.o\
				$(field_dir)/field_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_fluxes:		cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				field_class.o\
				parameters_kind.o\
				surrogate_class.o

test_dim2d_eq.o:		$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_dim2d_eq:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_state_eq_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_steadystate_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o

test_fv_operators.o:		$(pbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/dim2d_eq_class.o

test_fv_operators:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o\
				td_operators_class.o\
				fv_operators_class.o

test_bc_periodic.o:		$(pbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_periodic:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropbubble_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_vortex_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o

test_bc_reflection.o:		$(rbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_reflection:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_state_eq_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_steadystate_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				reflection_xy_module.o\
				bc_operators_class.o

test_wall_xy_module.o:		$(wbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(wbc_dir)/wall_xy_module.o

test_wall_xy_module:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropbubble_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_vortex_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				wall_prim_module.o\
				wall_xy_module.o\
				bc_operators_class.o


test_bc_wall.o:			$(wbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_bc_wall:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				wall_prim_module.o\
				wall_xy_module.o\
				bc_operators_class.o

test_rk3tvd.o:			$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(simpletest_dir)/dim2d_eq_class.o\
				$(rk3tvd_dir)/rk3tvd_class.o

test_rk3tvd:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				bc_abstract_class.o\
				bc_operators_class.o\
				td_operators_class.o\
				fv_operators_class.o\
				td_integrator_class.o\
				rk3tvd_steps_class.o\
				rk3tvd_class.o

test_rk3tvd_dim2d.o:		$(sbc_dir)/bc_operators_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(fv_dir)/fv_operators_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_class.o

test_rk3tvd_dim2d:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				periodic_xy_module.o\
				bc_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				td_operators_class.o\
				fv_operators_class.o\
				td_integrator_class.o\
				rk3tvd_class.o


test_nf90_operators.o:		$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(nf90_dir)/nf90_operators_wr_class.o

test_nf90_operators:		surrogate_class.o\
				parameters_kind.o\
				parameters_input.o\
				field_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_state_eq_module.o\
				dim2d_dropbubble_module.o\
				dim2d_steadystate_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				parameters_constant.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				nf90_operators_module.o\
				nf90_operators_wr_class.o\
				io_operators_module.o

test_dim2d_ic.o:		$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(nf90_dir)/nf90_operators_wr_class.o

test_dim2d_ic:			surrogate_class.o\
				parameters_input.o\
				parameters_kind.o\
				field_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o\
				dim2d_drop_retraction_module.o\
				dim2d_bubble_ascending_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				dim2d_state_eq_module.o\
				parameters_constant.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				nf90_operators_module.o\
				nf90_operators_wr_class.o\
				io_operators_module.o



test_field_par.o:		$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_kind.o

test_field_par:			field_par_class.o\
				field_class.o\
				surrogate_class.o\
				mpi_process_class.o\
				parameters_input.o\
				parameters_kind.o

test_mpi_mg_bc.o:		$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_mpi_mg_bc:			surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_process_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o

test_mpi_mg_bc_ext.o:		$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o

test_mpi_mg_bc_ext:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_process_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o

test_reflection_xy_par_module.o:$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rbc_dir)/reflection_xy_par_module.o	

test_reflection_xy_par_module:	parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				reflection_xy_module.o\
				reflection_xy_par_module.o

test_reflection_xy_par.o:	$(rbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_reflection_xy_par:		parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				reflection_xy_module.o\
				reflection_xy_par_module.o\
				bc_abstract_par_class.o\
				bc_operators_par_class.o

test_mpi_mg_construct.o:	$(sd_dir)/cg_operators_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_bc_dir)/mpi_mg_construct.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o


test_mpi_mg_construct:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_mg_construct.o

test_periodic_xy_par_module.o:	$(field_dir)/field_par_class.o\
				$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pbc_dir)/periodic_xy_par_module.o

test_periodic_xy_par_module:	parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o

test_periodic_xy_par.o:		$(pbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_periodic_xy_par:		parameters_kind.o\
				dim2d_parameters.o\
				surrogate_class.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				dim2d_prim_module.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				dim2d_homogeneous_module.o\
				dim2d_phase_separation_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_construct.o\
				mpi_mg_bc_ext_class.o\
				mpi_process_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o\
				bc_abstract_par_class.o\
				bc_operators_par_class.o

test_rk3tvd_dim2d_par.o:	$(sbc_dir)/bc_operators_par_class.o\
				$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(fv_dir)/fv_operators_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(rk3tvd_dir)/rk3tvd_par_class.o

test_rk3tvd_dim2d_par:		surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				dim2d_homogeneous_liquid.o\
				dim2d_phase_separation_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				mpi_process_class.o\
				field_par_class.o\
				mpi_mg_neighbours.o\
				mpi_mg_derived_types.o\
				mpi_mg_bc_class.o\
				mpi_mg_construct.o\
				mpi_mg_ini_bc_proc.o\
				mpi_mg_bc_ext_class.o\
				bc_abstract_par_class.o\
				mpi_tag_module.o\
				periodic_xy_par_module.o\
				bc_operators_par_class.o\
				td_operators_par_class.o\
				fv_operators_par_class.o\
				td_integrator_par_class.o\
				rk3tvd_steps_class.o\
				rk3tvd_par_class.o

test_nf90_operators_par.o:	$(sd_dir)/cg_operators_class.o\
				$(dim2d_dir)/dim2d_eq_class.o\
				$(field_dir)/field_par_class.o\
				$(mpi_dir)/mpi_process_class.o\
				$(nf90_dir)/nf90_operators_wr_par_class.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_nf90_operators_par:	surrogate_class.o\
				parameters_kind.o\
				parameters_constant.o\
				parameters_input.o\
				field_class.o\
				interface_primary.o\
				sd_operators_class.o\
				cg_operators_class.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_fluxes_module.o\
				dim2d_state_eq_module.o\
				dim2d_steadystate_module.o\
				dim2d_dropletbubble_module.o\
				phy_model_eq_class.o\
				dim2d_eq_class.o\
				mpi_process_class.o\
				field_par_class.o\
				io_operators_module.o\
				nf90_operators_module.o\
				nf90_operators_wr_par_class.o
