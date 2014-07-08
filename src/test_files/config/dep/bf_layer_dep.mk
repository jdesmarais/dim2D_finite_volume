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

$(bf_layer_dir)/bf_layer_remove_module.o:\
	$(bf_layer_dir)/bf_activation_module.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_layer_class.o:\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_layer_allocate_module.o\
	$(bf_layer_dir)/bf_layer_reallocate_module.o\
	$(bf_layer_dir)/bf_layer_merge_module.o\
	$(bf_layer_dir)/bf_layer_exchange_module.o\
	$(bf_layer_dir)/bf_layer_remove_module.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

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
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_element_class.o\
	$(param_dir)/parameters_kind.o\
	$(sbf_layer_dir)/sbf_list_class.o

$(nbf_layer_dir)/nbf_interface_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(nbf_layer_dir)/nbf_list_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(sbf_layer_dir)/sbf_list_class.o

$(bf_layer_dir)/bf_mainlayer_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_mainlayer_pointer_class.o:\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_class.o:\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_mainlayer_pointer_class.o\
	$(nbf_layer_dir)/nbf_interface_class.o\
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
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/dbf_element_class.o:\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/dbf_list_class.o:\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_module.o:\
	$(param_dir)/parameters_kind.o

$(dbf_layer_dir)/bf_detector_icr_list_class.o:\
	$(dbf_layer_dir)/bf_detector_module.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(dbf_layer_dir)/dbf_element_class.o\
	$(dbf_layer_dir)/dbf_list_class.o\
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

$(bf_layer_dir)/bf_activation_module.o:\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_icr_class.o:\
	$(bf_layer_dir)/bf_activation_module.o\
	$(dbf_layer_dir)/bf_detector_dcr_list_class.o\
	$(dbf_layer_dir)/bf_detector_dcr_list_N_class.o\
	$(dbf_layer_dir)/bf_detector_dcr_list_S_class.o\
	$(dbf_layer_dir)/bf_detector_dcr_list_E_class.o\
	$(dbf_layer_dir)/bf_detector_dcr_list_W_class.o\
	$(dbf_layer_dir)/bf_detector_icr_list_class.o\
	$(bf_layer_dir)/bf_layer_errors_module.o\
	$(bf_layer_dir)/bf_path_icr_class.o\
	$(bf_layer_dir)/bf_nbc_template_module.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(bf_layer_dir)/bf_interface_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o

$(bf_layer_dir)/bf_interface_dcr_class.o:\
	$(bf_layer_dir)/bf_interface_icr_class.o\
	$(bf_layer_dir)/bf_mainlayer_class.o\
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_input.o\
	$(param_dir)/parameters_kind.o\
	$(sbf_layer_dir)/sbf_list_class.o