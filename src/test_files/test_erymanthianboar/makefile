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
	$(bf_layer_dir)/bf_sublayer_class.o\
	$(param_dir)/parameters_kind.o

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
	$(pm_cdir)/pmodel_eq_class.o\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_layer_remove_prog.o:\
	$(cbf_layer_dir)/bf_remove_module.o\
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
	$(pm_cdir)/pmodel_eq_class.o\
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

test_bf_layer_nf90_operators_prog.o:\
	$(iobf_layer_dir)/bf_layer_nf90_operators_module.o\
	$(dim2d_dir)/pmodel_eq_class.o\
	$(bf_layer_dir)/parameters_bf_layer.o\
	$(param_dir)/parameters_constant.o\
	$(param_dir)/parameters_kind.o\
	$(param_dir)/parameters_input.o

#-----------------------------------------------------------------------
#dependencies for the test files executables
#-----------------------------------------------------------------------
test_bf_layer_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
	$(bc_dep)\
	test_bf_layer_module.o

test_bf_sublayer_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
	bf_sublayer_class.o

test_bf_mainlayer_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_bf_mainlayer_module.o

test_nbf_element_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
	bf_sublayer_class.o\
	bf_mainlayer_class.o\
	nbf_element_class.o\
	test_bf_layer_module.o

test_nbf_list_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
	bf_sublayer_class.o\
	nbf_element_class.o\
	bf_sublayer_pointer_class.o\
	sbf_list_class.o\
	nbf_list_class.o\
	bf_mainlayer_class.o\
	test_bf_layer_module.o\
	test_nbf_list_module.o

test_nbf_interface_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_layer_dep)\
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
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_interface_dep)\
	test_bf_layer_module.o

test_bf_path_icr_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_interface_dep)\
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
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_interface_dep)\
	dbf_element_class.o\
	dbf_list_class.o\
	bf_detector_module.o\
	bf_detector_icr_list_class.o\
	test_bf_layer_module.o

test_bf_interface_icr_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(simpletest_dep)\
	$(bf_interface_icr_dep)\
	test_bf_layer_module.o

test_bf_interface_icr_update_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(dim2d_dep)\
	$(bf_interface_icr_dep)\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_layer_remove_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(dim2d_dep)\
	$(bf_interface_dep)\
	test_bf_layer_module.o

test_bf_interface_dcr_update_prog:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(dim2d_dep)\
	$(bf_interface_dcr_dep)\
	test_bf_layer_module.o\
	test_cases_interface_update_module.o

test_bf_interface_dcr_update_prog2:\
	$(param_dep)\
	$(cg_dep)\
	$(fv_dep)\
	$(dim2d_dep)\
	$(bf_interface_dcr_dep)\
	test_bf_layer_module.o\
	bubble_class.o\
	rising_bubble_class.o\
	growing_bubble_class.o\
	test_case_class.o

test_bf_layer_nf90_operators_prog:\
	$(param_dep)\
	bf_layer_errors_module.o\
	parameters_bf_layer.o\
	bf_layer_nf90_operators_module.o\
	$(sd_dep)\
	$(dim2d_dep)

