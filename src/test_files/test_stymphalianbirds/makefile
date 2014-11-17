#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

test_openbc_local_operators_module.o:\
			$(bc_cdir)/bc_operators_class.o\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o\
			$(pm_cdir)/pmodel_eq_class.o

test_hedstrom_xy_bf.o:	test_openbc_local_operators_module.o

test_hedstrom_xy_bf:	$(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)\
			$(hedstrom_xy_dep)\
			test_openbc_local_operators_module.o


test_poinsot_xy_bf.o:	test_openbc_local_operators_module.o

test_poinsot_xy_bf:	$(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)\
			$(poinsot_xy_dep)\
			lodi_ns2d_class.o\
			test_openbc_local_operators_module.o

test_yoolodato_xy_bf.o:	test_openbc_local_operators_module.o

test_yoolodato_xy_bf:	$(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)\
			$(yoo_ns2d_dep)\
			test_openbc_local_operators_module.o


test_pmodel_nopt.o:	$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o

#test_pmodel_nopt:	$(param_dep)\
#			$(mt_dep)\
#			$(simpletest_dep)

#test_pmodel_nopt:	$(param_dep)\
#			$(mt_dep)\
#			$(wave1d_dep)

#test_pmodel_nopt:	$(param_dep)\
#			$(mt_dep)\
#			$(wave2d_dep)

test_pmodel_nopt:	$(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)

#test_pmodel_nopt:	$(param_dep)\
#			$(mt_dep)\
#			$(dim2d_dep)


test_td_operators_nopt.o:\
			$(bc_cdir)/bc_operators_class.o\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(td_cdir)/td_operators_class.o

test_td_operators_nopt: $(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)\
			$(poinsot_ns2d_dep)\
			$(fv_dep)\
			bf_layer_bc_procedure_module.o\
			bf_layer_bc_sections_class.o


test_bf_layer_bc_procedure.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_procedure:\
			parameters_kind.o\
			parameters_bf_layer.o\
			bf_layer_bc_procedure_module.o


test_bf_layer_bc_sections.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(cbf_layer_dir)/bf_layer_bc_sections_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_sections:\
			$(param_dep)\
			bf_layer_bc_procedure_module.o\
			bf_layer_bc_sections_class.o

test_bf_interior_bc_sections.o:\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(cbf_layer_dir)/bf_interior_bc_sections_module.o\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o

test_bf_interior_bc_sections:\
			$(param_dep)\
			bf_interior_bc_sections_module.o\
			bf_layer_bc_procedure_module.o

test_bf_compute.o:	$(bc_cdir)/bc_operators_class.o\
			$(cbf_layer_dir)/bf_compute_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(td_cdir)/td_operators_class.o\
			$(rk3tvd_dir)/rk3tvd_steps_module.o

test_bf_compute:	$(param_dep)\
			$(cg_dep)\
			$(ns2d_dep)\
			$(poinsot_xy_dep)\
			$(fv_dep)\
			lodi_ns2d_class.o\
			interface_integration_step.o\
			bf_layer_bc_sections_class.o\
			bf_layer_bc_procedure_module.o\
			bf_compute_class.o\
			bf_remove_module.o\
			bf_layer_errors_module.o\
			rk3tvd_steps_module.o

test_field_extended.o:	$(field_dir)/field_class.o\
			$(field_dir)/field_extended_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_field_extended:	$(sim_dep)\
			$(bf_interface_dcr_dep)\
			surrogate_class.o\
			field_abstract_class.o\
			field_class.o\
			field_extended_class.o


test_field_extended_timedev.o:\
			$(bf_layer_dir)/bf_interface_class.o\
			$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(field_dir)/field_class.o\
			$(field_dir)/field_extended_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(ns2d_pic)/ic_class.o

test_field_extended_timedev:\
			$(sim_dep)\
			$(bf_interface_dcr_dep)\
			surrogate_class.o\
			field_abstract_class.o\
			field_class.o\
			field_extended_class.o

test_fixed_field_extended.o:\
			$(bf_layer_dir)/bf_interface_class.o\
			$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(field_dir)/field_extended_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_fixed_field_extended:\
			$(sim_dep)\
			$(bf_interface_dcr_dep)\
			surrogate_class.o\
			field_abstract_class.o\
			field_extended_class.o

test_bf_layer_newgrdpt_procedure.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(cbf_layer_dir)/bf_layer_newgrdpt_procedure_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_kind.o

test_bf_layer_newgrdpt_procedure:\
			parameters_kind.o\
			parameters_constant.o\
			parameters_input.o\
			parameters_bf_layer.o\
			bf_layer_newgrdpt_procedure_module.o\
			bf_layer_bc_procedure_module.o

test_sd_operators_fd_n_module.o:\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/sd_operators_fd_n_module.o

test_sd_operators_fd_n_module:\
			parameters_kind.o\
			interface_primary.o\
			sd_operators_fd_module.o\
			sd_operators_fd_n_module.o

test_bf_newgrdpt.o:	$(cbf_layer_dir)/bf_newgrdpt_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wave2d_dir)/pmodel_eq_class.o\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/sd_operators_fd_n_module.o\
			$(wave2d_dir)/wave2d_parameters.o

test_bf_newgrdpt:	$(param_dep)\
			$(mt_dep)\
			$(wave2d_dep)\
			bf_layer_bc_procedure_module.o\
			bf_layer_newgrdpt_procedure_module.o\
			bf_newgrdpt_class.o\
			sd_operators_fd_n_module.o\
			n_coords_module.o
