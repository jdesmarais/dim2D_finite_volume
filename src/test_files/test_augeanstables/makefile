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
			interface_integration_step.o\
			bf_compute_class.o


test_bf_layer_prog.o:	$(bf_layer_dir)/parameters_bf_layer.o\
			$(bf_layer_dir)/bf_layer_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_layer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)


test_bf_mainlayer_prog.o:$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_mainlayer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)\
			bf_sublayer_class.o\
			bf_mainlayer_class.o

test_bf_interface_prog.o:$(bf_layer_dir)/bf_interface_class.o\
			$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_interface_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_interface_dep)

test_mattsson_operators.o:\
			$(sd_dir)/interface_primary.o\
			$(sd_dir)/sd_operators_abstract_class.o\
			$(mt_dir)/sd_operators_class.o\
			$(mt_dir)/sd_operators_x_oneside_L0_class.o\
			$(mt_dir)/sd_operators_x_oneside_L1_class.o\
			$(mt_dir)/sd_operators_x_oneside_R1_class.o\
			$(mt_dir)/sd_operators_x_oneside_R0_class.o\
			$(mt_dir)/sd_operators_y_oneside_L0_class.o\
			$(mt_dir)/sd_operators_y_oneside_L1_class.o\
			$(mt_dir)/sd_operators_y_oneside_R1_class.o\
			$(mt_dir)/sd_operators_y_oneside_R0_class.o\
			$(mt_dir)/sd_operators_fd_module.o\
			$(mt_dir)/sd_operators_fd_ncoords_module.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mattsson_operators:\
			$(param_dep)\
			$(mt_dep)\
			sd_operators_x_oneside_L0_class.o\
			sd_operators_x_oneside_L1_class.o\
			sd_operators_x_oneside_R1_class.o\
			sd_operators_x_oneside_R0_class.o\
			sd_operators_y_oneside_L0_class.o\
			sd_operators_y_oneside_L1_class.o\
			sd_operators_y_oneside_R1_class.o\
			sd_operators_y_oneside_R0_class.o\
			interface_primary.o\
			sd_operators_fd_module.o\
			sd_operators_fd_ncoords_module.o\
			dim2d_parameters.o\
			dim2d_prim_module.o

test_cg_operators.o:\
			$(sd_dir)/interface_primary.o\
			$(sd_dir)/sd_operators_abstract_class.o\
			$(cg_dir)/sd_operators_class.o\
			$(cg_dir)/sd_operators_x_oneside_L0_class.o\
			$(cg_dir)/sd_operators_x_oneside_L1_class.o\
			$(cg_dir)/sd_operators_x_oneside_R0_class.o\
			$(cg_dir)/sd_operators_x_oneside_R1_class.o\
			$(cg_dir)/sd_operators_y_oneside_L0_class.o\
			$(cg_dir)/sd_operators_y_oneside_L1_class.o\
			$(cg_dir)/sd_operators_y_oneside_R1_class.o\
			$(cg_dir)/sd_operators_y_oneside_R0_class.o\
			$(cg_dir)/sd_operators_fd_module.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_cg_operators:\
			$(param_dep)\
			$(cg_dep)\
			sd_operators_x_oneside_L0_class.o\
			sd_operators_x_oneside_L1_class.o\
			sd_operators_x_oneside_R0_class.o\
			sd_operators_x_oneside_R1_class.o\
			sd_operators_y_oneside_L0_class.o\
			sd_operators_y_oneside_L1_class.o\
			sd_operators_y_oneside_R1_class.o\
			sd_operators_y_oneside_R0_class.o\
			interface_primary.o\
			sd_operators_fd_module.o\
			dim2d_parameters.o\
			dim2d_prim_module.o

test_openbc_operators.o:\
			$(obc_dir)/openbc_operators_module.o\
			$(simpletest_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o\
			$(mt_dir)/sd_operators_x_oneside_L0_class.o\
			$(mt_dir)/sd_operators_x_oneside_L1_class.o\
			$(mt_dir)/sd_operators_x_oneside_R1_class.o\
			$(mt_dir)/sd_operators_x_oneside_R0_class.o\
			$(mt_dir)/sd_operators_y_oneside_L0_class.o\
			$(mt_dir)/sd_operators_y_oneside_L1_class.o\
			$(mt_dir)/sd_operators_y_oneside_R1_class.o\
			$(mt_dir)/sd_operators_y_oneside_R0_class.o

test_openbc_operators:	$(param_dep)\
			$(mt_dep)\
			$(simpletest_dep)\
			sd_operators_x_oneside_L0_class.o\
			sd_operators_x_oneside_L1_class.o\
			sd_operators_x_oneside_R0_class.o\
			sd_operators_x_oneside_R1_class.o\
			sd_operators_y_oneside_L0_class.o\
			sd_operators_y_oneside_L1_class.o\
			sd_operators_y_oneside_R1_class.o\
			sd_operators_y_oneside_R0_class.o\
			openbc_operators_module.o

test_hedstrom_xy.o:	$(hrobc_dir)/bc_operators_class.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(obc_dir)/openbc_operators_module.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o\
			$(wave1d_dir)/pmodel_eq_class.o\
			$(mt_dir)/sd_operators_fd_module.o

test_hedstrom_xy:	$(param_dep)\
			$(mt_dep)\
			$(wave1d_dep)\
			$(hedstrom_x_reflection_y_dep)

test_hedstrom_xy_corners.o:\
			$(hcobc_dir)/bc_operators_class.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o\
			$(wave2d_dir)/pmodel_eq_class.o\
			$(mt_dir)/sd_operators_fd_ncoords_module.o

test_hedstrom_xy_corners:\
			$(param_dep)\
			$(mt_dep)\
			$(wave2d_dep)\
			$(hedstrom_xy_corners_dep)

test_ns2d_eq_program.o:	$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/ns2d_parameters.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(ns2d_dir)/ns2d_fluxes_module.o\
			$(mt_dir)/sd_operators_class.o

test_ns2d_eq_program:	parameters_kind.o\
			ns2d_parameters.o\
			interface_primary.o\
			ns2d_prim_module.o\
			ns2d_fluxes_module.o\
			$(mt_dep)

#			openbc_operators_module.o\
#			sd_operators_abstract_class.o\
#			sd_operators_class.o\
#			sd_operators_x_oneside_L0_class.o\
#			sd_operators_x_oneside_L1_class.o\
#			sd_operators_x_oneside_R1_class.o\
#			sd_operators_x_oneside_R0_class.o\
#			sd_operators_y_oneside_L0_class.o\
#			sd_operators_y_oneside_L1_class.o\
#			sd_operators_y_oneside_R1_class.o\
#			sd_operators_y_oneside_R0_class.o\
#			sd_operators_fd_module.o\
#			hedstrom_xy_module.o