#yoo and lodato boundary conditions
$(yobc_dir)/bc_operators_class.o:\
			$(obc_dir)/bc_operators_openbc_class.o\
			$(sd_dir)/interface_primary.o\
			$(yobc_cdir)/lodi_edge_inflow_class.o\
			$(yobc_cdir)/lodi_edge_outflow_class.o\
			$(yobc_cdir)/lodi_corner_inflow_inflow_class.o\
			$(yobc_cdir)/lodi_corner_inflow_outflow_class.o\
			$(yobc_cdir)/lodi_corner_outflow_outflow_class.o\
			$(yobc_dir)/lodi_timedev_xy_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o


$(yobc_dir)/lodi_component_module.o:\
			$(param_dir)/parameters_constant.o

$(yobc_dir)/lodi_transverse_module.o:\
			$(sd_cdir)/sd_operators_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(yobc_dir)/lodi_timedev_xy_module.o:\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(yobc_cdir)/lodi_edge_inflow_class.o\
			$(yobc_cdir)/lodi_edge_outflow_class.o\
			$(yobc_cdir)/lodi_corner_inflow_inflow_class.o\
			$(yobc_cdir)/lodi_corner_inflow_outflow_class.o\
			$(yobc_cdir)/lodi_corner_outflow_outflow_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_dir)/lodi_edge_abstract_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_dir)/lodi_edge_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_edge_abstract_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_dir)/lodi_corner_abstract_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_dir)/lodi_corner_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_corner_abstract_class.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o


#Yoo-Lodato b.c. for the NS equations in 2D
#------------------------------------------
$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o:\
			$(ns2d_dir)/ns2d_parameters.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o

$(yobc_ns2d_dir)/lodi_edge_inflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_edge_class.o\
			$(yobc_dir)/lodi_component_module.o\
			$(yobc_dir)/lodi_transverse_module.o\
			$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/pmodel_eq_class.o

$(yobc_ns2d_dir)/lodi_edge_outflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_edge_class.o\
			$(yobc_dir)/lodi_component_module.o\
			$(yobc_dir)/lodi_transverse_module.o\
			$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/pmodel_eq_class.o

$(yobc_ns2d_dir)/lodi_corner_inflow_inflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_component_module.o\
			$(yobc_dir)/lodi_corner_class.o\
			$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_ns2d_dir)/lodi_corner_inflow_outflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_component_module.o\
			$(yobc_dir)/lodi_corner_class.o\
			$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(yobc_ns2d_dir)/lodi_corner_outflow_outflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(yobc_dir)/lodi_component_module.o\
			$(yobc_dir)/lodi_corner_class.o\
			$(yobc_ns2d_dir)/lodi_relaxation_coeff_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o