#hedstrom_xy open boundary conditions
$(hobc_dir)/hedstrom_xy_module.o:\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_fd_module.o

$(hobc_dir)/hedstrom_xy_anti_corner_flux_module.o:\
			$(bbf_layer_dir)/bf_layer_bc_checks_module.o\
			$(bbf_layer_dir)/bf_layer_bc_fluxes_module.o\
			$(bbf_layer_dir)/bf_layer_bc_sections_overlap_module.o\
			$(sbf_layer_dir)/bf_layer_extract_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o

$(hobc_dir)/hedstrom_xy_anti_corner_diag_flux_module.o:\
			$(bbf_layer_dir)/bf_layer_bc_checks_module.o\
			$(bbf_layer_dir)/bf_layer_bc_fluxes_module.o\
			$(bbf_layer_dir)/bf_layer_bc_sections_overlap_module.o\
			$(sbf_layer_dir)/bf_layer_extract_module.o\
			$(bbf_layer_dir)/bf_layer_bc_sections_overlap_module.o\
			$(hcobc_dir)/hedstrom_xy_corners_module.o\
			$(sd_dir)/interface_primary.o\
			$(sd_dir)/n_coords_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(sd_cdir)/sd_operators_fd_ncoords_module.o\
			$(sd_cdir)/sd_operators_n1_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_n1_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_n1_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_n1_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_n2_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_n2_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_n2_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_n2_oneside_R0_class.o

$(hobc_dir)/bc_operators_hedstrom_xy_class.o:\
			$(obc_dir)/bc_operators_openbc_normal_class.o\
			$(bc_dir)/bc_operators_nopt_module.o\
			$(bf_layer_dir)/bf_layer_errors_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(hobc_dir)/hedstrom_xy_anti_corner_flux_module.o\
			$(hobc_dir)/hedstrom_xy_anti_corner_diag_flux_module.o\
			$(sd_dir)/interface_primary.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

$(hobc_dir)/bc_operators_class.o:\
			$(hobc_dir)/bc_operators_hedstrom_xy_class.o\
			$(pm_cdir)/pmodel_eq_class.o