#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_openbc_operators.o:$(test_dir)/tools/check_data_module.o\
			$(obc_dir)/openbc_operators_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

test_openbc_operators:	$(param_dep)\
			$(mt_ext_dep)\
			$(dim2d_dep)\
			openbc_operators_module.o\
			check_data_module.o

#hedstrom_xy
test_hedstrom_xy.o:	$(test_dir)/tools/check_data_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(mt_dir)/sd_operators_fd_module.o

test_hedstrom_xy:	check_data_module.o\
			$(param_dep)\
			$(mt_ext_dep)\
			$(dim2d_dep)\
			openbc_operators_module.o\
			hedstrom_xy_module.o

test_hedstrom_xy_anti_corner_flux.o:\
			$(test_dir)/tools/check_data_module.o\
			$(hobc_dir)/hedstrom_xy_anti_corner_flux_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(obc_dir)/openbc_operators_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o

test_hedstrom_xy_anti_corner_flux:\
			check_data_module.o\
			$(param_dep)\
			$(mt_ext_dep)\
			$(dim2d_dep)\
			openbc_operators_module.o\
			hedstrom_xy_module.o\
			bf_layer_errors_module.o\
			bf_layer_extract_module.o\
			bf_layer_bc_sections_overlap_module.o\
			bf_layer_bc_checks_module.o\
			bf_layer_bc_fluxes_module.o\
			hedstrom_xy_anti_corner_flux_module.o


test_bc_operators_openbc.o:\
			$(hobc_dir)/bc_operators_class.o\
			$(test_dir)/tools/check_data_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(obc_dir)/openbc_operators_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_fd_module.o

test_bc_operators_openbc:\
			check_data_module.o\
			$(param_dep)\
			$(mt_ext_dep)\
			$(dim2d_dep)\
			$(hedstrom_xy_dep)

test_bc_operators_openbc_normal.o:\
			$(hobc_dir)/bc_operators_class.o\
			$(test_dir)/tools/check_data_module.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(obc_dir)/openbc_operators_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

test_bc_operators_openbc_normal:\
			check_data_module.o\
			$(param_dep)\
			$(mt_ext_dep)\
			$(dim2d_dep)\
			$(hedstrom_xy_dep)


# reflection_xy
test_reflection_xy.o:	$(rbc_dir)/bc_operators_reflection_xy_class.o\
			$(test_dir)/tools/check_data_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(wave2d_dir)/pmodel_eq_class.o

test_reflection_xy:	check_data_module.o\
			$(param_dep)\
			$(mt_dep)\
			$(wave2d_dep)\
			$(reflection_dep)


# periodic_xy
test_periodic_xy.o:	$(pbc_dir)/bc_operators_class.o\
			$(test_dir)/tools/check_data_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(wave2d_dir)/pmodel_eq_class.o

test_periodic_xy:	check_data_module.o\
			$(param_dep)\
			$(mt_dep)\
			$(wave2d_dep)\
			$(periodic_dep)

# wall_xy
ridders_method_fcts_module.o:	$(wbc_dir)/ridders_method_module.o\
				$(param_dir)/parameters_kind.o

test_ridders_method.o:		$(test_dir)/tools/check_data_module.o\
				$(wbc_dir)/ridders_method_module.o\
				ridders_method_fcts_module.o\
				$(param_dir)/parameters_kind.o

test_ridders_method:		check_data_module.o\
				ridders_method_module.o\
				ridders_method_fcts_module.o\
				parameters_kind.o

test_wall_xy_equilibrium.o:	$(test_dir)/tools/check_data_module.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(dim2d_dir)/dim2d_state_eq_module.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_kind.o\
				ridders_method_fcts_module.o\
				$(mt_dir)/sd_operators_class.o\
				$(mt_dir)/sd_operators_fd_module.o\
				$(wbc_dir)/wall_xy_equilibrium_module.o

test_wall_xy_equilibrium:	$(param_dep)\
				$(mt_dep)\
				n_coords_module.o\
				check_data_module.o\
				dim2d_parameters.o\
				dim2d_prim_module.o\
				dim2d_state_eq_module.o\
				ridders_method_module.o\
				wall_xy_equilibrium_module.o\
				ridders_method_fcts_module.o

test_wall_xy.o:			$(wbc_dir)/bc_operators_wall_xy_class.o\
				$(dim2d_dir)/dim2d_parameters.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(dim2d_dir)/pmodel_eq_class.o

test_wall_xy:			$(param_dep)\
				$(mt_dep)\
				sd_operators_x_oneside_L0_class.o\
				sd_operators_x_oneside_L1_class.o\
				sd_operators_x_oneside_R1_class.o\
				sd_operators_x_oneside_R0_class.o\
				n_coords_module.o\
				check_data_module.o\
				$(wall_xy_dep)\
				$(dim2d_dep)
