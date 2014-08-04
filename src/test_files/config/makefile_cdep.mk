param_dep=		parameters_kind.o\
			parameters_constant.o\
			parameters_input.o

sd_dep=			interface_primary.o\
			sd_operators_abstract_class.o\
			sd_operators_class.o

cg_dep=			$(sd_dep)\
			sd_operators_fd_module.o

mt_dep=			$(sd_dep)\
			sd_operators_fd_module.o

pm_dep=			pmodel_eq_abstract_class.o\
			pmodel_eq_class.o

simpletest_dep=		$(pm_dep)

wave2d_dep=		wave2d_parameters.o\
			wave2d_prim_module.o\
			$(pm_dep)

dim2d_flux_dep=		dim2d_parameters.o\
			dim2d_prim_module.o\
			dim2d_fluxes_module.o

dim2d_ic_dep=		dim2d_dropbubble_module.o\
			dim2d_vortex_module.o\
			dim2d_state_eq_module.o\
			dim2d_bubble_ascending_module.o\
			dim2d_drop_collision_module.o\
			dim2d_drop_retraction_module.o\
			dim2d_homogeneous_module.o\
			dim2d_phase_separation_module.o\
			dim2d_steadystate_module.o

dim2d_dep=		$(dim2d_flux_dep)\
			$(dim2d_ic_dep)\
			parameters_bf_layer.o\
			$(pm_dep)

bc_dep=			bc_operators_abstract_class.o\
			bc_operators_default_class.o\
			bc_operators_class.o

reflection_dep=		$(bc_dep)\
			reflection_xy_module.o

periodic_dep=		$(bc_dep)

wall_xy_dep=		$(bc_dep)\
			wall_prim_module.o\
			wall_xy_module.o

wall_x_reflection_dep=	$(bc_dep)\
			reflection_xy_module.o\
			wall_prim_module.o\
			wall_xy_module.o

hedstrom_xy_dep=	$(bc_dep)\
			openbc_operators_module.o\
			sd_operators_x_oneside_L0_class.o\
			sd_operators_x_oneside_L1_class.o\
			sd_operators_x_oneside_R1_class.o\
			sd_operators_x_oneside_R0_class.o\
			sd_operators_y_oneside_L0_class.o\
			sd_operators_y_oneside_L1_class.o\
			sd_operators_y_oneside_R1_class.o\
			sd_operators_y_oneside_R0_class.o\
			sd_operators_fd_module.o

td_dep=			td_operators_abstract_class.o\
			td_operators_class.o

fv_dep=			$(td_dep)


td_par_dep=		td_operators_abstract_par_class.o\
			td_operators_par_class.o

fv_par_dep=		$(td_par_dep)


ti_dep=			td_integrator_abstract_class.o\
			td_integrator_class.o

ti_par_dep=		td_integrator_abstract_par_class.o\
			td_integrator_par_class.o

rk_dep=			$(ti_dep)\
			rk3tvd_steps_module.o

rk_par_dep=		$(ti_par_dep)\
			rk3tvd_steps_module.o

io_dep=			io_operators_abstract_class.o\
			io_operators_class.o

nf90_dep=		$(io_dep)\
			io_operators_module.o\
			nf90_operators_module.o

io_par_dep=		io_operators_abstract_par_class.o\
			io_operators_par_class.o

nf90_par_dep=		$(io_par_dep)\
			io_operators_module.o\
			nf90_operators_module.o

mpi_dep=		mpi_process_class.o

mpi_mg_dep=		mpi_mg_neighbours.o\
			mpi_mg_derived_types.o\
			mpi_mg_bc_class.o

mpi_mg_ext_dep=		$(mpi_mg_dep)\
			mpi_mg_construct.o\
			mpi_mg_ini_bc_proc.o\
			mpi_mg_bc_ext_class.o

mpi_bc_dep=		$(mpi_mg_ext_dep)\
			mpi_tag_module.o\
			mpi_requests_module.o

reflection_xy_par_dep=	$(mpi_dep)\
			$(mpi_bc_dep)\
			reflection_xy_module.o\
			reflection_xy_par_module.o\
			bc_operators_abstract_par_class.o\
			bc_operators_par_class.o

periodic_xy_par_dep=	$(mpi_dep)\
			$(mpi_bc_dep)\
			periodic_xy_par_module.o\
			bc_operators_abstract_par_class.o\
			bc_operators_par_class.o

wall_xy_par_dep=	$(mpi_dep)\
			$(mpi_bc_dep)\
			wall_prim_module.o\
			wall_xy_module.o\
			wall_xy_par_module.o\
			bc_operators_abstract_par_class.o\
			bc_operators_par_class.o

wall_x_refl_y_par_dep=	$(mpi_dep)\
			$(mpi_bc_dep)\
			reflection_xy_module.o\
			wall_prim_module.o\
			wall_xy_module.o\
			wall_x_reflection_y_par_module.o\
			wall_xy_par_module.o\
			bc_operators_abstract_par_class.o\
			bc_operators_par_class.o

bf_layer_dep=		parameters_bf_layer.o\
			bf_remove_module.o\
			bf_compute_class.o\
			bf_layer_errors_module.o\
			bf_layer_allocate_module.o\
			bf_layer_reallocate_module.o\
			bf_layer_merge_module.o\
			bf_layer_exchange_module.o\
			bf_layer_nf90_operators_module.o\
			bf_layer_class.o

bf_interface_dep=	$(bf_layer_dep)\
			bf_sublayer_class.o\
			bf_mainlayer_class.o\
			bf_mainlayer_pointer_class.o\
			nbf_element_class.o\
			bf_sublayer_pointer_class.o\
			sbf_list_class.o\
			nbf_list_class.o\
			nbf_interface_class.o\
			bf_interface_class.o

bf_interface_icr_dep=	$(bf_interface_dep)\
			dbf_element_class.o\
			dbf_list_class.o\
			bf_detector_module.o\
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
			bf_interface_icr_class.o

bf_interface_dcr_dep=	$(bf_interface_icr_dep)\
			bf_interface_dcr_class.o
