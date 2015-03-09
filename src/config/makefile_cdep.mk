param_dep=		parameters_kind.o\
			parameters_constant.o\
			parameters_input.o

sd_dep=			interface_primary.o\
			sd_operators_abstract_class.o\
			sd_operators_class.o

sd_oneside_dep=		sd_operators_x_oneside_L0_class.o\
			sd_operators_x_oneside_L1_class.o\
			sd_operators_x_oneside_R1_class.o\
			sd_operators_x_oneside_R0_class.o\
			sd_operators_y_oneside_L0_class.o\
			sd_operators_y_oneside_L1_class.o\
			sd_operators_y_oneside_R1_class.o\
			sd_operators_y_oneside_R0_class.o

sd_oneside_n_dep=	sd_operators_n_class.o\
			sd_operators_fd_ncoords_module.o\
			sd_operators_n1_oneside_L0_class.o\
			sd_operators_n1_oneside_L1_class.o\
			sd_operators_n1_oneside_R1_class.o\
			sd_operators_n1_oneside_R0_class.o\
			sd_operators_n2_oneside_L0_class.o\
			sd_operators_n2_oneside_L1_class.o\
			sd_operators_n2_oneside_R1_class.o\
			sd_operators_n2_oneside_R0_class.o

cg_dep=			$(sd_dep)\
			sd_operators_fd_module.o

mt_dep=			$(sd_dep)\
			sd_operators_fd_module.o

mt_ext_dep=		$(mt_dep)\
			$(sd_oneside_dep)

pm_dep=			pmodel_eq_abstract_class.o\
			pmodel_eq_default_class.o\
			pmodel_eq_class.o

simpletest_dep=		$(pm_dep)

wave1d_dep=		wave1d_parameters.o\
			wave1d_prim_module.o\
			$(pm_dep)

wave2d_dep=		n_coords_module.o\
			wave2d_parameters.o\
			wave2d_prim_module.o\
			$(pm_dep)

ns2d_dep=		ic_abstract_class.o\
			ic_class.o\
			ns2d_parameters.o\
			ns2d_prim_module.o\
			ns2d_ncoords_module.o\
			ns2d_fluxes_module.o\
			ns2d_steadystate_module.o\
			ns2d_peak_module.o\
			ns2d_vortex_module.o\
			$(pm_dep)

ns_vdw2d_dep=		ns_vdw2d_prim_module.o

dim2d_flux_dep=		dim2d_parameters.o\
			dim2d_prim_module.o\
			dim2d_fluxes_module.o

dim2d_ic_dep=		ic_abstract_class.o\
			ic_class.o

dim2d_dep=		$(ns_vdw2d_dep)\
			$(dim2d_flux_dep)\
			$(dim2d_ic_dep)\
			n_coords_module.o\
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
			bf_layer_errors_module.o\
			bf_layer_extract_module.o\
			bf_layer_bc_sections_overlap_module.o\
			bf_layer_bc_checks_module.o\
			bf_layer_bc_anticorner_module.o\
			bc_operators_openbc_class.o\
			openbc_operators_module.o\
			bc_operators_openbc_normal_class.o\
			sd_operators_fd_module.o\
			$(sd_oneside_dep)\
			$(sd_oneside_n_dep)\
			hedstrom_xy_module.o\
			hedstrom_xy_corners_module.o\
			hedstrom_xy_anti_corner_flux_module.o\
			hedstrom_xy_anti_corner_diag_flux_module.o

hedstrom_xy_corners_dep=\
			$(bc_dep)\
			bc_operators_nopt_module.o\
			bc_operators_openbc_class.o\
			openbc_operators_module.o\
			bc_operators_openbc_normal_class.o\
			$(sd_oneside_dep)\
			sd_operators_fd_module.o\
			hedstrom_xy_module.o\
			hedstrom_xy_corners_module.o

hedstrom_x_reflection_y_dep=\
			$(hedstrom_xy_dep)\
			reflection_xy_module.o

poinsot_xy_dep=		$(bc_dep)\
			bc_operators_nopt_module.o\
			bc_operators_openbc_class.o\
			openbc_operators_module.o\
			bc_operators_openbc_normal_class.o\
			$(sd_oneside_dep)\
			sd_operators_fd_module.o\
			lodi_abstract_class.o\
			lodi_inflow_class.o\
			lodi_outflow_class.o\
			lodi_xy_module.o

poinsot_ns2d_dep=	$(poinsot_xy_dep)\
			lodi_ns2d_class.o


yoo_xy_dep=		$(bc_dep)\
			bc_operators_nopt_module.o\
			bc_operators_openbc_class.o\
			openbc_operators_module.o\
			$(sd_oneside_dep)\
			sd_operators_fd_module.o\
			lodi_edge_abstract_class.o\
			lodi_corner_abstract_class.o\
			lodi_transverse_module.o\
			lodi_component_module.o\
			lodi_timedev_xy_module.o

yoo_ns2d_dep=		$(yoo_xy_dep)\
			lodi_relaxation_coeff_module.o\
			lodi_edge_class.o\
			lodi_edge_inflow_class.o\
			lodi_edge_outflow_class.o\
			lodi_corner_class.o\
			lodi_corner_inflow_inflow_class.o\
			lodi_corner_inflow_outflow_class.o\
			lodi_corner_outflow_outflow_class.o

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
			nf90_operators_module.o\
			nf90_operators_read_module.o

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

bf_newgrdpt_dep=	bf_newgrdpt_class.o\
			n_coords_module.o


#bf_layer objects
bf_layer_basic_dep=	bf_layer_extract_module.o\
			bf_layer_basic_class.o\
			bf_layer_errors_module.o

bf_layer_print_dep=	$(bf_layer_basic_dep)\
			bf_layer_nf90_operators_module.o\
			bf_layer_print_class.o

bf_layer_grdpts_id_update_dep=	$(bf_layer_print_dep)\
				bf_suspicious_bc_interior_pt_module.o\
				bf_bc_crenel_module.o\
				bf_layer_grdpts_id_update_class.o

bf_layer_sync_dep=	$(bf_layer_grdpts_id_update_dep)\
			bf_layer_exchange_module.o\
			bf_layer_sync_class.o

bf_layer_dyn_dep=	$(bf_layer_sync_dep)\
			bf_layer_allocate_module.o\
			bf_layer_reallocate_module.o\
			bf_layer_merge_module.o\
			bf_layer_dyn_class.o

bf_compute_dep=		bf_newgrdpt_procedure_module.o\
			bf_newgrdpt_verification_module.o\
			bf_newgrdpt_class.o\
			bf_compute_class.o\
			bf_compute_newgrdpt_class.o

bf_bc_sections_dep=	bf_layer_bc_procedure_module.o\
			bf_layer_bc_sections_overlap_module.o\
			bf_layer_bc_sections_class.o

bf_layer_time_dep=	$(bf_layer_dyn_dep)\
			$(bf_compute_dep)\
			$(bf_bc_sections_dep)\
			bf_layer_time_class.o

bf_layer_newgrdpt_dep=	$(bf_layer_time_dep)\
			bf_layer_newgrdpt_class.o

bf_layer_dep=		$(bf_layer_newgrdpt_dep)\
			bf_remove_module.o\
			bf_layer_class.o

bf_sublayer_dep=	$(bf_layer_dep)\
			bf_sublayer_class.o


#mainlayer_interface objects
mainlayer_interface_basic_dep=	$(bf_sublayer_dep)\
				mainlayer_interface_basic_class.o

mainlayer_interface_sync_dep=	$(mainlayer_interface_basic_dep)\
				mainlayer_interface_sync_class.o

mainlayer_interface_dyn_dep=	$(mainlayer_interface_sync_dep)\
				mainlayer_interface_dyn_class.o


#bf_mainlayer objects
bf_mainlayer_basic_dep= $(bf_sublayer_dep)\
			bf_mainlayer_basic_class.o

bf_mainlayer_print_dep=	$(bf_mainlayer_basic_dep)\
			bf_mainlayer_print_class.o

bf_mainlayer_sync_dep=	$(bf_mainlayer_print_dep)\
			bf_interior_bc_sections_module.o\
			bf_mainlayer_sync_class.o

bf_mainlayer_dep=	$(bf_mainlayer_sync_dep)\
			bf_mainlayer_class.o

bf_mainlayer_pointer_dep=$(bf_mainlayer_dep)\
			bf_mainlayer_pointer_class.o


#bf_interface objects
bf_interface_basic_dep=	$(bf_mainlayer_pointer_dep)\
			$(mainlayer_interface_dyn_dep)\
			bf_interface_basic_class.o

bf_interface_print_dep= $(bf_interface_basic_dep)\
			bf_interface_print_class.o

bf_interface_sync_dep=	$(bf_interface_print_dep)\
			bf_interface_sync_class.o


bf_interface_dep=	$(bf_layer_dep)\
			bf_restart_module.o\
			bf_sublayer_class.o\
			bf_mainlayer_class.o\
			bf_mainlayer_pointer_class.o\
			bf_interior_bc_sections_module.o\
			nbf_element_class.o\
			bf_sublayer_pointer_class.o\
			sbf_list_class.o\
			nbf_list_class.o\
			nbf_interface_class.o\
			nbf_interface_newgrdpt_class.o\
			bf_interface_class.o\
			nf90_operators_module.o\
			nf90_operators_read_module.o

bf_interface_icr_dep=	$(bf_interface_dep)\
			bf_detector_module.o\
			bf_detector_icr_list_class.o\
			bf_detector_dcr_param_class.o\
			bf_path_icr_class.o\
			bf_nbc_template_module.o\
			bf_interface_icr_class.o

bf_interface_dcr_dep=	$(bf_interface_icr_dep)\
			bf_interface_dcr_class.o
