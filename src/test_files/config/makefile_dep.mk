#parameter_dir						
$(param_dir)/parameters_kind.o:

$(param_dir)/parameters_constant.o:\
			$(param_dir)/parameters_kind.o

$(param_dir)/parameters_input.o:\
			$(param_dir)/parameters_constant.o


#fields
$(field_dir)/surrogate_class.o:

$(field_dir)/field_abstract_class.o:\
			$(bc_cdir)/bc_operators_class.o\
			$(field_dir)/surrogate_class.o\
			$(ti_dir)/interface_integration_step.o\
			$(io_dir)/cmd_operators_class.o\
			$(io_cdir)/io_operators_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(td_cdir)/td_operators_class.o

$(field_dir)/field_abstract_par_class.o:\
			$(bc_cdir)/bc_operators_par_class.o\
			$(ti_dir)/interface_integration_step.o\
			$(io_cdir)/io_operators_par_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(field_dir)/surrogate_class.o\
			$(td_cdir)/td_operators_par_class.o

$(field_dir)/field_class.o:\
			$(field_dir)/field_abstract_class.o\
			$(param_dir)/parameters_kind.o\
			$(ti_cdir)/td_integrator_class.o

$(field_dir)/field_extended_class.o:\
			$(bf_layer_dir)/bf_interface_dcr_class.o\
			$(field_dir)/field_abstract_class.o\
			$(ti_dir)/interface_integration_step.o\
			$(param_dir)/parameters_kind.o\
			$(ti_cdir)/td_integrator_class.o

$(field_dir)/field_par_class.o:\
			$(field_dir)/field_abstract_par_class.o\
			$(param_dir)/parameters_kind.o\
			$(ti_cdir)/td_integrator_par_class.o


#physical models
$(phy_eq_dir)/pmodel_eq_abstract_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_class.o

$(phy_eq_dir)/pmodel_eq_default_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/pmodel_eq_abstract_class.o

$(phy_eq_dir)/ic_abstract_class.o:\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

#simple test equations
$(simpletest_dir)/pmodel_eq_class.o:\
			$(sd_dir)/interface_primary.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/pmodel_eq_default_class.o\
			$(sd_cdir)/sd_operators_class.o

#wave1d equations
$(wave1d_dir)/wave1d_parameters.o:\
			$(param_dir)/parameters_kind.o

$(wave1d_dir)/wave1d_prim_module.o:\
			$(param_dir)/parameters_kind.o

$(wave1d_dir)/pmodel_eq_class.o:\
			$(sd_dir)/interface_primary.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/pmodel_eq_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(wave1d_dir)/wave1d_parameters.o\
			$(wave1d_dir)/wave1d_prim_module.o

#wave2d equations
$(wave2d_dir)/wave2d_parameters.o:\
			$(param_dir)/parameters_kind.o

$(wave2d_dir)/wave2d_prim_module.o:\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/n_coords_module.o

$(wave2d_dir)/wave2d_ncoords_module.o:\
			$(sd_dir)/interface_primary.o\
			$(wave2d_dir)/wave2d_parameters.o\
			$(wave2d_dir)/wave2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(wave2d_dir)/pmodel_eq_class.o:\
			$(sd_dir)/interface_primary.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/pmodel_eq_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(wave2d_dir)/wave2d_parameters.o\
			$(wave2d_dir)/wave2d_prim_module.o


#navier-stokes equations
$(ns2d_dir)/ns2d_parameters.o:\
			$(param_dir)/parameters_kind.o

$(ns2d_dir)/ns2d_prim_module.o:\
			$(ns2d_dir)/ns2d_parameters.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(ns2d_dir)/ns2d_ncoords_module.o:\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(param_dir)/parameters_kind.o

$(ns2d_dir)/ns2d_fluxes_module.o:\
			$(ns2d_dir)/ns2d_parameters.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_class.o

#$(ns2d_ic)/ns2d_steadystate_module.o:\
#			$(param_dir)/parameters_kind.o
#
#$(ns2d_ic)/ns2d_vortex_module.o:\
#			$(ns2d_dir)/ns2d_parameters.o\
#			$(param_dir)/parameters_kind.o
#
#$(ns2d_ic)/ns2d_peak_module.o:\
#			$(ns2d_dir)/ns2d_parameters.o\
#			$(param_dir)/parameters_kind.o

$(ns2d_sic)/ic_class.o:	$(phy_eq_dir)/ic_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o

$(ns2d_pic)/ic_class.o:	$(ns2d_dir)/ns2d_parameters.o\
			$(phy_eq_dir)/ic_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(ns2d_vic)/ic_class.o:	$(ns2d_dir)/ns2d_parameters.o\
			$(phy_eq_dir)/ic_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(ns2d_sxic)/ic_class.o:$(ns2d_dir)/ns2d_parameters.o\
			$(phy_eq_dir)/ic_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(ns2d_syic)/ic_class.o:$(ns2d_dir)/ns2d_parameters.o\
			$(phy_eq_dir)/ic_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o


$(ns2d_dir)/pmodel_eq_class.o:\
			$(ic_cdir)/ic_class.o\
			$(sd_dir)/interface_primary.o\
			$(sd_cdir)/sd_operators_class.o\
			$(ns2d_dir)/ns2d_ncoords_module.o\
			$(ns2d_dir)/ns2d_parameters.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(ns2d_dir)/ns2d_fluxes_module.o\
			$(ns2d_ic)/ns2d_steadystate_module.o\
			$(ns2d_ic)/ns2d_vortex_module.o\
			$(ns2d_ic)/ns2d_peak_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(phy_eq_dir)/pmodel_eq_default_class.o




#boundary conditions
$(bc_dir)/bc_operators_abstract_class.o:\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o

$(bc_dir)/bc_operators_default_class.o:\
			$(bc_dir)/bc_operators_abstract_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o

$(bc_dir)/bc_operators_abstract_par_class.o:\
			$(sd_cdir)/sd_operators_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(bc_dir)/bc_operators_nopt_module.o:\
			$(cbf_layer_dir)/bf_layer_sync_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o


#periodic boundary conditions
$(pbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(pbc_dir)/periodic_xy_par_module.o:\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(pbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_operators_abstract_par_class.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pbc_dir)/periodic_xy_par_module.o

#reflection boundary conditions
$(rbc_dir)/reflection_xy_module.o:\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(rbc_dir)/reflection_xy_par_module.o:\
			$(pm_cdir)/pmodel_eq_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(rbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o

$(rbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_operators_abstract_par_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(rbc_dir)/reflection_xy_par_module.o

#wall boundary conditions
$(wbc_dir)/wall_prim_module.o:\
			$(dim2d_dir)/dim2d_parameters.o\
			$(param_dir)/parameters_kind.o

$(wbc_dir)/wall_xy_module.o:\
			$(sd_cdir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_prim_module.o

$(wbc_dir)/wall_xy_par_module.o:\
			$(sd_cdir)/sd_operators_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

$(wbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

$(wbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_operators_abstract_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(wbc_dir)/wall_xy_module.o\
			$(wbc_dir)/wall_xy_par_module.o

#wall_x reflection_y boundary conditions
$(wrbc_dir)/wall_x_reflection_y_par_module.o:\
			$(sd_cdir)/sd_operators_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

$(wrbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(wbc_dir)/wall_xy_module.o

$(wrbc_dir)/bc_operators_par_class.o:\
			$(bc_dir)/bc_operators_abstract_par_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(sd_cdir)/sd_operators_class.o\
			$(wrbc_dir)/wall_x_reflection_y_par_module.o\
			$(wbc_dir)/wall_xy_module.o\
			$(wbc_dir)/wall_xy_par_module.o

#open boundary conditions
include $(AUGEANSTABLES_CONFIG)/dep/bf_layer_dep.mk

$(obc_dir)/bc_operators_openbc_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(bc_dir)/bc_operators_nopt_module.o\
			$(bf_layer_dir)/bf_layer_errors_module.o\
			$(sd_dir)/interface_primary.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

$(obc_dir)/bc_operators_openbc_normal_class.o:\
			$(obc_dir)/bc_operators_openbc_class.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

$(obc_dir)/openbc_operators_module.o:\
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



#hedstrom_xy_corners open boundary conditions
$(hcobc_dir)/hedstrom_xy_corners_module.o:\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(hcobc_dir)/bc_operators_class.o:\
			$(obc_dir)/bc_operators_openbc_normal_class.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(hcobc_dir)/hedstrom_xy_corners_module.o\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_fd_ncoords_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

#hedstrom_x_reflection_y boundary conditions
$(hrobc_dir)/bc_operators_class.o:\
			$(obc_dir)/bc_operators_openbc_normal_class.o\
			$(hobc_dir)/hedstrom_xy_module.o\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

#poinsot open boundary conditions
$(pobc_dir)/lodi_abstract_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(pobc_dir)/lodi_xy_module.o:\
			$(sd_dir)/interface_primary.o\
			$(pobc_cdir)/lodi_inflow_class.o\
			$(pobc_cdir)/lodi_outflow_class.o\
			$(obc_dir)/openbc_operators_module.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_fd_module.o

$(pobc_dir)/bc_operators_class.o:\
			$(obc_dir)/bc_operators_openbc_normal_class.o\
			$(pobc_dir)/lodi_xy_module.o\
			$(pobc_cdir)/lodi_inflow_class.o\
			$(pobc_cdir)/lodi_outflow_class.o\
			$(sd_dir)/interface_primary.o\
			$(obc_dir)/openbc_operators_module.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(sd_cdir)/sd_operators_x_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_x_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_x_oneside_R0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L0_class.o\
			$(sd_cdir)/sd_operators_y_oneside_L1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R1_class.o\
			$(sd_cdir)/sd_operators_y_oneside_R0_class.o

#poinsot boundary conditions for NS equations
$(pobc_dir)/lodi_class.o:\
			$(sd_dir)/interface_primary.o\
			$(pobc_dir)/lodi_abstract_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(pobc_ns2d_dir)/lodi_inflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(pobc_dir)/lodi_class.o\
			$(ns2d_dir)/ns2d_parameters.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/pmodel_eq_class.o

$(pobc_ns2d_dir)/lodi_outflow_class.o:\
			$(sd_dir)/interface_primary.o\
			$(pobc_dir)/lodi_class.o\
			$(ns2d_dir)/ns2d_parameters.o\
			$(ns2d_dir)/ns2d_prim_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/pmodel_eq_class.o

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




#time discretization methods
$(td_dir)/td_operators_abstract_class.o:\
			$(bc_cdir)/bc_operators_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(td_dir)/td_operators_par_class.o:\
			$(sbc_dir)/bc_operators_par_class.o\
			$(sd_dir)/cg_operators_class.o\
			$(field_dir)/field_par_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(spm_dir)/dim2d_eq_class.o

$(fv_dir)/td_operators_class.o:\
			$(cbf_layer_dir)/bf_layer_bc_sections_class.o\
			$(bc_cdir)/bc_operators_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(td_dir)/td_operators_abstract_class.o

$(fv_dir)/td_operators_par_class.o:\
			$(bc_cdir)/bc_operators_par_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(td_dir)/td_operators_abstract_par_class.o


#time integration methods
$(ti_dir)/interface_integration_step.o:\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(ti_dir)/td_integrator_abstract_class.o:\
			$(field_dir)/field_abstract_class.o\
			$(param_dir)/parameters_kind.o

$(ti_dir)/td_integrator_abstract_par_class.o:\
			$(field_dir)/field_abstract_par_class.o\
			$(param_dir)/parameters_kind.o

$(rk3tvd_dir)/td_integrator_class.o:\
			$(field_dir)/field_abstract_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rk3tvd_dir)/rk3tvd_steps_module.o\
			$(ti_dir)/td_integrator_abstract_class.o

$(rk3tvd_dir)/td_integrator_par_class.o:\
			$(field_dir)/field_abstract_par_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rk3tvd_dir)/rk3tvd_steps_module.o\
			$(ti_dir)/td_integrator_abstract_par_class.o


#io operators
$(io_dir)/cmd_operators_class.o:\
			$(bf_layer_dir)/parameters_bf_layer.o

$(io_dir)/io_operators_module.o:

$(io_dir)/io_operators_abstract_class.o:\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(io_dir)/io_operators_abstract_par_class.o:\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(nf90_dir)/nf90_operators_module.o:\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(nf90_dir)/nf90_operators_read_module.o:\
			$(nf90_dir)/nf90_operators_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(nf90_dir)/io_operators_class.o:\
			$(io_dir)/io_operators_module.o\
			$(io_dir)/io_operators_abstract_class.o\
			$(nf90_dir)/nf90_operators_module.o\
			$(nf90_dir)/nf90_operators_read_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o

$(nf90_dir)/io_operators_par_class.o:\
			$(io_dir)/io_operators_module.o\
			$(io_dir)/io_operators_abstract_par_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(nf90_dir)/nf90_operators_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o


#mpi processes
$(mpi_dir)/mpi_tag_module.o:

$(mpi_dir)/mpi_requests_module.o:\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_tag_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(mpi_dir)/mpi_process_class.o:\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o


#mpi messenger bc
$(mpi_bc_dir)/mpi_mg_neighbours.o:\
			$(param_dir)/parameters_constant.o

$(mpi_bc_dir)/mpi_mg_derived_types.o:\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

$(mpi_bc_dir)/mpi_mg_bc_class.o:\
			$(mpi_bc_dir)/mpi_mg_neighbours.o\
			$(mpi_bc_dir)/mpi_mg_derived_types.o\
			$(param_dir)/parameters_input.o

$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o:\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o

$(mpi_bc_dir)/mpi_mg_construct.o:\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o

$(mpi_bc_dir)/mpi_mg_bc_ext_class.o:\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(mpi_bc_dir)/mpi_mg_construct.o\
			$(mpi_bc_dir)/mpi_mg_ini_bc_proc.o\
			$(param_dir)/parameters_input.o


#simulations
sim_dim2d.o:		$(field_dir)/field_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

sim_dim2d_bf.o:		$(field_dir)/field_extended_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

sim_dim2d_par.o:	$(field_dir)/field_par_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o


#dependencies for the executable code
sim_dim2d:		$(sim_dep)\
			cmd_operators_class.o\
			surrogate_class.o\
			field_abstract_class.o\
			field_class.o\

sim_dim2d_bf:		$(sim_dep)\
			$(bf_interface_dcr_dep)\
			cmd_operators_class.o\
			surrogate_class.o\
			field_abstract_class.o\
			field_extended_class.o

sim_dim2d_par:		$(sim_par_dep)\
			surrogate_class.o\
			field_abstract_par_class.o\
			field_par_class.o

#include $(sim_dim2d_dep)
#include $(sim_dim2d_par_dep)

