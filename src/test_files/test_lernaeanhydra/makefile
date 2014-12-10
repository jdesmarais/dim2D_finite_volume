#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_cg_operators.o:	$(dim2d_dir)/dim2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

test_cg_operators:	$(param_dep)\
			$(sd_dep)\
			dim2d_prim_module.o


test_dim2d_prim.o:	$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_fd_module.o

test_dim2d_prim:	$(param_dep)\
			dim2d_parameters.o\
			dim2d_prim_module.o\
			sd_operators_fd_module.o


test_dim2d_fluxes.o:	$(cg_dir)/sd_operators_class.o\
			$(dim2d_dir)/dim2d_fluxes_module.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_dim2d_fluxes:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_flux_dep)


test_dim2d_eq.o:	$(cg_dir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_dim2d_eq:		$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)


test_fv_operators.o:	$(pbc_dir)/bc_operators_class.o\
			$(cg_dir)/sd_operators_class.o\
			$(fv_dir)/td_operators_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(simpletest_dir)/pmodel_eq_class.o

test_fv_operators:	$(param_dep)\
			$(sd_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(td_dep)


test_bc_periodic.o:	$(pbc_dir)/bc_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bc_periodic:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(periodic_dep)


test_bc_reflection.o:	$(rbc_dir)/bc_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bc_reflection:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(reflection_dep)


test_wall_xy_module.o:	$(wbc_dir)/bc_operators_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/wall_xy_module.o

test_wall_xy_module:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(wall_xy_dep)


test_bc_wall.o:		$(wbc_dir)/bc_operators_class.o\
			$(cg_dir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bc_wall:		$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(wall_xy_dep)


test_rk3tvd.o:		$(field_dir)/field_abstract_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rk3tvd_dir)/td_integrator_class.o

test_rk3tvd:		$(param_dep)\
			$(sd_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			$(rk_dep)\
			field_abstract_class.o


test_rk3tvd_dim2d.o:	$(field_dir)/field_abstract_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rk3tvd_dir)/td_integrator_class.o

test_rk3tvd_dim2d:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			$(nf90_dep)\
			$(rk_dep)\
			field_abstract_class.o


test_nf90_operators.o:	$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(nf90_dir)/io_operators_class.o

test_nf90_operators:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(nf90_dep)

test_dim2d_ic.o:	$(dim2d_dir)/pmodel_eq_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(nf90_dir)/io_operators_class.o

test_dim2d_ic:		$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(nf90_dep)

#test_field_par.o:	$(field_dir)/field_par_class.o\
#			$(mpi_dir)/mpi_process_class.o\
#			$(param_dir)/parameters_kind.o
#
#test_field_par:		field_par_class.o\
#			field_class.o\
#			surrogate_class.o\
#			mpi_process_class.o\
#			parameters_input.o\
#			parameters_kind.o

test_mpi_mg_bc.o:	$(mpi_dir)/mpi_process_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_mg_bc:		$(param_dep)\
			$(mpi_dep)\
			$(mpi_mg_dep)

test_mpi_mg_bc_ext.o:	$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o

test_mpi_mg_bc_ext:	$(param_dep)\
			$(mpi_dep)\
			$(mpi_mg_ext_dep)

test_reflection_xy_par_module.o:\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rbc_dir)/reflection_xy_module.o\
			$(rbc_dir)/reflection_xy_par_module.o	

test_reflection_xy_par_module:\
			$(param_dep)\
			$(mpi_dep)\
			$(mpi_bc_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			reflection_xy_module.o\
			reflection_xy_par_module.o

test_reflection_xy_par.o:\
			$(rbc_dir)/bc_operators_par_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_reflection_xy_par:	\
			$(param_dep)\
			$(mpi_dep)\
			$(mpi_bc_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(reflection_xy_par_dep)

test_mpi_mg_construct.o:$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_bc_dir)/mpi_mg_construct.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_mpi_mg_construct:	$(param_dep)\
			$(mpi_dep)\
			$(mpi_mg_ext_dep)

test_periodic_xy_par_module.o:\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pbc_dir)/periodic_xy_par_module.o

test_periodic_xy_par_module:\
			$(param_dep)\
			$(mpi_dep)\
			$(mpi_bc_dep)\
			periodic_xy_par_module.o


test_periodic_xy_par.o:	$(pbc_dir)/bc_operators_par_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_periodic_xy_par:	$(param_dep)\
			$(mpi_dep)\
			$(mpi_bc_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(periodic_xy_par_dep)

test_rk3tvd_dim2d_par.o:$(field_dir)/field_abstract_par_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(rk3tvd_dir)/td_integrator_par_class.o

test_rk3tvd_dim2d_par:	$(param_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(periodic_xy_par_dep)\
			$(td_par_dep)\
			$(nf90_par_dep)\
			field_abstract_par_class.o\
			$(rk_par_dep)\


test_nf90_operators_par.o:\
			$(dim2d_dir)/pmodel_eq_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(nf90_dir)/io_operators_par_class.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_nf90_operators_par:\
			$(param_dep)\
			$(mpi_dep)\
			$(sd_dep)\
			$(dim2d_dep)\
			$(nf90_par_dep)


test_ns2d_eq.o:		$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(ns2d_dir)/pmodel_eq_class.o

test_ns2d_eq:		$(param_dep)\
			$(mt_dep)\
			$(ns2d_dep)