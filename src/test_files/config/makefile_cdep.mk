param_dep=		parameters_kind.o\
			parameters_constant.o\
			parameters_input.o

sd_dep=			interface_primary.o\
			sd_operators_abstract_class.o\
			sd_operators_class.o

pm_dep=			pmodel_eq_abstract_class.o\
			pmodel_eq_class.o

simpletest_dep=		$(pm_dep)

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


td_dep=			td_operators_abstract_class.o\
			td_operators_class.o

fv_dep=			$(td_dep)


ti_dep=			td_integrator_abstract_class.o\
			td_integrator_class.o

rk_dep=			$(ti_dep)


io_dep=			io_operators_abstract_class.o\
			io_operators_class.o

nf90_dep=		$(io_dep)\
			io_operators_module.o\
			nf90_operators_module.o


mpi_dep=		mpi_process_class.o

mpi_mg_dep=		mpi_mg_neighbours.o\
			mpi_mg_derived_types.o\
			mpi_mg_bc_class.o
