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