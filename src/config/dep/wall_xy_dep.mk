$(wbc_dir)/ridders_method_module.o:\
			$(test_dir)/tools/check_data_module.o\
			$(param_dir)/parameters_kind.o

$(wbc_dir)/wall_xy_equilibrium_module.o:\
			$(test_dir)/tools/check_data_module.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(wbc_dir)/ridders_method_module.o\
			$(sd_cdir)/sd_operators_class.o\
			$(sd_cdir)/sd_operators_fd_module.o

$(wsbc_dir)/bc_operators_class.o:\
			$(bc_dir)/bc_operators_default_class.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_prim_module.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(sd_cdir)/sd_operators_fd_module.o\
			$(wbc_dir)/wall_xy_equilibrium_module.o