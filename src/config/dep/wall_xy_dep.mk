$(wbc_dir)/ridders_method_module.o:\
			$(test_dir)/tools/check_data_module.o\
			$(param_dir)/parameters_kind.o

$(wbc_dir)/wall_xy_equilibrium_module.o:\
			$(test_dir)/tools/check_data_module.o\
			$(dim2d_dir)/dim2d_parameters.o\
			$(dim2d_dir)/dim2d_state_eq_module.o\
			$(wbc_dir)/ridders_method_module.o\
			$(param_dir)/parameters_kind.o