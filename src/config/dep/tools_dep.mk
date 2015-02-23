#check data module
$(test_dir)/tools/check_data_module.o:\
	$(param_dir)/parameters_kind.o

#verify x symmetry
test_verify_x_symmetry.o:\
			$(io_dir)/cmd_operators_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o\
			$(nf90_dir)/nf90_operators_module.o\
			$(nf90_dir)/nf90_operators_read_module.o\
			$(pm_cdir)/pmodel_eq_class.o

test_verify_x_symmetry:	cmd_operators_class.o\
			$(param_dep)\
			$(sd_cdep)\
			$(pm_cdep)\
			nf90_operators_module.o\
			nf90_operators_read_module.o

cmd_verify_symmetry_class.o:\
			$(bf_layer_dir)/parameters_bf_layer.o

#verify y symmetry
test_verify_y_symmetry.o:\
			cmd_verify_symmetry_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(param_dir)/parameters_input.o\
			$(nf90_dir)/nf90_operators_module.o\
			$(nf90_dir)/nf90_operators_read_module.o\
			$(pm_cdir)/pmodel_eq_class.o

test_verify_y_symmetry:	cmd_verify_symmetry_class.o\
			$(param_dep)\
			parameters_bf_layer.o\
			$(sd_cdep)\
			$(pm_cdep)\
			nf90_operators_module.o\
			nf90_operators_read_module.o

#verify x symmetry for detectors
test_dct_verify_x_symmetry.o:\
			$(bf_layer_dir)/bf_restart_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(io_dir)/cmd_operators_class.o\
			$(param_dir)/parameters_kind.o

test_dct_verify_x_symmetry:\
			$(param_dep)\
			parameters_bf_layer.o\
			cmd_operators_class.o\
			bf_restart_module.o

#verify y symmetry for detectors
test_dct_verify_y_symmetry.o:\
			$(bf_layer_dir)/bf_restart_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(io_dir)/cmd_operators_class.o\
			$(param_dir)/parameters_kind.o\
			$(nf90_dir)/nf90_operators_module.o\
			$(nf90_dir)/nf90_operators_read_module.o\
			$(pm_cdir)/pmodel_eq_class.o

test_dct_verify_y_symmetry:\
			$(param_dep)\
			$(sd_cdep)\
			$(pm_cdep)\
			parameters_bf_layer.o\
			cmd_operators_class.o\
			bf_restart_module.o\
			nf90_operators_module.o\
			nf90_operators_read_module.o

