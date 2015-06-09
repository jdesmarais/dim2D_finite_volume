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
			$(bf_layer_dir)/parameters_bf_layer.o\
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
