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

$(mpi_dir)/mpi_interface_class.o:\
			$(mpi_bc_dir)/mpi_mg_bc_ext_class.o\
			$(mpi_dir)/mpi_process_class.o\
			$(mpi_dir)/mpi_requests_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o


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