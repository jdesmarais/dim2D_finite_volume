#space discretization operators
$(sd_dir)/interface_primary.o:\
			$(param_dir)/parameters_kind.o

$(sd_dir)/sd_operators_fd_n_module.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_fd_module.o

$(sd_dir)/sd_operators_abstract_class.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o


#cockburn and gau operators
$(cg_dir)/sd_operators_fd_module.o:\

$(cg_dir)/sd_operators_x_oneside_R1_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

$(cg_dir)/sd_operators_x_oneside_R0_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

$(cg_dir)/sd_operators_y_oneside_L0_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

$(cg_dir)/sd_operators_y_oneside_L1_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

$(cg_dir)/sd_operators_y_oneside_R1_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o

$(cg_dir)/sd_operators_y_oneside_R0_class.o:\
			$(cg_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o\
			$(cg_dir)/sd_operators_class.o


#mattsson operators
$(mt_dir)/sd_operators_fd_module.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o

$(mt_dir)/sd_operators_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/sd_operators_abstract_class.o

$(mt_dir)/sd_operators_x_oneside_L0_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_x_oneside_L1_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_x_oneside_R1_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_x_oneside_R0_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_y_oneside_L0_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_y_oneside_L1_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_y_oneside_R1_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_y_oneside_R0_class.o:\
			$(mt_dir)/sd_operators_fd_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(mt_dir)/sd_operators_class.o

$(mt_dir)/sd_operators_fd_ncoords_module.o:\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_kind.o


$(mt_dir)/sd_operators_n_class.o:\
			$(sd_cdir)/sd_operators_fd_ncoords_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(sd_dir)/sd_operators_abstract_class.o

$(mt_dir)/sd_operators_n1_oneside_L0_class.o:\
			$(sd_cdir)/sd_operators_fd_ncoords_module.o\
			$(sd_dir)/interface_primary.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_kind.o\
			$(sd_cdir)/sd_operators_n_class.o