#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

#bf_layer_sync
test_verify_x_symmetry.o:	$(test_dir)/tools/check_data_module.o\
				$(io_dir)/cmd_operators_class.o\
				$(nf90_dir)/nf90_operators_module.o\
				$(nf90_dir)/nf90_operators_read_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o

test_verify_x_symmetry:		$(param_dep)\
				parameters_bf_layer.o\
				$(mt_dep)\
				$(dim2d_dep)\
				cmd_operators_class.o\
				nf90_operators_module.o\
				nf90_operators_read_module.o\
				check_data_module.o