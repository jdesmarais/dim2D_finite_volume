#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

test_geometry_update.o:		$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(field_dir)/field_extended_class.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o

test_geometry_update:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(rk_dep)\
				$(nf90_dep)\
				$(field_extended_dep)\
				check_data_module.o