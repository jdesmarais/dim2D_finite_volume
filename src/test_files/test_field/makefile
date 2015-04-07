#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk

test_field_extended.o:		$(bc_cdir)/bc_operators_class.o\
				$(field_dir)/field_class.o\
				$(field_dir)/field_extended_class.o\
				$(bf_layer_dir)/bf_sublayer_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(sd_cdir)/sd_operators_class.o\
				$(td_cdir)/td_operators_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o

test_field_extended:		$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(hedstrom_xy_dep)\
				$(fv_dep)\
				$(rk_dep)\
				$(nf90_dep)\
				$(field_dep)\
				$(field_extended_dep)\
				check_data_module.o\
				rk3tvd_steps_module.o


test_field_par.o:		$(bc_cdir)/bc_operators_class.o\
				$(field_dir)/field_class.o\
				$(field_dir)/field_par_class.o\
				$(test_dir)/tools/check_data_module.o\
				$(bf_layer_dir)/parameters_bf_layer.o\
				$(param_dir)/parameters_constant.o\
				$(param_dir)/parameters_input.o\
				$(param_dir)/parameters_kind.o\
				$(pm_cdir)/pmodel_eq_class.o\
				$(sd_cdir)/sd_operators_class.o\
				$(td_cdir)/td_operators_class.o\
				$(rk3tvd_dir)/rk3tvd_steps_module.o

test_field_par:			$(param_dep)\
				$(mt_dep)\
				$(wave2d_dep)\
				$(reflection_dep)\
				$(fv_dep)\
				$(rk_dep)\
				$(nf90_dep)\
				$(field_dep)\
				$(field_par_dep)\
				check_data_module.o\
				rk3tvd_steps_module.o