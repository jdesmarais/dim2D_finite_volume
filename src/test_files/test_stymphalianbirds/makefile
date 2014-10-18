#-----------------------------------------------------------------------
#makefile commands
#-----------------------------------------------------------------------
include $(AUGEANSTABLES_CONFIG)/makefile_header.mk


test_bf_compute_prog.o:	$(bc_cdir)/bc_operators_class.o\
			$(cbf_layer_dir)/bf_compute_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o\
			$(pm_cdir)/pmodel_eq_class.o\
			$(sd_cdir)/sd_operators_class.o\
			$(td_cdir)/td_operators_class.o

test_bf_compute_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			bf_compute_class.o\
			bf_remove_module.o\
			bf_layer_errors_module.o

test_bf_layer_bc_procedure.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_procedure:\
			parameters_kind.o\
			parameters_bf_layer.o\
			bf_layer_bc_procedure_module.o


test_bf_layer_bc_sections.o:\
			$(cbf_layer_dir)/bf_layer_bc_procedure_module.o\
			$(cbf_layer_dir)/bf_layer_bc_sections_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o

test_bf_layer_bc_sections:\
			$(param_dep)\
			bf_layer_bc_procedure_module.o\
			bf_layer_bc_sections_class.o


test_bf_layer_prog.o:	$(bf_layer_dir)/parameters_bf_layer.o\
			$(bf_layer_dir)/bf_layer_class.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_layer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)


test_bf_mainlayer_prog.o:$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_mainlayer_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_layer_dep)\
			bf_sublayer_class.o\
			bf_mainlayer_class.o

test_bf_interface_prog.o:$(bf_layer_dir)/bf_interface_class.o\
			$(bf_layer_dir)/bf_mainlayer_class.o\
			$(bf_layer_dir)/bf_sublayer_class.o\
			$(bf_layer_dir)/parameters_bf_layer.o\
			$(param_dir)/parameters_constant.o\
			$(param_dir)/parameters_input.o\
			$(param_dir)/parameters_kind.o

test_bf_interface_prog:	$(param_dep)\
			$(cg_dep)\
			$(simpletest_dep)\
			$(periodic_dep)\
			$(fv_dep)\
			interface_integration_step.o\
			$(bf_interface_dep)

