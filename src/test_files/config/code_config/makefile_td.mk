#finite volume time discretization
ifeq ($(strip $(td_choice)), finitevolume_choice)
	td_cdir=$(fv_dir)
	sim_dep+=$(fv_dep)
	sim_par_dep+=$(fv_par_dep)
endif
