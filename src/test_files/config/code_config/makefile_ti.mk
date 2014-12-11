#runge-kutta 3rd order TVD
ifeq ($(strip $(ti_choice)), rk3tvd_choice)
	ti_cdir=$(rk3tvd_dir)
	sim_dep+=$(rk_dep)
	sim_par_dep+=$(rk_par_dep)
endif
