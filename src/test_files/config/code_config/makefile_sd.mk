ifeq ($(strip $(sd_choice)), cg_choice)
	sd_cdir=$(cg_dir)
	sim_dep+=$(cg_dep)
	sim_par_dep+=$(cg_dep)
endif
ifeq ($(strip $(sd_choice)), mt_choice)
	sd_cdir=$(mt_dir)
	sim_dep+=$(mt_dep)
	sim_par_dep+=$(mt_dep)
endif