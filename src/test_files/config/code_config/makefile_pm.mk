ifeq ($(strip $(pm_choice)), simpletest_choice)
	pm_cdir=$(simpletest_dir)
	sim_dep+=$(simpletest_dep)
	sim_par_dep+=$(simpletest_dep)
endif
ifeq ($(strip $(pm_choice)), wave1d_choice)
	pm_cdir=$(wave1d_dir)
	sim_dep+=$(wave1d_dep)
	sim_par_dep+=$(wave1d_dep)
endif
ifeq ($(strip $(pm_choice)), wave2d_choice)
	pm_cdir=$(wave2d_dir)
	sim_dep+=$(wave2d_dep)
	sim_par_dep+=$(wave2d_dep)
endif
ifeq ($(strip $(pm_choice)), ns2d_choice)
	pm_cdir=$(ns2d_dir)
	sim_dep+=$(ns2d_dep)
	sim_par_dep+=$(ns2d_dep)
endif
ifeq ($(strip $(pm_choice)), dim2d_choice)
	pm_cdir=$(dim2d_dir)
	sim_dep+=$(dim2d_dep)
	sim_par_dep+=$(dim2d_dep)
endif