#netcdf i/o
ifeq ($(strip $(io_choice)), nf90_choice)
	io_cdir=$(nf90_dir)
	sim_dep+=$(nf90_dep)
	sim_par_dep+=$(nf90_par_dep)
endif