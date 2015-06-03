ifeq ($(strip $(bc_choice)), periodic_xy_choice)
	bc_cdir=$(pbc_dir)
	sim_dep+=$(periodic_dep)
	sim_par_dep+=$(periodic_dep)
endif
ifeq ($(strip $(bc_choice)), reflection_xy_choice)
	bc_cdir=$(rbc_dir)
	sim_dep+=$(reflection_dep)
	sim_par_dep+=$(reflection_dep)
endif
ifeq ($(strip $(bc_choice)), wall_xy_choice)
	bc_cdir=$(wbc_dir)
	sim_dep+=$(wall_xy_dep)
	sim_par_dep+=$(wall_xy_dep)
endif
ifeq ($(strip $(bc_choice)), wall_x_reflection_y_choice)
	bc_cdir=$(wrbc_dir)
	sim_dep+=$(wall_x_reflection_dep)
	sim_par_dep+=$(wall_x_reflection_dep)
endif
ifeq ($(strip $(bc_choice)), wall_x_simplified_choice)
	bc_cdir=$(wsbc_dir)
	sim_dep+=$(wall_x_simplified_dep)
	sim_par_dep+=$(wall_x_simplified_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_xy_choice)
	bc_cdir=$(hobc_dir)
	sim_dep+=$(hedstrom_xy_dep)
	sim_par_dep+=$(hedstrom_xy_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_xy_corners_choice)
	bc_cdir=$(hcobc_dir)
	sim_dep+=$(hedstrom_xy_corners_dep)
	sim_par_dep+=$(hedstrom_xy_corners_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_x_reflection_y_choice)
	bc_cdir=$(hrobc_dir)
	sim_dep+=$(hedstrom_x_reflection_y_dep)
	sim_par_dep+=$(hedstrom_x_reflection_y_dep)
endif
ifeq ($(strip $(bc_choice)), poinsot_xy_choice)
	bc_cdir=$(pobc_dir)
	sim_dep+=$(poinsot_ns2d_dep)
	sim_par_dep+=$(poinsot_ns2d_dep)

	ifeq ($(strip $(pm_choice)), ns2d_choice)
		pobc_cdir=$(pobc_ns2d_dir)
	endif
endif
ifeq ($(strip $(bc_choice)), yoolodato_xy_choice)
	bc_cdir=$(yobc_dir)

	ifeq ($(strip $(pm_choice)), ns2d_choice)

		sim_dep+=$(yoo_ns2d_dep)
		sim_par_dep+=$(yoolodato_ns2d_dep)

		yobc_cdir=$(yobc_ns2d_dir)
	endif
endif
