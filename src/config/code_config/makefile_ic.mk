#Navier-Stokes initial conditions
ifeq ($(strip $(pm_choice)), ns2d_choice)
	ifeq ($(strip $(ic_choice)), steady_state)
		ic_cdir=$(ns2d_sic)
	endif

	ifeq ($(strip $(ic_choice)), peak)
		ic_cdir=$(ns2d_pic)
	endif

	ifeq ($(strip $(ic_choice)), vortex)
		ic_cdir=$(ns2d_vic)
	endif

	ifeq ($(strip $(ic_choice)), sym_x)
		ic_cdir=$(ns2d_sxic)
	endif

	ifeq ($(strip $(ic_choice)), sym_y)
		ic_cdir=$(ns2d_syic)
	endif
endif

#Diffuse interface model initial conditions
ifeq ($(strip $(pm_choice)), dim2d_choice)
	ifeq ($(strip $(ic_choice)), bubble_ascending)
		ic_cdir=$(dim2d_ic)/bubble_ascending

		dim2d_ic_dep+=	dim2d_state_eq_module.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o
	endif

	ifeq ($(strip $(ic_choice)), bubble_transported)
		ic_cdir=$(dim2d_ic)/bubble_transported

		dim2d_ic_dep+=	dim2d_state_eq_module.o\
				dim2d_dropbubble_module.o
	endif

	ifeq ($(strip $(ic_choice)), drop_collision)
		ic_cdir=$(dim2d_ic)/drop_collision

		dim2d_ic_dep+=	dim2d_state_eq_module.o\
				dim2d_dropbubble_module.o\
				dim2d_vortex_module.o
	endif

	ifeq ($(strip $(ic_choice)), drop_retraction)
		ic_cdir=$(dim2d_ic)/drop_retraction

		dim2d_ic_dep+=	dim2d_state_eq_module.o\
				dim2d_dropbubble_module.o
	endif

	ifeq ($(strip $(ic_choice)), homogeneous_liquid)
		ic_cdir=$(dim2d_ic)/homogeneous_liquid

		dim2d_ic_dep+=	dim2d_state_eq_module.o
	endif

	ifeq ($(strip $(ic_choice)), phase_separation)
		ic_cdir=$(dim2d_ic)/phase_separation

		dim2d_ic_dep+=	dim2d_state_eq_module.o
	endif

	ifeq ($(strip $(ic_choice)), steady_state)
		ic_cdir=$(dim2d_ic)/steady_state

		dim2d_ic_dep+=	dim2d_state_eq_module.o
	endif

endif