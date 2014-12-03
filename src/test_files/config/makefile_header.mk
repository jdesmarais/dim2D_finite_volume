#-----------------------------------------------------------------------
#cluster
#-----------------------------------------------------------------------
cluster = $(AUGEANSTABLES_CLUSTER)


#-----------------------------------------------------------------------
#user main options for profiling, optimizing ...
#-----------------------------------------------------------------------
profile       = $(AUGEANSTABLES_PROFILE)       #ddt/totalview debugging
analyse       = $(AUGEANSTABLES_ANALYSE)       #gprof
trace         = $(AUGEANSTABLES_TRACE)         #vampir tracing
report        = $(AUGEANSTABLES_OPTREPORT)     #optimization report
profileGuided = $(AUGEANSTABLES_PROFILE_GUIDED)#profile guided
scalasca      = false			       #scalasca instrumentation


#-----------------------------------------------------------------------
#user main options
#-----------------------------------------------------------------------
src_dir	   = $(augeanstables)/src
config_dir = $(AUGEANSTABLES_CONFIG)
dep_dir	   = $(AUGEANSTABLES_CONFIG)/dep

sd_choice = mt_choice            #space discretization choice
pm_choice = simpletest_choice 	 #physical model choice
ic_choice = vortex               #initial conditions choice
bc_choice = hedstrom_xy_choice   #boundary condition choice
td_choice = finitevolume_choice  #time discretization choice
ti_choice = rk3tvd_choice        #time integration choice
io_choice = nf90_choice          #writer choice

#-----------------------------------------------------------------------
#source files directories
#-----------------------------------------------------------------------
#define the folder paths for the source code compilation
#-----------------------------------------------------------------------
include $(config_dir)/makefile_folders.mk
include $(config_dir)/makefile_cdep.mk


#-----------------------------------------------------------------------
#folder options for compilation
#-----------------------------------------------------------------------
sim_dep     = $(param_dep)
sim_par_dep = $(param_dep)

#space discretisation
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

#physical model
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

#initial conditions
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

#boundary conditions
ifeq ($(strip $(bc_choice)), periodic_xy_choice)
	bc_cdir=$(pbc_dir)
	sim_dep+=$(periodic_dep)
	sim_par_dep+=$(periodic_xy_par_dep)
endif
ifeq ($(strip $(bc_choice)), reflection_xy_choice)
	bc_cdir=$(rbc_dir)
	sim_dep+=$(reflection_dep)
	sim_par_dep+=$(reflection_xy_par_dep)
endif
ifeq ($(strip $(bc_choice)), wall_xy_choice)
	bc_cdir=$(wbc_dir)
	sim_dep+=$(wall_xy_dep)
	sim_par_dep+=$(wall_xy_par_dep)
endif
ifeq ($(strip $(bc_choice)), wall_x_reflection_y_choice)
	bc_cdir=$(wrbc_dir)
	sim_dep+=$(wall_x_reflection_dep)
	sim_par_dep+=$(wall_x_refl_y_par_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_xy_choice)
	bc_cdir=$(hobc_dir)
	sim_dep+=$(hedstrom_xy_dep)
	sim_par_dep+=$(hedstrom_xy_par_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_xy_corners_choice)
	bc_cdir=$(hcobc_dir)
	sim_dep+=$(hedstrom_xy_corners_dep)
	sim_par_dep+=$(hedstrom_xy_corners_par_dep)
endif
ifeq ($(strip $(bc_choice)), hedstrom_x_reflection_y_choice)
	bc_cdir=$(hrobc_dir)
	sim_dep+=$(hedstrom_x_reflection_y_dep)
	sim_par_dep+=$(hedstrom_x_reflection_y_par_dep)
endif
ifeq ($(strip $(bc_choice)), poinsot_xy_choice)
	bc_cdir=$(pobc_dir)
	sim_dep+=$(poinsot_ns2d_dep)
	sim_par_dep+=$(poinsot_xy_par_dep)

	ifeq ($(strip $(pm_choice)), ns2d_choice)
		pobc_cdir=$(pobc_ns2d_dir)
	endif
endif
ifeq ($(strip $(bc_choice)), yoolodato_xy_choice)
	bc_cdir=$(yobc_dir)

	ifeq ($(strip $(pm_choice)), ns2d_choice)

		sim_dep+=$(yoo_ns2d_dep)
		sim_par_dep+=$(yoolodato_xy_par_dep)

		yobc_cdir=$(yobc_ns2d_dir)
	endif
endif

#time discretization
ifeq ($(strip $(td_choice)), finitevolume_choice)
	td_cdir=$(fv_dir)
	sim_dep+=$(fv_dep)
	sim_par_dep+=$(fv_par_dep)
endif

#time integration
ifeq ($(strip $(ti_choice)), rk3tvd_choice)
	ti_cdir=$(rk3tvd_dir)
	sim_dep+=$(rk_dep)
	sim_par_dep+=$(rk_par_dep)
endif

#writer choice
ifeq ($(strip $(io_choice)), nf90_choice)
	io_cdir=$(nf90_dir)
	sim_dep+=$(nf90_dep)
	sim_par_dep+=$(nf90_par_dep)
endif

#-----------------------------------------------------------------------
#compiler and options
#-----------------------------------------------------------------------
PREP	=
FC	= $(AUGEANSTABLES_COMPILER)
FFLAGS	= 
LDFLAGS =
INCLUDE =
LIBS	=

#compiler options------------------------------------------------
include $(config_dir)/compiler_config/options_ifort.mk
include $(config_dir)/compiler_config/options_gfortran.mk

#path for the mpi, hdf5, netcdf libraries depending on the cluster used-
include $(config_dir)/libraries_config/libraries_bolt.mk
include $(config_dir)/libraries_config/libraries_cartesius.mk
include $(config_dir)/libraries_config/libraries_newcluster.mk

#options handeling depending on the user main options: debug, trace,
#and the type of compiler used (ifort, gfortran, mpiifort...)
include $(config_dir)/compiler_config/options_choice.mk


#-----------------------------------------------------------------------
#source code files of the current dir
#-----------------------------------------------------------------------
SRC	= $(wildcard *.f)	#fortran source files
SRCS	= $(SRC)		#fortran program files
OBJS	= $(SRC:.f=.o)		#objects (.o and .mod)
CMD	= $(SRCS:.f=)		#executables


#----------------------------------------------------------------------
#general rules
#----------------------------------------------------------------------
include $(config_dir)/makefile_rules.mk


#----------------------------------------------------------------------
#main code dependencies (w/o the test files)
#----------------------------------------------------------------------
include $(config_dir)/makefile_dep.mk
