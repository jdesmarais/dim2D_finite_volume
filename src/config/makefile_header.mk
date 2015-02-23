#-----------------------------------------------------------------------
#cluster
#-----------------------------------------------------------------------
cluster = $(AUGEANSTABLES_CLUSTER)


#-----------------------------------------------------------------------
#user main options for profiling, optimizing ...
#-----------------------------------------------------------------------
profile       = $(AUGEANSTABLES_PROFILE)        #ddt/totalview debugging
analyse       = $(AUGEANSTABLES_ANALYSE)        #gprof
trace         = $(AUGEANSTABLES_TRACE)          #vampir tracing
report        = $(AUGEANSTABLES_OPTREPORT)      #optimization report
profileGuided = $(AUGEANSTABLES_PROFILE_GUIDED) #profile guided
scalasca      = false			        #scalasca instrumentation


#-----------------------------------------------------------------------
#user main options
#-----------------------------------------------------------------------
src_dir	   = $(augeanstables)/src
config_dir = $(AUGEANSTABLES_CONFIG)
dep_dir	   = $(AUGEANSTABLES_CONFIG)/dep

sd_choice = mt_choice                    #space discretization choice
pm_choice = dim2d_choice                #physical model choice
ic_choice = bubble_transported                         #initial conditions choice
bc_choice = hedstrom_xy_choice           #boundary condition choice
td_choice = finitevolume_choice          #time discretization choice
ti_choice = rk3tvd_choice                #time integration choice
io_choice = nf90_choice                  #writer choice

sd_cdep=
pm_cdep=

#-----------------------------------------------------------------------
#source files directories
#-----------------------------------------------------------------------
#define the folder paths for the source code compilation
#-----------------------------------------------------------------------
include $(config_dir)/makefile_folders.mk
include $(config_dir)/makefile_cdep.mk


#-----------------------------------------------------------------------
#folder options for compilation
# e.g. depending on the boundary condition used, the folder selected
#      for the compilation changes
#-----------------------------------------------------------------------
sim_dep     = $(param_dep)
sim_par_dep = $(param_dep)

include $(config_dir)/code_config/makefile_sd.mk  #space discretisation
include $(config_dir)/code_config/makefile_pm.mk  #physical model
include $(config_dir)/code_config/makefile_ic.mk  #initial conditions
include $(config_dir)/code_config/makefile_bc.mk  #boundary conditions
include $(config_dir)/code_config/makefile_td.mk  #time discretization
include $(config_dir)/code_config/makefile_ti.mk  #time integration
include $(config_dir)/code_config/makefile_io.mk  #writer choice


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

#path for the mpi, hdf5, netcdf libraries depending on the cluster used
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

include $(config_dir)/dep/tools_dep.mk		#tools for tests
include $(config_dir)/dep/sd_dep.mk         	#space discretization operators
include $(config_dir)/dep/dim2d_dep.mk      	#dim2D physical model
include $(config_dir)/dep/hedstrom_xy_dep.mk	#hedstrom_xy b.c.



