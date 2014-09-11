bc_dir			= $(src_dir)/boundary_conditions
field_dir      	        = $(src_dir)/field
io_dir			= $(src_dir)/io_operators
mpi_dir			= $(src_dir)/mpi_operators
param_dir    		= $(src_dir)/parameters
phy_eq_dir        	= $(src_dir)/physical_models
sd_dir			= $(src_dir)/sd_operators
ti_dir			= $(src_dir)/td_integrator
td_dir        		= $(src_dir)/td_operators
test_dir                = $(src_dir)/test_files

config_dir		= $(test_dir)/config

cg_dir			= $(sd_dir)/cg_operators
mt_dir			= $(sd_dir)/mattsson_operators

rbc_dir		       	= $(bc_dir)/reflection_xy
pbc_dir		        = $(bc_dir)/periodic_xy
wbc_dir			= $(bc_dir)/wall_xy
wrbc_dir		= $(bc_dir)/wall_x_reflection_y
obc_dir                 = $(bc_dir)/open_bc

nf90_dir		= $(io_dir)/nf90_operators

mpi_bc_dir 		= $(mpi_dir)/mpi_messenger_bc

simpletest_dir          = $(phy_eq_dir)/simpletest
wave1d_dir              = $(phy_eq_dir)/wave1d
wave2d_dir              = $(phy_eq_dir)/wave2d
ns2d_dir		= $(phy_eq_dir)/ns2d
dim2d_dir               = $(phy_eq_dir)/dim2d


rk3tvd_dir		= $(ti_dir)/rungekutta3rdtvd

fv_dir			= $(td_dir)/finitevolume

ns2d_ic			= $(ns2d_dir)/ns2d_ic
dim2d_ic		= $(dim2d_dir)/dim2d_ic

ns2d_sic		= $(ns2d_ic)/steady_state
ns2d_pic		= $(ns2d_ic)/peak
ns2d_vic		= $(ns2d_ic)/vortex
ns2d_sxic		= $(ns2d_ic)/sym_x
ns2d_syic		= $(ns2d_ic)/sym_y


hobc_dir                = $(obc_dir)/hedstrom_xy
hcobc_dir               = $(obc_dir)/hedstrom_xy_corners
hrobc_dir		= $(obc_dir)/hedstrom_x_reflection_y

pobc_dir		= $(obc_dir)/poinsot_xy
pobc_ns2d_dir		= $(pobc_dir)/poinsot_ns2d

yobc_dir		= $(obc_dir)/yoo_lodato_xy
yobc_ns2d_dir		= $(yobc_dir)/yoo_lodato_ns2d

bf_layer_dir		= $(obc_dir)/buffer_layer
nbf_layer_dir		= $(bf_layer_dir)/bf_layer_neighbors
dbf_layer_dir           = $(bf_layer_dir)/bf_layer_detectors
sbf_layer_dir           = $(bf_layer_dir)/bf_sublayer_lists
iobf_layer_dir		= $(bf_layer_dir)/bf_layer_io_operators
cbf_layer_dir		= $(bf_layer_dir)/bf_layer_computations