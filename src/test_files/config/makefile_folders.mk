#for buffer layer
bf_layer_dir		= $(src_dir)/buffer_layer
nbf_layer_dir		= $(bf_layer_dir)/bf_layer_neighbors
dbf_layer_dir           = $(bf_layer_dir)/bf_layer_detectors
sbf_layer_dir           = $(bf_layer_dir)/bf_sublayer_lists

#main folders
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


#boundary_conditions folder
rbc_dir		       	= $(bc_dir)/reflection_xy
pbc_dir		        = $(bc_dir)/periodic_xy
wbc_dir			= $(bc_dir)/wall_xy
wrbc_dir		= $(bc_dir)/wall_x_reflection_y

#io_operator folder
nf90_dir		= $(io_dir)/nf90_operators

#mpi_operators folder
mpi_bc_dir		= $(mpi_dir)/mpi_messenger_bc			

#phy_model_equations folder
dim2d_dir               = $(phy_eq_dir)/dim2d
simpletest_dir          = $(phy_eq_dir)/simpletest

#td_integrator folder
rk3tvd_dir		= $(ti_dir)/rungekutta3rdtvd

#td_operators folder
fv_dir			= $(td_dir)/finitevolume

#test_files folder
dep_dir 		= $(test_dir)/dep

#dim2d_folder
dim2d_ic		= $(dim2d_dir)/dim2d_ic

##unused folders
#sim_dir			= $(src_dir)/simulations
#tile_operators_dir      = $(src_dir)/tile_operators
#
#io_operator_dir    	= $(io_interface_dir)/io_operators
#restart_op_dir    	= $(io_interface_dir)/restart_operators
#cmd_operators_dir	= $(io_interface_dir)/cmd_operators
#
#mpi_comm_dir		= $(mpi_operators_dir)/mpi_comm
#mpi_pr_dir              = $(mpi_operators_dir)/mpi_process
#mpi_mg_wr_dir           = $(mpi_operators_dir)/mpi_mg_wr
#
#bc_proc_direction_dir   = $(phy_model_bc_dir)/bc_procedure_direction
#bc_proc_order_dir   	= $(phy_model_bc_dir)/bc_procedure_order
#
#threevartest_dir	= $(phy_eq_dir)/threevartest
#fourvartest_dir		= $(phy_eq_dir)/fourvartest
#heat2d_dir              = $(phy_eq_dir)/heat2d
#
#profile_dir		= $(simulation_dir)/profile
#
#mpi_mg_bc_bl_dir	= $(mpi_mg_bc_dir)/mpi_mg_bc_blocking
#
#dim2d_param_dir		= $(dim2d_dir)/dim2d_parameters
#dim2d_state_dir		= $(dim2d_dir)/dim2d_state_eq
