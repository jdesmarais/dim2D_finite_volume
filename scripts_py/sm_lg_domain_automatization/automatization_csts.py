#constant variables
#------------------------------------------------------------
bc_size                           = 2    #number of grid points for the boundary layer
precision                         = 4    #only the four digits after the comma are kept (for dx,dy,dt,x_min,...)
nb_pts_inside_interface_required  = 10   #number of grid points to resolve the interface shape
md_threshold_ac_default           = 0    #the mass density criterion is desactivated by default for the openbc undermined
md_threshold_default              = 0.0  #the mass density criterion is desactivated by default for the openbc undermined
nb_pts_in_interface_default       = 10   #number of grid points required to resolve the interface gradient
ratio_bubble_interface_default    = 2    #ratio between the bubble and the interface at the initialization
CFL_constant_default              = 0.1  #CFL constant used by default
ratio_interface_influence_default = 2.0  #ratio used to determine at which distance the interface is supposed not to have an influence
total_nb_files_default            = 1000 #total number of files written
