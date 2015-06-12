import os

#constant variables
#------------------------------------------------------------
bc_size                           = 2    #number of grid points for the boundary layer
precision                         = 4    #only the four digits after the comma are kept (for dx,dy,dt,x_min,...)
nb_pts_inside_interface_required  = 10   #number of grid points to resolve the interface shape
md_threshold_ac_default           = 0    #the mass density criterion is desactivated by default for the openbc undermined
md_threshold_default              = 0.0  #the mass density criterion is desactivated by default for the openbc undermined
ic_perturbation_ac_default        = 0    #the addition of gaussian perturbation at the beginning is desactivated by default
ic_perturbation_amp_default       = 0.0  #the amplitude maximum of the perturbation is set to zero by default
li_perturbation_ac_default        = 0    #the addition of perturbation in the interface length at the beginning is desactivated by default
li_perturbation_amp_default       = 0.0  #the amplitude maximum of the perturbation is set to zero by default
bc_perturbation_T0_ac_default     = 0    #the addition of perturation to the far-field temperatureat is desactivated by default
bc_perturbation_T0_amp_default    = 0.0  #the amplitude maximum of the perturbation is set to zero by default
bc_perturbation_vx0_ac_default    = 0    #the addition of perturation to the far-field x-velocity is desactivated by default
bc_perturbation_vx0_amp_default   = 0.0  #the amplitude maximum of the perturbation is set to zero by default
bc_perturbation_vy0_ac_default    = 0    #the addition of perturation to the far-field y-velocity is desactivated by default
bc_perturbation_vy0_amp_default   = 0.0  #the amplitude maximum of the perturbation is set to zero by default
nb_pts_in_interface_default       = 10   #number of grid points required to resolve the interface gradient
ratio_bubble_interface_default    = 2    #ratio between the bubble and the interface at the initialization
CFL_constant_default              = 0.05 #CFL constant used by default
ratio_interface_influence_default = 4.0  #ratio used to determine at which distance the interface is supposed not to have an influence
dct_distance_default              = 4    #default distance between the edge and the detectors

total_nb_files_default            = 1000  #total number of files written
initial_reflection_study          = False #study of the initial refletions


changeParameterPath = os.path.join(os.getenv('augeanstables'),
                                   'src',
                                   'config',
                                   'change_parameter.sh')

getParameterPath    = os.path.join(os.getenv('augeanstables'),
                                   'src',
                                   'config',
                                   'get_parameter.sh')
