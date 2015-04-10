      !> @file
      !> module with the constant used in the implementation of the
      !> buffer layer concepts
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module with the constant used in the implementation of the
      !> buffer layer concepts
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module parameters_bf_layer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,bc_size

        use parameters_kind, only :
     $       ikind

        implicit none


        !if a function returns successfully
        !------------------------------------------------------------
        logical, parameter :: BF_SUCCESS=.true.
        integer, parameter :: NEWGRDPT_NO_ERROR=0
        integer, parameter :: NEWGRDPT_PROC_NBGRDPTS_ERROR=1
        integer, parameter :: NEWGRDPT_DATA_NBGRDPTS_ERROR=2
        integer, parameter :: NEWGRDPT_NO_PREVIOUS_DATA_ERROR=3


        !identification of the grid points in the buffer layer
        !------------------------------------------------------------
        integer, parameter :: no_pt = 0
        integer, parameter :: interior_pt = 1
        integer, parameter :: bc_interior_pt = 2
        integer, parameter :: bc_pt = 3


        !indices identifying the beginning of the buffer layers
        !------------------------------------------------------------
        integer, parameter :: align_N = ny-bc_size+1
        integer, parameter :: align_S = bc_size
        integer, parameter :: align_E = nx-bc_size+1
        integer, parameter :: align_W = bc_size


        !convention for the determination of neighbor1 and neighbor2
        !------------------------------------------------------------
        !bf_neighbors(N,1) : cardinal coordinate corresponding to the
        !                    neighbor of type 1 for the main layer N
        !bf_neighbors(N,2) : cardinal coordinate corresponding to the
        !                    neighbor of type 2 for the main layer N
        !------------------------------------------------------------
        integer, dimension(4,2), parameter :: bf_neighbors =
     $       RESHAPE((/W,W,S,S,E,E,N,N/),(/4,2/))


        !convention for the determination of the type of neighbor
        !------------------------------------------------------------
        !bf_neighbors_id(N) : type of neighbor if the cardinal
        !                     coordinate is N
        !------------------------------------------------------------
        integer, dimension(4), parameter :: bf_neighbors_id = (/2,1,2,1/)


        !convention for the boundary layer procedures
        !------------------------------------------------------------
        integer, parameter :: no_bc_procedure_type=0
        integer, parameter :: SW_corner_type=1
        integer, parameter :: SE_corner_type=2
        integer, parameter :: NW_corner_type=3
        integer, parameter :: NE_corner_type=4
        integer, parameter :: S_edge_type=5
        integer, parameter :: E_edge_type=6
        integer, parameter :: W_edge_type=7
        integer, parameter :: N_edge_type=8
        integer, parameter :: SE_edge_type=9
        integer, parameter :: SW_edge_type=10
        integer, parameter :: NE_edge_type=11
        integer, parameter :: NW_edge_type=12


        !convention for the overlap of bc_sections
        !-------------------------------------------------------------
        integer, parameter :: no_overlap = 0
        integer, parameter :: N_overlap  = 1
        integer, parameter :: S_overlap  = 2
        integer, parameter :: E_overlap  = 3
        integer, parameter :: W_overlap  = 4
        integer, parameter :: NE_overlap = 5
        integer, parameter :: NW_overlap = 6
        integer, parameter :: SE_overlap = 7
        integer, parameter :: SW_overlap = 8
        integer, parameter :: NS_overlap = 9
        integer, parameter :: EW_overlap = 10


        !convention for the overlap of grid-points between corners
        !or between anti-corners
        !-------------------------------------------------------------
        integer, parameter :: cptnot_type     = 0
        integer, parameter :: cptnormal_type  = 1
        integer, parameter :: cptoverlap_type = 2

        integer, parameter :: cpt1normal_and_cpt4normal   = 0
        integer, parameter :: cpt1normal_and_cpt4not      = 1
        integer, parameter :: cpt1normal_and_cpt4overlap  = 2
        integer, parameter :: cpt1not_and_cpt4normal      = 3
        integer, parameter :: cpt1not_and_cpt4not         = 4
        integer, parameter :: cpt1not_and_cpt4overlap     = 5
        integer, parameter :: cpt1overlap_and_cpt4normal  = 6
        integer, parameter :: cpt1overlap_and_cpt4not     = 7
        integer, parameter :: cpt1overlap_and_cpt4overlap = 8

        integer, parameter :: cpt2normal_and_cpt3normal   =  0
        integer, parameter :: cpt2normal_and_cpt3not      = 11
        integer, parameter :: cpt2normal_and_cpt3overlap  = 12
        integer, parameter :: cpt2not_and_cpt3normal      = 13
        integer, parameter :: cpt2not_and_cpt3not         = 14
        integer, parameter :: cpt2not_and_cpt3overlap     = 15
        integer, parameter :: cpt2overlap_and_cpt3normal  = 16
        integer, parameter :: cpt2overlap_and_cpt3not     = 17
        integer, parameter :: cpt2overlap_and_cpt3overlap = 18


        !convention for the identification of the mainlayer interfaces
        !-------------------------------------------------------------
        integer, parameter :: NE_interface_type = 1
        integer, parameter :: NW_interface_type = 2
        integer, parameter :: SE_interface_type = 3
        integer, parameter :: SW_interface_type = 4


        !convention for new grid-point procedures
        !-------------------------------------------------------------
        integer, parameter :: no_gradient_type        = 0
        integer, parameter :: gradient_I_type         = 1
        integer, parameter :: gradient_L0_type        = 2
        integer, parameter :: gradient_R0_type        = 3
        integer, parameter :: gradient_xLR0_yI_type   = 4
        integer, parameter :: gradient_xI_yLR0_type   = 5
        integer, parameter :: gradient_xLR0_yLR0_type = 6


        !default positions for the increasing detectors
        !------------------------------------------------------------
        integer       , parameter :: dct_icr_distance  = bc_size
        integer(ikind), parameter :: dct_icr_N_default = ny-(bc_size+dct_icr_distance)+1
        integer(ikind), parameter :: dct_icr_S_default = bc_size+dct_icr_distance
        integer(ikind), parameter :: dct_icr_E_default = nx-(bc_size+dct_icr_distance)+1
        integer(ikind), parameter :: dct_icr_W_default = bc_size+dct_icr_distance


        !default strategy to update the position of the increasing
        !detectors
        ! - dct_velocity_strategy : the detectors are transported by
        !                           local velocity field
        ! - dct_bc_dir_strategy   : the detectors moves towards the
        !                           boundary
        !------------------------------------------------------------
        integer       , parameter :: dct_velocity_strategy = 0
        integer       , parameter :: dct_bc_dir_strategy   = 1
        integer       , parameter :: dct_update_strategy   = dct_bc_dir_strategy


        !default parameters for the decreasing detectors
        !------------------------------------------------------------
        !search_dcr  : radius expressed as number of grid
        !              points to check around the line
        !              for removing a buffer layer
        !------------------------------------------------------------
        integer, parameter :: search_dcr = 4


        !extra checks when determining the bc_procedures
        !------------------------------------------------------------
        !bc_procedure_extra_checks : when determining the boundary
        !                            procedures for the buffer layer,
        !                            assumptions can be made for
        !                            optimizations. However, for extra
        !                            security, i.e. to detect wrong updates
        !                            of the grdpts_id b/f they interfere
        !                            with the computations, extra checks
        !                            can be performed
        !------------------------------------------------------------
        integer, parameter :: bc_procedure_extra_checks = .true.

      end module parameters_bf_layer
