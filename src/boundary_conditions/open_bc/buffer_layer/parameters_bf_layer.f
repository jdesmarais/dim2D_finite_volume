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

        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : nx, ny, bc_size
        use parameters_kind    , only : ikind

        implicit none

        !if a function returns successfully
        logical, parameter :: BF_SUCCESS=.true.

        !identification of the grid points in the buffer layer
        integer, parameter :: no_pt = 0
        integer, parameter :: interior_pt = 1
        integer, parameter :: bc_interior_pt = 2
        integer, parameter :: bc_pt = 3

        !indices identifying the beginning of the buffer layers
        integer, parameter :: align_N = ny-bc_size+1
        integer, parameter :: align_S = bc_size
        integer, parameter :: align_E = nx-bc_size+1
        integer, parameter :: align_W = bc_size


        !convention for the determination of neighbor1 and neighbor2
        !bf_neighbors(N,1) : cardinal coordinate corresponding to the
        !                    neighbor of type 1 for the main layer N
        !bf_neighbors(N,2) : cardinal coordinate corresponding to the
        !                    neighbor of type 2 for the main layer N        
        integer, dimension(4,2), parameter :: bf_neighbors =
     $       RESHAPE((/W,W,S,S,E,E,N,N/),(/4,2/))

        !convention for the determination of the type of neighbor
        !bf_neighbors_id(N) : type of neighbor if the cardinal
        !                     coordinate is N
        integer, dimension(4), parameter :: bf_neighbors_id = (/2,1,2,1/)


        !default positions for the increasing detectors
        integer       , parameter :: dct_icr_distance  = 2*bc_size
        integer(ikind), parameter :: dct_icr_N_default = ny-(bc_size+dct_icr_distance)+1
        integer(ikind), parameter :: dct_icr_S_default = bc_size+dct_icr_distance
        integer(ikind), parameter :: dct_icr_E_default = nx-(bc_size+dct_icr_distance)+1
        integer(ikind), parameter :: dct_icr_W_default = bc_size+dct_icr_distance

      end module parameters_bf_layer
