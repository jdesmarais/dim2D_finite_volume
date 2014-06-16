      module parameters_bf_layer

        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : nx, ny, bc_size

        implicit none

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



        logical, parameter :: clockwise=.false.
        logical, parameter :: counter_clockwise=.not.(clockwise)

      end module parameters_bf_layer
