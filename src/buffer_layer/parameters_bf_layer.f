      module parameters_bf_layer

        use parameters_input, only : nx, ny, bc_size

        implicit none

        integer, parameter :: no_pt = 0
        integer, parameter :: interior_pt = 1
        integer, parameter :: bc_interior_pt = 2
        integer, parameter :: bc_pt = 3
        !integer, parameter :: exchange_pt = 4

        !indices identifying the beginning of the buffer layers
        integer, parameter :: align_N = ny-bc_size+1
        integer, parameter :: align_S = bc_size
        integer, parameter :: align_E = nx-bc_size+1
        integer, parameter :: align_W = bc_size

        logical, parameter :: clockwise=.false.
        logical, parameter :: counter_clockwise=.not.(clockwise)

      end module parameters_bf_layer
