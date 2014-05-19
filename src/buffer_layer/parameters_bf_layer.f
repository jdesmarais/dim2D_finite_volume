      module parameters_bf_layer

        implicit none

        integer, parameter :: no_pt = 0
        integer, parameter :: interior_pt = 1
        integer, parameter :: bc_interior_pt = 2
        integer, parameter :: bc_pt = 3
        !integer, parameter :: exchange_pt = 4

        logical, parameter :: clockwise=.false.
        logical, parameter :: counter_clockwise=.not.(clockwise)

      end module parameters_bf_layer
