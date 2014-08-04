      module wave1d_parameters

        use parameters_kind, only : rkind


        !initial conditions parameters
        real(rkind), parameter :: amplitude = 1.0d0
        real(rkind), parameter :: period    = 0.5d0
        real(rkind), parameter :: x_center  = 0.0d0
        real(rkind), parameter :: y_center  = 0.0d0


        !governing equations parameters
        real(rkind), parameter :: c   = 1.0d0


      end module wave1d_parameters
