      module wave2d_parameters

        use parameters_kind, only : rkind


        !initial conditions parameters
        real(rkind), parameter :: amplitude = 1.0d0
        real(rkind), parameter :: period    = 2.0d0*ACOS(-1.0d0) !0.5d0
        real(rkind), parameter :: x_center  = 0.0d0
        real(rkind), parameter :: y_center  = 0.0d0


        !governing equations parameters
        real(rkind), parameter :: c   = 0.5d0
        real(rkind), parameter :: c_x = 0.0d0
        real(rkind), parameter :: c_y = 0.0d0
        real(rkind), parameter :: epsilon = 0.01d0


      end module wave2d_parameters
