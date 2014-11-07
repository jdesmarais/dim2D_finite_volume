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


        !forcing term for the position
        real(rkind), parameter :: amplitude_force = 2.0d0
        real(rkind), parameter :: period_force    = period/c**2
        real(rkind), parameter :: x_center_force  = x_center
        real(rkind), parameter :: y_center_force  = y_center

      end module wave2d_parameters
