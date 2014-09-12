      module ns2d_parameters

        use parameters_kind, only : rkind

        implicit none


        !0.01 User defined variables
        !===========================
        !real(rkind), parameter :: length_c    = 1.0e-6 !length scale [m]

c$$$        real(rkind), parameter :: viscous_r   = -2.0d0/3.0d0
c$$$        real(rkind), parameter :: Re          = 10.0d0        
c$$$        real(rkind), parameter :: Pr          = 1.0d0
c$$$        real(rkind), parameter :: gamma       = 5.0d0/3.0d0
c$$$        real(rkind), parameter :: mach_infty  = 0.2d0 !1.0d0
c$$$        real(rkind), parameter :: gravity     = 0.03d0 !time_c/u_c*dim2d_g

        real(rkind), parameter :: viscous_r   = -2.0d0/3.0d0
        real(rkind), parameter :: Re          = 100000.d0
        real(rkind), parameter :: Pr          = 100.d0
        real(rkind), parameter :: gamma       = 5.0d0/3.0d0
        real(rkind), parameter :: mach_infty  = 0.1d0
        real(rkind), parameter :: gravity     = 0.03d0 !time_c/u_c*dim2d_g

        real(rkind), parameter :: epsilon     = 1.0d0/Re

      end module ns2d_parameters
