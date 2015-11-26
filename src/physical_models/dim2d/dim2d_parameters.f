      !> @file
      !> Definition of the constants appearing in the Diffuse Interface
      !> Model governing equations: viscosity ratio, reynolds number...
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> Definition of the constants appearing in the Diffuse Interface
      !> Model governing equations: viscosity ratio, Reynolds number,
      !> Prandtl number, Weber number, reduced heat capacity ...
      !
      !> @date
      !> 09_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_parameters

        use parameters_input, only :
     $       gravity_amp

        use parameters_kind, only :
     $       rkind


        private
        public ::
     $       rho_c, T_c, u_c, length_c, time_c,
     $       viscous_r, re, pr, we, cv_r, gravity,
     $       epsilon, zeta


        !============================================================
        ! input quantities for water
        !============================================================

        !>@brief
        !> Van der Waals constant modelling the pressure corrections,
        !> unit=[J.m3.mol-2]
        !> \f[ a \f]
        !------------------------------------------------------------
        real(rkind), parameter :: dim2d_a      = 0.5536

        !>@brief
        !> Van der Waals constant modelling the covolume, unit=[m3.mol-1]
        !> \f[ b \f]
        !------------------------------------------------------------
        real(rkind), parameter :: dim2d_b      = 0.03049e-3

        !>@brief
        !> molar mass of the fluid, unit=[kg.mol-1]
        !> \f[ M \f]
        !---------------------------------------------------------        
        real(rkind), parameter :: dim2d_M      = 18.0e-3        
        
        !>@brief
        !> universal gas constant, unit=[J.mol-1.K-1]
        !> \f[ R \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_R      = 8.314      

        !>@brief
        !> viscosity coefficient, unit=[Pa.s]
        !> \f[ \mu \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_mu     = 55.25e-6   

        !>@brief
        !> viscosity coefficient, unit=[Pa.s]
        !> \f[ \eta \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_nu     = -36.83e-6  

        !>@brief
        !> capillarity coefficient, unit=[J.kg-2.m5]
        !> \f[ K \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_K      = 6.81e-15   

        !>@brief
        !> heat capacity at constant volume, unit=[J.K-1.kg-1]
        !> \f[ c_v \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_cv     = 1410.      

        !>@brief
        !> thermal conductivity, unit=[J.s-1.K-1.m-1]
        !> \f[ \lambda \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_lambda = 0.6

        !>@brief
        !> universal acceleration, unit=[m.s-2]
        !> \f[ g \f]
        !---------------------------------------------------------         
        real(rkind), parameter :: dim2d_g      = 9.81         


        !============================================================
        ! intermediate variables
        !============================================================

        !>@brief
        !> critical density, unit=[kg.m3]
        !> \f[ \rho_c = \displaystyle{\frac{M}{3 b}}\f]
        !------------------------------------------------------------
        real(rkind), parameter :: rho_c    = dim2d_M/(3.0d0*dim2d_b)

        !@brief
        !> critical pressure, unit=[Pa]
        !> \f[ P_c = \displaystyle{\frac{a}{27 b^2}}\f]
        !------------------------------------------------------------
        real(rkind), parameter :: p_c      = dim2d_a/(27.0d0*dim2d_b**2)

        !@brief
        !> critical temperature, unit=[K]
        !> \f[ T_c = \displaystyle{\frac{8 a}{27 R b}}\f]
        !------------------------------------------------------------
        real(rkind), parameter :: T_c      = 8.0d0*dim2d_a/(27.0d0*dim2d_b*dim2d_R)

        !@brief
        !> critical velocity, unit=[m.s-1]
        !> \f[ u_c = \displaystyle{\sqrt{\frac{P_c}{\rho_c}}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: u_c      = (p_c/rho_c)**0.5

        !@brief
        !> length scale, unit=[m]
        !> \f[ L_c \f]
        !------------------------------------------------------------
        real(rkind), parameter :: length_c = 1.0e-6

        !@brief
        !> time scale, unit=[s]
        !> \f[ t_c = \displaystyle{\frac{L_c}{u_c}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: time_c   = length_c/u_c


        !============================================================
        ! reduced parameters used in the governing equations
        !============================================================
        !>@brief
        !> viscous ratio, unit=[-]
        !> \f[ \displaystyle{\frac{\eta}{\mu}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: viscous_r = dim2d_nu/dim2d_mu

        !>@brief
        !> Reynolds number, unit=[-]    
        !> \f[ Re = \displaystyle{\frac{\rho_c u_c L_c}{\mu}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: Re = rho_c*u_c*length_c/dim2d_mu

        !>@brief
        !> Weber number, unit=[-]
        !> \f[ We = \displaystyle{\frac{L_c^2 u_c^2}{\rho_c K}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: We = (length_c**2*u_c**2)/(rho_c*dim2d_K)

        !>@brief
        !> Prandtl number, unit=[-]    
        !> \f[ Pr = \displaystyle{\frac{\mu c_v}{\lambda}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: Pr = dim2d_mu*dim2d_cv/dim2d_lambda

        !>@brief
        !> heat capacity reduced, unit=[-]    
        !> \f[ \tilde{c_v} = \displaystyle{\frac{M c_v}{R}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: cv_r = dim2d_M*dim2d_cv/dim2d_R

        !>@brief
        !> universal acceleration reduced, unit=[-]    
        !> \f[ \tilde{g} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: gravity = gravity_amp !time_c/u_c*dim2d_g

        !--------------------------------------------------------
        !cv=dim2d_M*dim2d_cv/dim2d_R
        ! for tests
        !--------------------------------------------------------
        !real(rkind), parameter :: viscous_r = -1.5d0
        !real(rkind), parameter :: Re = 5.0d0
        !real(rkind), parameter :: We = 10.0d0
        !real(rkind), parameter :: Pr = 20.0d0
        !real(rkind), parameter :: cv_r = 2.5d0
        !real(rkind), parameter :: gravity = 9.81d0

        !>@brief
        !> inverse of the Reynolds number
        !> \f[ \epsilon = \displaystyle{\frac{1}{Re}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: epsilon = 1.0d0/Re

        !>@brief
        !> inverse of the Weber number
        !> \f[ \zeta = \displaystyle{\frac{1}{We}} \f]
        !------------------------------------------------------------
        real(rkind), parameter :: zeta    = 1.0d0/We


      end module dim2d_parameters



