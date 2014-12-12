      !> @file
      !> definition of the constants appearing in the Diffuse Interface
      !> Model governing equations: viscosity ratio, reynolds number...
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> definition of the constants appearing in the Diffuse Interface
      !> Model governing equations: viscosity ratio, Reynolds number,
      !> Prandtl number, Weber number, reduced heat capacity ...
      !
      !> @date
      ! 09_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_parameters

      
        use parameters_kind, only : rkind


        private
        public ::
     $       rho_c, T_c, u_c, length_c, time_c,
     $       viscous_r, re, pr, we, cv_r, gravity,
     $       epsilon, zeta


        ! input quantities for water
        !---------------------------------------------------------
        !@param dim2d_a       Van der Waals constant   [J.m3.mol-2]     
        !@param dim2d_b       Van der Waals constant   [m3.mol-1]       
        !@param dim2d_M       molar mass of the fluid  [kg.mol-1]       
        !@param dim2d_R       universal gas constant   [J.mol-1.K-1]    
        !@param dim2d_mu      viscosity coefficient    [Pa.s]
        !@param dim2d_nu      viscosity coefficient    [Pa.s] 
        !@param dim2d_K       capillarity coefficient  [J.kg-2.m5]      
        !@param dim2d_cv      heat capacity at vol cst [J.K-1.kg-1]   
        !@param dim2d_lambda  thermal conductivity     [J.s-1.K-1.m-1]
        !@param dim2d_g       universal acceleration   [m.s-2]
        !---------------------------------------------------------
        real(rkind), parameter :: dim2d_a      = 0.5536     
        real(rkind), parameter :: dim2d_b      = 0.03049e-3 
        real(rkind), parameter :: dim2d_M      = 18.0e-3   
        real(rkind), parameter :: dim2d_R      = 8.314      
        real(rkind), parameter :: dim2d_mu     = 55.25e-6   
        real(rkind), parameter :: dim2d_nu     = -36.83e-6  
        real(rkind), parameter :: dim2d_K      = 6.81e-15   
        real(rkind), parameter :: dim2d_cv     = 1410.      
        real(rkind), parameter :: dim2d_lambda = 0.6
        real(rkind), parameter :: dim2d_g      = 9.81


        ! intermediate variables
        !--------------------------------------------------------
        !@param rho_c          critical density           [kg.m3]
        !@param p_c            critical pressure          [Pa]   
        !@param T_c            critical temperature       [K]    
        !@param u_c            critical velocity          [m.s-1]
        !@param length_c       length scale               [m]    
        !@param time_c         time scale                 [s]
        !--------------------------------------------------------
        real(rkind), parameter :: rho_c   = dim2d_M/(3*dim2d_b)
        real(rkind), parameter :: p_c     = dim2d_a/(27*dim2d_b**2)
        real(rkind), parameter :: T_c     = 8*dim2d_a/(27*dim2d_b*dim2d_R)
        real(rkind), parameter :: u_c      = (p_c/rho_c)**0.5
        real(rkind), parameter :: length_c = 1.0e-6
        real(rkind), parameter :: time_c   = length_c/u_c


        ! parameters initialized
        !--------------------------------------------------------
        !@param viscous_r      viscous ratio                  [-]
        !@param Re             Reynolds number                [-]    
        !@param We             Weber number                   [-]    
        !@param Pr             Prandtl number                 [-]    
        !@param cv_r           heat capacity reduced          [-]
        !@param gravity        universal acceleration reduced [-]
        !--------------------------------------------------------
        !real(rkind), parameter :: viscous_r= dim2d_nu/dim2d_mu
        !real(rkind), parameter :: Re = rho_c*u_c*length_c/dim2d_mu
        !real(rkind), parameter :: We =(length_c**2*u_c**2)/(rho_c*dim2d_K)
        !real(rkind), parameter :: Pr = dim2d_mu*dim2d_cv/dim2d_lambda
        !real(rkind), parameter :: cv_r = dim2d_M*dim2d_cv/dim2d_R
        !real(rkind), parameter :: gravity = 0.03d0 !time_c/u_c*dim2d_g
       
        real(rkind), parameter :: viscous_r= -1.5
        real(rkind), parameter :: Re = 5
        real(rkind), parameter :: We = 10
        real(rkind), parameter :: Pr = 20
        real(rkind), parameter :: cv_r = 2.5
        real(rkind), parameter :: gravity = 9.81

        real(rkind), parameter :: epsilon = 1.0d0/Re
        real(rkind), parameter :: zeta    = 1.0d0/We

      end module dim2d_parameters



