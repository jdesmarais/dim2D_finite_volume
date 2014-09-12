      !> @file
      !> module implementing the subroutines computing the relaxation
      !> terms for the non-reflecting boundary conditions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> relaxation coefficients for non-reflecting yoo and lodato b.c.
      !
      !> @date
      !> 05_09_2014 - initial version   - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_relaxation_coeff_module

        use ns2d_parameters, only :
     $       gamma, mach_infty

        use parameters_constant, only :
     $       right

        use parameters_kind, only :
     $       rkind

        use parameters_input, only :
     $       sigma_P


        implicit none

        
        private
        public ::
     $       get_local_mach,
     $       get_relaxation_lodiT,
     $       get_relaxation_normal_velocity,
     $       get_relaxation_trans_velocity,
     $       get_relaxation_temperature,
     $       get_relaxation_pressure
        

        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the local Mach number \f$ M = \frac{\sqrt{u**2+v**2}}{c}\f$
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param u
        !> velocity component along the x-axis
        !
        !>@param v
        !> velocity component along the y-axis
        !
        !>@return M_local
        !> local Mach number
        !---------------------------------------------------------------
        function get_local_mach(u,v,c) result(M_local)
        
          implicit none

          real(rkind), intent(in) :: u
          real(rkind), intent(in) :: v
          real(rkind), intent(in) :: c
          real(rkind)             :: M_local

          M_local = SQRT(u**2+v**2)/c

        end function get_local_mach
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the relaxation coefficient for the LODI transverse
        !> components
        !> \f$ M \f$
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param M_local
        !> local Mach number
        !
        !>@return relaxationCoeff
        !> relaxation coefficient for the LODI transverse components
        !---------------------------------------------------------------
        function get_relaxation_lodiT(
     $       M_local)
     $       result(relaxationCoeff)

          implicit none

          real(rkind), intent(in) :: M_local
          real(rkind)             :: relaxationCoeff

          relaxationCoeff = M_local          

        end function get_relaxation_lodiT

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the relaxation coefficient for the velocity normal
        !> to the boundary
        !> \f$ \eta_{u_n} \frac{1}{l_n} \frac{M_{n,\infty}}{M^3_\infty}
        !> (1 - M^2_{n,\infty}) \f$
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param l_domain_n
        !> length of the domain in the direction normal to the boundary
        !
        !>@param M_un
        !> Mach number for the velocity normal to the boundary
        !
        !>@return relaxationCoeff
        !> relaxation coefficient for the velocity normal to the boundary
        !---------------------------------------------------------------
        function get_relaxation_normal_velocity(
     $       l_domain_n,
     $       M_un,
     $       side)
     $       result(relaxationCoeff)

          implicit none

          real(rkind), intent(in) :: l_domain_n
          real(rkind), intent(in) :: M_un
          logical    , intent(in) :: side
          real(rkind)             :: relaxationCoeff

          logical :: side_s


          if(rkind.eq.8) then
             relaxationCoeff =
     $            sigma_P*
     $            (1.0d0/l_domain_n)*
     $            (M_un/mach_infty**3)*
     $            (1.0d0-M_un**2)

          else
             relaxationCoeff =
     $            sigma_P*
     $            (1.0/l_domain_n)*
     $            (M_un/mach_infty**3)*
     $            (1.0-M_un**2)

          end if

          !if(side.eqv.right) then
          !   relaxationCoeff = -relaxationCoeff
          !end if
          side_s = side

        end function get_relaxation_normal_velocity


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the relaxation coefficient for the velocity transverse
        !> to the boundary
        !> \f$ \eta_{u_t} \frac{1}{l_n} \frac{1}{M_\infty}
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param l_domain_n
        !> length of the domain in the direction normal to the boundary
        !
        !>@param M_local
        !> local Mach number
        !
        !>@return relaxationCoeff
        !> relaxation coefficient for the velocity transverse to the
        !> boundary
        !---------------------------------------------------------------
        function get_relaxation_trans_velocity(
     $       l_domain_n,
     $       M_local)
     $       result(relaxationCoeff)

          implicit none

          real(rkind), intent(in) :: l_domain_n
          real(rkind), intent(in) :: M_local
          real(rkind)             :: relaxationCoeff
        
          
          relaxationCoeff = M_local/(l_domain_n*mach_infty)

        end function get_relaxation_trans_velocity


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the relaxation coefficient for the temperature
        !> \f$ \eta_{u_T} \frac{1}{l_n} \frac{1}{\gamma M_\infty^3}
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param l_domain_n
        !> length of the domain in the direction normal to the boundary
        !
        !>@param M_local
        !> local Mach number
        !
        !>@return relaxationCoeff
        !> relaxation coefficient for the temperature
        !---------------------------------------------------------------
        function get_relaxation_temperature(
     $       l_domain_n,
     $       M_local)
     $       result(relaxationCoeff)

          implicit none

          real(rkind), intent(in) :: l_domain_n
          real(rkind), intent(in) :: M_local
          real(rkind)             :: relaxationCoeff
        
          
          relaxationCoeff = M_local/(l_domain_n*gamma*mach_infty**3)

        end function get_relaxation_temperature


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the relaxation coefficient for the temperature
        !> \f$ \eta_{u_T} \frac{1}{l_n} \frac{1}{\gamma M_\infty^3}
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param l_domain_n
        !> length of the domain in the direction normal to the boundary
        !
        !>@param M_un
        !> Mach number for the velocity normal to the boundary
        !
        !>@return relaxationCoeff
        !> relaxation coefficient for the temperature
        !---------------------------------------------------------------
        function get_relaxation_pressure(
     $       l_domain_n,
     $       M_un)
     $       result(relaxationCoeff)

          implicit none

          real(rkind), intent(in) :: l_domain_n
          real(rkind), intent(in) :: M_un
          real(rkind)             :: relaxationCoeff
        
          
          if(rkind.eq.8) then
             relaxationCoeff = sigma_P/(l_domain_n*mach_infty)*
     $            (1.0d0-M_un**2)
          else
             relaxationCoeff = sigma_P/(l_domain_n*mach_infty)*
     $            (1.0-M_un**2)
          end if

c$$$          if(rkind.eq.8) then
c$$$             relaxationCoeff = sigma_P/(l_domain_n*mach_infty)*
c$$$     $            (1.0d0-mach_infty**2)
c$$$          else
c$$$             relaxationCoeff = sigma_P/(l_domain_n*mach_infty)*
c$$$     $            (1.0-mach_infty**2)
c$$$          end if

        end function get_relaxation_pressure

      end module lodi_relaxation_coeff_module
