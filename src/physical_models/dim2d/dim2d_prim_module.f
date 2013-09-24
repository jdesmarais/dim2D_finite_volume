      !> @file
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the Diffuse Interface
      !> Model governing equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the Diffuse Interface
      !> Model governing equations
      !
      !> @date
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_prim_module

        use dim2d_parameters, only : cv_r, we
        use field_class     , only : field
        use parameters_kind , only : ikind, rkind

        implicit none

        private
        public :: mass_density, momentum_x, momentum_y, total_energy,
     $            velocity_x, velocity_y,
     $            classical_pressure, temperature_eff,
     $            classical_pressure_xwork, classical_pressure_ywork,
     $            qx_transport_x, qy_transport_x,
     $            qx_transport_y, qy_transport_y,
     $            energy_transport_x, energy_transport_y,
     $            capillarity_pressure,
     $            capillarity_pressure_xwork, capillarity_pressure_ywork

        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the mass density \f$ \rho \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function mass_density(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,1)

        end function mass_density


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum density along the x-axis
        !> \f$ \rho u_x \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_x \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function momentum_x(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,2)

        end function momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum density along the y-axis
        !> \f$ \rho u_y \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_y \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function momentum_y(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,3)

        end function momentum_y

      
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the total energy density \f$\rho E\f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$\rho E\f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function total_energy(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,4)

        end function total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the velocity along the x-axis
        !> \f$ u_x \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ u_x \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_x(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,2)/field_used%nodes(i,j,1)

        end function velocity_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the velocity along the y-axis
        !> \f$ u_y \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ u_y \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_y(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,3)/field_used%nodes(i,j,1)

        end function velocity_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical pressure
        !> \f$ \frac{3}{(3-\rho) c_v}
        !> \left( \rho E - \frac{1}{2} \rho (u_x^2 + u_y^2) + 3\rho^2
        !> \right) - 3 \rho^2 \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ P \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function classical_pressure(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then
             var=3.0d0/((3.0d0-field_used%nodes(i,j,1))*cv_r)*(
     $            field_used%nodes(i,j,4)
     $            - 0.5d0*field_used%nodes(i,j,1)*(
     $            (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $            (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $            + 3.0d0*field_used%nodes(i,j,1)**2)
     $            - 3.0d0*field_used%nodes(i,j,1)**2
          else
             var=3./((3.-field_used%nodes(i,j,1))*cv_r)*(
     $            field_used%nodes(i,j,4)
     $            - 1./2.*field_used%nodes(i,j,1)*(
     $            (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $            (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $            + 3*field_used%nodes(i,j,1)**2)
     $            - 3*field_used%nodes(i,j,1)**2
          end if

        end function classical_pressure


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the effective temperature
        !> \f$ \frac{1}{\rho}
        !> \left( \rho E - \frac{1}{2} \rho (u_x^2 + u_y^2)
        !> - \frac{1}{2 We} {|\nabla \rho|}^2 + 3 \rho^2 \right)
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ T_{\textrm{eff}} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function temperature_eff(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then

          var=1.0d0/(field_used%nodes(i,j,1))*(
     $           field_used%nodes(i,j,4)
     $           - 0.5d0*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           - 0.5d0/we*((
     $              (field_used%nodes(i+1,j,1)-field_used%nodes(i-1,j,1))
     $              /(2.0d0*field_used%dx))**2+(
     $              (field_used%nodes(i,j+1,1)-field_used%nodes(i,j-1,1))
     $              /(2.0d0*field_used%dy))**2)
     $           + 3.0d0*field_used%nodes(i,j,1)**2)

          else
             var=1./(field_used%nodes(i,j,1))*(
     $           field_used%nodes(i,j,4)
     $           - 1./2.*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           - 1./(2*we)*((
     $              (field_used%nodes(i+1,j,1)-field_used%nodes(i-1,j,1))
     $              /(2*field_used%dx))**2+(
     $              (field_used%nodes(i,j+1,1)-field_used%nodes(i,j-1,1))
     $              /(2*field_used%dy))**2)
     $           + 3*field_used%nodes(i,j,1)**2)
          end if

        end function temperature_eff


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical pressure work along the x-axis
        !> \f$ \left(\frac{3}{(3-\rho) c_v}
        !> \left[ \rho E - \frac{1}{2} \rho (u_x^2 + u_y^2) + 3\rho^2
        !> \right] - 3 \rho^2 \right) u_x \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> work of \f$ P \f$ along the x-axis evaluated at [i,j]
        !---------------------------------------------------------------
        function classical_pressure_xwork(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then

             var=(3.0d0/((3.0d0-field_used%nodes(i,j,1))*cv_r)*(
     $           field_used%nodes(i,j,4)
     $           - 0.5d0*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           + 3.0d0*field_used%nodes(i,j,1)**2)
     $         - 3.0d0*field_used%nodes(i,j,1)**2)*
     $         field_used%nodes(i,j,2)/field_used%nodes(i,j,1)

          else
             var=(3./((3.-field_used%nodes(i,j,1))*cv_r)*(
     $           field_used%nodes(i,j,4)
     $           - 1./2.*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           + 3*field_used%nodes(i,j,1)**2)
     $         - 3*field_used%nodes(i,j,1)**2)*
     $         field_used%nodes(i,j,2)/field_used%nodes(i,j,1)

          end if

        end function classical_pressure_xwork


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical pressure work along the y-axis
        !> \f$ \left(\frac{3}{(3-\rho) c_v}
        !> \left[ \rho E - \frac{1}{2} \rho (u_x^2 + u_y^2) + 3\rho^2
        !> \right] - 3 \rho^2 \right) u_y \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> work of \f$ P \f$ along the y-axis evaluated at [i,j]
        !---------------------------------------------------------------
        function classical_pressure_ywork(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then
             var=(3.0d0/((3.0d0-field_used%nodes(i,j,1))*cv_r)*(
     $           field_used%nodes(i,j,4)
     $           - 0.5d0*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           + 3.0d0*field_used%nodes(i,j,1)**2)
     $         - 3.0d0*field_used%nodes(i,j,1)**2)*
     $         field_used%nodes(i,j,3)/field_used%nodes(i,j,1)

          else
             var=(3./((3.-field_used%nodes(i,j,1))*cv_r)*(
     $           field_used%nodes(i,j,4)
     $           - 1./2.*field_used%nodes(i,j,1)*(
     $              (field_used%nodes(i,j,2)/field_used%nodes(i,j,1))**2+
     $              (field_used%nodes(i,j,3)/field_used%nodes(i,j,1))**2)
     $           + 3*field_used%nodes(i,j,1)**2)
     $         - 3*field_used%nodes(i,j,1)**2)*
     $         field_used%nodes(i,j,3)/field_used%nodes(i,j,1)
          end if

        end function classical_pressure_ywork


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_x transported along the x-axis
        !> \f$ \rho u_x^2 \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_x^2 \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qx_transport_x(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,2)**2/field_used%nodes(i,j,1)

        end function qx_transport_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_y transported along the x-axis
        !> \f$ \rho u_y u_x \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_y u_x \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qy_transport_x(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,3)*
     $         field_used%nodes(i,j,2)/
     $         field_used%nodes(i,j,1)

        end function qy_transport_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_x transported along the y-axis
        !> \f$ \rho u_x u_y \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_x u_y \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qx_transport_y(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,2)*
     $         field_used%nodes(i,j,3)/
     $         field_used%nodes(i,j,1)

        end function qx_transport_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_y transported along the y-axis
        !> \f$ \rho u_y^2 \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho u_y^2 \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qy_transport_y(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,3)**2/field_used%nodes(i,j,1)

        end function qy_transport_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the energy transported along the x-axis
        !> \f$ \rho E u_x \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho E u_x \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function energy_transport_x(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,4)*
     $         field_used%nodes(i,j,2)/
     $         field_used%nodes(i,j,1)

        end function energy_transport_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the energy transported along the y-axis
        !> \f$ \rho E u_y \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho E u_y \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function energy_transport_y(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,4)*
     $         field_used%nodes(i,j,3)/
     $         field_used%nodes(i,j,1)

        end function energy_transport_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the first term in the capillarity pressure
        !> \f$ \frac{1}{3-\rho} \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \frac{1}{3-\rho} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function capillarity_pressure(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then
             var=1.0d0/(3.0d0-field_used%nodes(i,j,1))
          else
             var=1./(3.-field_used%nodes(i,j,1))
          end if

        end function capillarity_pressure


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the first term in the capillarity pressure work
        !> along the x-axis \f$ \frac{u_x}{3-\rho} \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \frac{u_x}{3-\rho} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function capillarity_pressure_xwork(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then
             var=field_used%nodes(i,j,2)/
     $          (field_used%nodes(i,j,1)*(3.0d0-field_used%nodes(i,j,1)))
          else
             var=field_used%nodes(i,j,2)/
     $          (field_used%nodes(i,j,1)*(3.-field_used%nodes(i,j,1)))
          end if

        end function capillarity_pressure_xwork


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the first term in the capillarity pressure work
        !> along the y-axis \f$ \frac{u_y}{3-\rho} \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \frac{u_y}{3-\rho} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function capillarity_pressure_ywork(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          if(rkind.eq.8) then
             var=field_used%nodes(i,j,3)/
     $          (field_used%nodes(i,j,1)*(3.0d0-field_used%nodes(i,j,1)))
          else
             var=field_used%nodes(i,j,3)/
     $          (field_used%nodes(i,j,1)*(3.-field_used%nodes(i,j,1)))
          end if

        end function capillarity_pressure_ywork

      end module dim2d_prim_module
