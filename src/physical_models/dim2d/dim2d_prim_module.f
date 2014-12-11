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

        use dim2d_parameters, only :
     $       cv_r,
     $       we

        use interface_primary, only :
     $       get_primary_var,
     $       gradient_x_proc,
     $       gradient_y_proc

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public ::
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       velocity_x,
     $       velocity_y,
     $       classical_pressure,
     $       classical_temperature_eff,
     $       capillarity_temperature_eff,
     $       temperature_eff,
     $       classical_pressure_xwork,
     $       classical_pressure_ywork,
     $       qx_transport_x,
     $       qy_transport_x,
     $       qx_transport_y,
     $       qy_transport_y,
     $       energy_transport_x,
     $       energy_transport_y,
     $       capillarity_pressure,
     $       capillarity_pressure_xwork,
     $       capillarity_pressure_ywork,
     $       speed_of_sound,
     $       compute_jacobian_prim_to_cons,
     $       compute_jacobian_cons_to_prim


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
        !>@param nodes
        !> array with the grid point data
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
        function mass_density(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function momentum_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)

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
        !>@param nodes
        !> array with the grid point data
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
        function momentum_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)

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
        !>@param nodes
        !> array with the grid point data
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
        function total_energy(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,4)

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
        !>@param nodes
        !> array with the grid point data
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
        function velocity_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)/nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function velocity_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)/nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function classical_pressure(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var=3.0d0/((3.0d0-nodes(i,j,1))*cv_r)*(
     $            nodes(i,j,4)
     $            - 0.5d0*nodes(i,j,1)*(
     $            (nodes(i,j,2)/nodes(i,j,1))**2+
     $            (nodes(i,j,3)/nodes(i,j,1))**2)
     $            + 3.0d0*nodes(i,j,1)**2)
     $            - 3.0d0*nodes(i,j,1)**2
          else
             var=3./((3.-nodes(i,j,1))*cv_r)*(
     $            nodes(i,j,4)
     $            - 1./2.*nodes(i,j,1)*(
     $            (nodes(i,j,2)/nodes(i,j,1))**2+
     $            (nodes(i,j,3)/nodes(i,j,1))**2)
     $            + 3*nodes(i,j,1)**2)
     $            - 3*nodes(i,j,1)**2
          end if

        end function classical_pressure


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the effective temperature
        !> \f$ \frac{1}{\rho}
        !> \left( \rho E - \frac{1}{2} \rho (u_x^2 + u_y^2)
        !>  + 3 \rho^2 \right)
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
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
        function classical_temperature_eff(
     $     nodes,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then

          var=1.0d0/(nodes(i,j,1))*(
     $           nodes(i,j,4)
     $           - 0.5d0*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3.0d0*nodes(i,j,1)**2)

          else
             var=1./(nodes(i,j,1))*(
     $           nodes(i,j,4)
     $           - 0.5*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3*nodes(i,j,1)**2)
          end if

        end function classical_temperature_eff


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the effective temperature
        !> \f$ \frac{1}{2 \rho} {|\nabla \rho|}^2 \f$
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
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
        function capillarity_temperature_eff(
     $     nodes,i,j,dx,dy,gradient_x,gradient_y)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind)                               :: var

          if(rkind.eq.8) then

             var=0.5d0/(nodes(i,j,1))*(
     $            (gradient_x(nodes,i,j,mass_density,dx))**2 +
     $            (gradient_y(nodes,i,j,mass_density,dy))**2
     $            )

          else

             var=0.5/(nodes(i,j,1))*(
     $            (gradient_x(nodes,i,j,mass_density,dx))**2 +
     $            (gradient_y(nodes,i,j,mass_density,dy))**2
     $            )

          end if

        end function capillarity_temperature_eff


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
        !>@param nodes
        !> array with the grid point data
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
        function temperature_eff(
     $     nodes,i,j,dx,dy,gradient_x,gradient_y)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind)                               :: var

          if(rkind.eq.8) then

          var=1.0d0/(nodes(i,j,1))*(
     $           nodes(i,j,4)
     $           - 0.5d0*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           - 0.5d0/we*(
     $              (gradient_x(nodes,i,j,mass_density,dx))**2
     $            + (gradient_y(nodes,i,j,mass_density,dy))**2)
     $           + 3.0d0*nodes(i,j,1)**2)

          else
             var=1./(nodes(i,j,1))*(
     $           nodes(i,j,4)
     $           - 0.5*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           - 0.5/we*(
     $              (gradient_x(nodes,i,j,mass_density,dx))**2
     $            + (gradient_y(nodes,i,j,mass_density,dy))**2)
     $           + 3*nodes(i,j,1)**2)
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
        !>@param nodes
        !> array with the grid point data
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
        function classical_pressure_xwork(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then

             var=(3.0d0/((3.0d0-nodes(i,j,1))*cv_r)*(
     $           nodes(i,j,4)
     $           - 0.5d0*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3.0d0*nodes(i,j,1)**2)
     $         - 3.0d0*nodes(i,j,1)**2)*
     $         nodes(i,j,2)/nodes(i,j,1)

          else
             var=(3./((3.-nodes(i,j,1))*cv_r)*(
     $           nodes(i,j,4)
     $           - 1./2.*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3*nodes(i,j,1)**2)
     $         - 3*nodes(i,j,1)**2)*
     $         nodes(i,j,2)/nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function classical_pressure_ywork(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var=(3.0d0/((3.0d0-nodes(i,j,1))*cv_r)*(
     $           nodes(i,j,4)
     $           - 0.5d0*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3.0d0*nodes(i,j,1)**2)
     $         - 3.0d0*nodes(i,j,1)**2)*
     $         nodes(i,j,3)/nodes(i,j,1)

          else
             var=(3./((3.-nodes(i,j,1))*cv_r)*(
     $           nodes(i,j,4)
     $           - 1./2.*nodes(i,j,1)*(
     $              (nodes(i,j,2)/nodes(i,j,1))**2+
     $              (nodes(i,j,3)/nodes(i,j,1))**2)
     $           + 3*nodes(i,j,1)**2)
     $         - 3*nodes(i,j,1)**2)*
     $         nodes(i,j,3)/nodes(i,j,1)
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
        !>@param nodes
        !> array with the grid point data
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
        function qx_transport_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)**2/nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function qy_transport_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)*
     $         nodes(i,j,2)/
     $         nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function qx_transport_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)*
     $         nodes(i,j,3)/
     $         nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function qy_transport_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)**2/nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function energy_transport_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,4)*
     $         nodes(i,j,2)/
     $         nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function energy_transport_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,4)*
     $         nodes(i,j,3)/
     $         nodes(i,j,1)

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
        !>@param nodes
        !> array with the grid point data
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
        function capillarity_pressure(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var=1.0d0/(3.0d0-nodes(i,j,1))
          else
             var=1./(3.-nodes(i,j,1))
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
        !>@param nodes
        !> array with the grid point data
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
        function capillarity_pressure_xwork(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var=nodes(i,j,2)/
     $          (nodes(i,j,1)*(3.0d0-nodes(i,j,1)))
          else
             var=nodes(i,j,2)/
     $          (nodes(i,j,1)*(3.-nodes(i,j,1)))
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
        !>@param nodes
        !> array with the grid point data
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
        function capillarity_pressure_ywork(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var=nodes(i,j,3)/
     $          (nodes(i,j,1)*(3.0d0-nodes(i,j,1)))
          else
             var=nodes(i,j,3)/
     $          (nodes(i,j,1)*(3.-nodes(i,j,1)))
          end if

        end function capillarity_pressure_ywork


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the speed of sound given the data at the grid
        !> point location
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !>
        !>@param var
        !> \f$ c = \frac{\sqrt{2} \rho}{a d} \f$ where
        !> \f$ a = \frac{1}{\sqrt{1+b^2}}\f$
        !> \f$ b = \sqrt{\frac{P + 3 \rho^2}{c_V ( P+ \rho^2 (-3+2 \rho))}}\f$
        !> \f$ d = \sqrt{\rho^3} \sqrt{\frac{3-\rho}{P+\rho^2(-3+2\rho)}} \f$
        !---------------------------------------------------------------
        function speed_of_sound(nodes) result(var)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind)                            :: var


          real(rkind) :: P
          real(rkind) :: a
          real(rkind) :: b
          real(rkind) :: d
          real(rkind) :: e


          if(rkind.eq.8) then

             P = 3.0d0/((3.0d0-nodes(1))*cv_r)*(
     $            nodes(4)
     $            - 0.5d0/nodes(1)*(nodes(2)**2+nodes(3)**2)
     $            + 3.0d0*nodes(1)**2)
     $            - 3.0d0*nodes(1)**2

             e = P+nodes(1)**2*(-3.0d0+2.0d0*nodes(1))
             b = SQRT((P+3.0d0*nodes(1)**2)/(cv_r*e))
             a = 1.0d0/SQRT(1.0d0+b**2)
             d = SQRT(nodes(1)**3)*SQRT((3.0d0-nodes(1))/e)

             var = SQRT(3.0d0)*nodes(1)/(a*d)

          else

             P = 3.0/((3.0-nodes(1))*cv_r)*(
     $            nodes(4)
     $            - 0.5/nodes(1)*(nodes(2)**2+nodes(3)**2)
     $            + 3.0*nodes(1)**2)
     $            - 3.0*nodes(1)**2
             e = P+nodes(1)**2*(-3.0+2.0*nodes(1))
             b = SQRT((P+3.0*nodes(1)**2)/(cv_r*e))
             a = 1.0/SQRT(1.0+b**2)
             d = SQRT(nodes(1)**3)*SQRT((3.0-nodes(1))/e)

             var = SQRT(3.0d0)*nodes(1)/(a*d)

          end if

        end function speed_of_sound


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for primitive to
        !> to conservative variables
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacPrimCons
        !> jacobian matrix for primitive to conservative
        !> variables \f$ \frac{\partial p}{\partial v} \f$
        !--------------------------------------------------------------
        function compute_jacobian_prim_to_cons(nodes)
     $     result(jacPrimCons)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacPrimCons

          real(rkind) :: ux
          real(rkind) :: uy
          real(rkind) :: ratio
          real(rkind) :: dPdrho

          ux    = nodes(2)/nodes(1)
          uy    = nodes(3)/nodes(1)

          if(rkind.eq.8) then             

             ratio = 3.0d0/(cv_r*(-3.0d0+nodes(1)))

             dPdrho = -3.0d0*(nodes(2)**2*(-3.0d0 + 2.0d0*nodes(1)) + 
     $            nodes(3)**2*(-3.0d0 + 2.0d0*nodes(1)) + 
     $            2.0d0*nodes(1)**2*((3.0d0*(-6.0d0 + nodes(1)) + 
     $            2.0d0*cv_r*(-3.0d0 + nodes(1))**2)*nodes(1) - nodes(4)))/
     $            (2.0d0*cv_r*(-3.0d0 + nodes(1))**2*nodes(1)**2)

             jacPrimCons(1,1) = 1.0d0
             jacPrimCons(2,1) = 0.0d0
             jacPrimCons(3,1) = 0.0d0
             jacPrimCons(4,1) = 0.0d0

             jacPrimCons(1,2) = - ux/nodes(1)
             jacPrimCons(2,2) = 1.0d0/nodes(1)
             jacPrimCons(3,2) = 0.0d0
             jacPrimCons(4,2) = 0.0d0

             jacPrimCons(1,3) = - uy/nodes(1)
             jacPrimCons(2,3) = 0.0d0
             jacPrimCons(3,3) = 1.0d0/nodes(1)
             jacPrimCons(4,3) = 0.0d0

             jacPrimCons(1,4) = dPdrho
             jacPrimCons(2,4) = ratio*ux
             jacPrimCons(3,4) = ratio*uy
             jacPrimCons(4,4) = -ratio

          else

             ratio = 3.0/(cv_r*(-3.0+nodes(1)))

             dPdrho = (-3.0*(nodes(2)**2*(-3.0 + 2.0*nodes(1)) + 
     $            nodes(2)**2*(-3.0 + 2.0*nodes(1)) + 
     $            2.0*nodes(1)**2*((3.0*(-6.0 + nodes(1)) + 
     $            2.0*cv_r*(-3.0 + nodes(1))**2)*nodes(1) - nodes(4))))/
     $            (2.0*cv_r*(-3.0 + nodes(1))**2*nodes(1)**2)

             jacPrimCons(1,1) = 1.0
             jacPrimCons(2,1) = 0.0
             jacPrimCons(3,1) = 0.0
             jacPrimCons(4,1) = 0.0

             jacPrimCons(1,2) = - ux/nodes(1)
             jacPrimCons(2,2) = 1.0/nodes(1)
             jacPrimCons(3,2) = 0.0
             jacPrimCons(4,2) = 0.0

             jacPrimCons(1,3) = - uy/nodes(1)
             jacPrimCons(2,3) = 0.0
             jacPrimCons(3,3) = 1.0/nodes(1)
             jacPrimCons(4,3) = 0.0

             jacPrimCons(1,4) = dPdrho
             jacPrimCons(2,4) = ratio*ux
             jacPrimCons(3,4) = ratio*uy
             jacPrimCons(4,4) = -ratio

          end if

        end function compute_jacobian_prim_to_cons


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for conservative to
        !> to primitive variables
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacConsPrim
        !> jacobian matrix for conservative to primitive
        !> variables \f$ \frac{\partial v}{\partial p} \f$
        !--------------------------------------------------------------
        function compute_jacobian_cons_to_prim(nodes)
     $     result(jacConsPrim)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacConsPrim


          real(rkind) :: ux
          real(rkind) :: uy
          real(rkind) :: dEdrho

          ux    = nodes(2)/nodes(1)
          uy    = nodes(3)/nodes(1)


          if(rkind.eq.8) then

             dEdrho = (-3.0d0*nodes(2)**2 - 3.0d0*nodes(3)**2 + 
     $            2.0d0*nodes(1)**2*((-3.0d0*(-6.0d0 + nodes(1)) - 
     $            2.0d0*cv_r*(-3.0d0 + nodes(1))**2)*nodes(1) + nodes(4)))/
     $            (2.0d0*(-3.0d0 + nodes(1))*nodes(1)**2)
             
             jacConsPrim(1,1) = 1.0d0
             jacConsPrim(2,1) = 0.0d0
             jacConsPrim(3,1) = 0.0d0
             jacConsPrim(4,1) = 0.0d0
                        
             jacConsPrim(1,2) = ux
             jacConsPrim(2,2) = nodes(1)
             jacConsPrim(3,2) = 0.0d0
             jacConsPrim(4,2) = 0.0d0
                        
             jacConsPrim(1,3) = uy
             jacConsPrim(2,3) = 0.0d0
             jacConsPrim(3,3) = nodes(1)
             jacConsPrim(4,3) = 0.0d0
                        
             jacConsPrim(1,4) = dEdrho
             jacConsPrim(2,4) = nodes(2)
             jacConsPrim(3,4) = nodes(3)
             jacConsPrim(4,4) = cv_r*(1.0d0-nodes(1)/3.0d0)

          else

             dEdrho = (-3.0*nodes(2)**2 - 3.0*nodes(3)**2 + 
     $            2.0*nodes(1)**2*((-3.0*(-6.0 + nodes(1)) - 
     $            2.0*cv_r*(-3.0 + nodes(1))**2)*nodes(1) + nodes(4)))/
     $            (2.0*(-3.0 + nodes(1))*nodes(1)**2)
             
             jacConsPrim(1,1) = 1.0
             jacConsPrim(2,1) = 0.0
             jacConsPrim(3,1) = 0.0
             jacConsPrim(4,1) = 0.0
                        
             jacConsPrim(1,2) = ux
             jacConsPrim(2,2) = nodes(1)
             jacConsPrim(3,2) = 0.0
             jacConsPrim(4,2) = 0.0
                        
             jacConsPrim(1,3) = uy
             jacConsPrim(2,3) = 0.0
             jacConsPrim(3,3) = nodes(1)
             jacConsPrim(4,3) = 0.0
                        
             jacConsPrim(1,4) = dEdrho
             jacConsPrim(2,4) = nodes(2)
             jacConsPrim(3,4) = nodes(3)
             jacConsPrim(4,4) = cv_r*(1.0-nodes(1)/3.0)

          end if

        end function compute_jacobian_cons_to_prim

      end module dim2d_prim_module
