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
     $       gradient_proc

        use n_coords_module, only :
     $       get_n1_coord,
     $       get_n2_coord

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
     $       velocity_n1,
     $       velocity_n2,
     $       classical_pressure,
     $       classical_pressure_local,
     $       pressure,
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
     $       compute_jacobian_cons_to_prim,
     $       compute_x_timedev_from_LODI_vector_dim2d,
     $       compute_y_timedev_from_LODI_vector_dim2d,
     $       compute_timedev_from_LODI_vectors_dim2d


        contains

        
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the mass density, \f$ \rho \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
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
        !> compute the momentum density along the x-axis,
        !> \f$ \rho u \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ \rho u \f$ evaluated at [i,j]
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
        !> compute the momentum density along the y-axis,
        !> \f$ \rho v \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ \rho v \f$ evaluated at [i,j]
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
        !> compute the total energy density, \f$\rho E\f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
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
        !> compute the velocity along the x-axis,
        !> \f$ u \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ u \f$ evaluated at [i,j]
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
        !> compute the velocity along the y-axis,
        !> \f$ v \f$
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ v \f$ evaluated at [i,j]
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
        !> compute the velocity along the n1-axis
        !> \f[ u_{n_1} = \frac{\sqrt{2}}{2} (u - v) \f]
        !
        !> @date
        !> 05_02_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ u_{n_1} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_n1(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = get_n1_coord(
     $         nodes(i,j,2)/nodes(i,j,1),
     $         nodes(i,j,3)/nodes(i,j,1))

        end function velocity_n1


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the velocity along the n2-axis
        !> \f[ u_{n_2} = \frac{\sqrt{2}}{2} (u + v) \f]
        !
        !> @date
        !> 05_02_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ u_{n_2} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_n2(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = get_n2_coord(
     $         nodes(i,j,2)/nodes(i,j,1),
     $         nodes(i,j,3)/nodes(i,j,1))

        end function velocity_n2


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical pressure
        !> \f[ P^{\textrm{cl}} = \frac{3}{(3-\rho) c_v}
        !> \left( \rho E - \frac{1}{2} \rho (u^2 + v^2) + 3\rho^2
        !> \right) - 3 \rho^2 \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ P^{\textrm{cl}} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function classical_pressure(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=classical_pressure_local(nodes(i,j,:))

        end function classical_pressure



        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical pressure
        !> \f[ P^{\textrm{cl}} = \frac{3}{(3-\rho) c_v}
        !> \left( \rho E - \frac{1}{2} \rho (u^2 + v^2) + 3\rho^2
        !> \right) - 3 \rho^2 \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@return
        !> \f$ P^{\textrm{cl}} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function classical_pressure_local(nodes) result(pressure)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind)                            :: pressure

          if(rkind.eq.8) then
             pressure = 3.0d0/((3.0d0-nodes(1))*cv_r)*(
     $            nodes(4)
     $            - 0.5d0*nodes(1)*(
     $            (nodes(2)/nodes(1))**2+
     $            (nodes(3)/nodes(1))**2)
     $            + 3.0d0*nodes(1)**2)
     $            - 3.0d0*nodes(1)**2
             
          else

             pressure = 3.0/((3.0-nodes(1))*cv_r)*(
     $            nodes(4)
     $            - 0.5*nodes(1)*(
     $            (nodes(2)/nodes(1))**2+
     $            (nodes(3)/nodes(1))**2)
     $            + 3.0*nodes(1)**2)
     $            - 3.0*nodes(1)**2

          end if

        end function classical_pressure_local


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the thermodynamic pressure from
        !> the Van der Waals state equation
        !> \f[ P = \frac{8 \rho T}{3-\rho} - 3 \rho^2 \f]
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature, \f$ T \f$
        !
        !>@param mass_density
        !> mass_density, \f$ \rho \f$
        !
        !>@return
        !> pressure, \f$ P \f$
        !--------------------------------------------------------------
        function pressure(temperature,mass_density)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind), intent(in) :: mass_density
          real(rkind)             :: pressure

          pressure =
     $         8.0d0*mass_density*temperature/(3.0d0-mass_density) -
     $         3.0d0*mass_density**2

        end function pressure


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the classical part of the effective temperature
        !> \f[ T^{\textrm{cl}}_{\textrm{eff}} = \frac{1}{\rho}
        !> \left( \rho E - \frac{1}{2} \rho (u^2 + v^2)
        !>  + 3 \rho^2 \right) \f]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ T_{\textrm{eff}}^{\textrm{cl}} \f$ evaluated at [i,j]
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
        !> compute the capillary part of the effective temperature
        !> \f[ T_{\textrm{eff}}^{\textrm{ca}} =
        !> \frac{1}{2 \rho} {|\nabla \rho|}^2 \f]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing along the x-direction
        !
        !>@param dy
        !> grid spacing along the y-direction
        !
        !>@param gradient_x
        !> procedure to compute the gradient in the x-direction
        !
        !>@param gradient_y
        !> procedure to compute the gradient in the y-direction
        !
        !>@return
        !> \f$ T_{\textrm{eff}}^{\textrm{ca}} \f$ evaluated at [i,j]
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
          procedure(gradient_proc)                  :: gradient_x
          procedure(gradient_proc)                  :: gradient_y
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
        !> \f[ T_{\textrm{eff}} = \frac{1}{\rho}
        !> \left( \rho E - \frac{1}{2} \rho (u^2 + v^2)
        !> - \frac{1}{2 We} {|\nabla \rho|}^2 + 3 \rho^2 \right) \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing along the x-direction
        !
        !>@param dy
        !> grid spacing along the y-direction
        !
        !>@param gradient_x
        !> procedure to compute the gradient in the x-direction
        !
        !>@param gradient_y
        !> procedure to compute the gradient in the y-direction
        !
        !>@return
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
          procedure(gradient_proc)                  :: gradient_x
          procedure(gradient_proc)                  :: gradient_y
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
        !> \f[ \left(\frac{3}{(3-\rho) c_v}
        !> \left[ \rho E - \frac{1}{2} \rho (u^2 + v^2) + 3\rho^2
        !> \right] - 3 \rho^2 \right) u \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
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
        !> \f[ \left(\frac{3}{(3-\rho) c_v}
        !> \left[ \rho E - \frac{1}{2} \rho (u^2 + v^2) + 3\rho^2
        !> \right] - 3 \rho^2 \right) v \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
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
        !> \f[ \rho u^2 \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho u^2 \f$ evaluated at [i,j]
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
        !> \f[ \rho v u \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho v u \f$ evaluated at [i,j]
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
        !> \f[ \rho u v \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho u v \f$ evaluated at [i,j]
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
        !> \f[ \rho v^2 \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho v^2 \f$ evaluated at [i,j]
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
        !> \f[ \rho E u \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho E u \f$ evaluated at [i,j]
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
        !> \f[ \rho E v \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \rho E v \f$ evaluated at [i,j]
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
        !> \f[ \frac{1}{3-\rho} \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
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
        !> along the x-axis \f[ \frac{u}{3-\rho} \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \frac{u}{3-\rho} \f$ evaluated at [i,j]
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
        !> along the y-axis \f[ \frac{v}{3-\rho} \f]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ \frac{v}{3-\rho} \f$ evaluated at [i,j]
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
        !> compute the speed of sound given the conservative
        !> variables at the grid point location
        !> \f[ c = \frac{\sqrt{2} \rho}{a d} \f] where
        !> \f[ a = \frac{1}{\sqrt{1+b^2}}\f]
        !> \f[ b = \sqrt{\frac{P + 3 \rho^2}{c_V ( P+ \rho^2 (-3+2 \rho))}}\f]
        !> \f[ d = \sqrt{\rho^3} \sqrt{\frac{3-\rho}{P+\rho^2(-3+2\rho)}} \f]
        
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !>
        !>@return
        !> speed of sound, \f$ c \f$
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

          P = classical_pressure_local(nodes)

          if(rkind.eq.8) then

             e = P+nodes(1)**2*(-3.0d0+2.0d0*nodes(1))
             b = SQRT((P+3.0d0*nodes(1)**2)/(cv_r*e))
             a = 1.0d0/SQRT(1.0d0+b**2)
             d = SQRT(nodes(1)**3)*SQRT((3.0d0-nodes(1))/e)

             var = SQRT(3.0d0)*nodes(1)/(a*d)

          else

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
        !> \f[ J^p_v =
        !> \begin{pmatrix}
        !>    1 & 0 & 0 & 0 \\\
        !>    - u/\rho & 1/\rho & 0 & 0 \\\
        !>    - v/\rho & 1/\rho & 0 & 0 \\\
        !>    -\displaystyle{\frac{3(u^2 + v^2 + 12 \rho) + 2 c_v(P + 9(\rho-2)\rho) }{2 c_v (-3 + \rho)}} &
        !>    \displaystyle{\frac{3 u}{c_v (-3 + \rho)}} &
        !>    \displaystyle{\frac{3 v}{c_v (-3 + \rho)}} &
        !>  - \displaystyle{\frac{3}{c_v (-3 + \rho)}  }
        !> \end{pmatrix}\f]
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@return
        !> jacobian matrix for primitive to conservative
        !> variables
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
        !> \f[ J^v_p =
        !> \begin{pmatrix}
        !>    1 & 0 & 0 & 0 \\\
        !>    u & \rho & 0 & 0 \\\
        !>    v & 0    & \rho & 0 \\\
        !>    \displaystyle{\frac{1}{2}(u^2 + v^2) - 6 \rho - \frac{1}{3} c_v ( P + 9(\rho-2)\rho)} &
        !>    \rho u &
        !>    \rho v &
        !>    \displaystyle{c_v \left( 1 - \frac{\rho}{3} \right)}
        !> \end{pmatrix}\f]
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@return
        !> jacobian matrix for conservative to primitive
        !> variables
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !> \f[ \frac{\partial v}{\partial t} = - J^v_p \cdot
        !> \begin{pmatrix}
        !>     \displaystyle{\frac{1}{c^2} \left[ \mathscr{L}_2^1 + \frac{1}{2}( \mathscr{L}_3^1 + \mathscr{L}_4^1 ) \right]} \\\
        !>     \displaystyle{\frac{1}{2 \rho c} ( \mathscr{L}_4^1 - \mathscr{L}_3^1 )} \\\
        !>     \mathscr{L}_1^1 \\\
        !>     \displaystyle{\frac{1}{2} ( \mathscr{L}_3^1 + \mathscr{L}_4^1 )}
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param lodi
        !> lodi vector in the x-direction, \f$ (\mathscr{L}_1^1, \mathscr{L}_2^1, \mathscr{L}_3^1, \mathscr{L}_4^1)^T \f$
        !
        !>@return
        !> contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !--------------------------------------------------------------
        function compute_x_timedev_from_LODI_vector_dim2d(
     $     nodes, lodi) result(timedev)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: lodi
          real(rkind), dimension(ne)             :: timedev

          real(rkind)                   :: c
          real(rkind), dimension(ne,ne) :: jacConsPrim

          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then
             timedev(1) = - 1.0d0/c**2*(lodi(2)+0.5d0*(lodi(3)+lodi(4)))
             timedev(2) = - 0.5d0/(nodes(1)*c)*(-lodi(3)+lodi(4))
             timedev(3) = - lodi(1)
             timedev(4) = - 0.5d0*(lodi(3)+lodi(4))
          else
             timedev(1) = - 1.0/c**2*(lodi(2)+0.5*(lodi(3)+lodi(4)))
             timedev(2) = - 0.5/(nodes(1)*c)*(-lodi(3)+lodi(4))
             timedev(3) = - lodi(1)
             timedev(4) = - 0.5*(lodi(3)+lodi(4))
          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)

        end function compute_x_timedev_from_LODI_vector_dim2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the y-direction to the time derivatives of
        !> the conservative variables
        !> \f[ \frac{\partial v}{\partial t} = - J^v_p \cdot
        !> \begin{pmatrix}
        !>     \displaystyle{\frac{1}{c^2} \left[ \mathscr{L}_2^2 + \frac{1}{2}( \mathscr{L}_3^2 + \mathscr{L}_4^2 ) \right]} \\\
        !>     \mathscr{L}_1^2 \\\
        !>     \displaystyle{\frac{1}{2 \rho c} ( \mathscr{L}_4^2 - \mathscr{L}_3^2 )} \\\
        !>     \displaystyle{\frac{1}{2} ( \mathscr{L}_3^2 + \mathscr{L}_4^2 )}
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param lodi
        !> LODI vector in the y-direction, \f$ (\mathscr{L}_1^2, \mathscr{L}_2^2, \mathscr{L}_3^2, \mathscr{L}_4^2)^T \f$
        !
        !>@return
        !> contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !--------------------------------------------------------------
        function compute_y_timedev_from_LODI_vector_dim2d(
     $     nodes, lodi) result(timedev)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: lodi
          real(rkind), dimension(ne)             :: timedev

          real(rkind)                   :: c
          real(rkind), dimension(ne,ne) :: jacConsPrim

          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then
             timedev(1) = - 1.0d0/c**2*(lodi(2)+0.5d0*(lodi(3)+lodi(4)))
             timedev(2) = - lodi(1)
             timedev(3) = - 0.5d0/(nodes(1)*c)*(-lodi(3)+lodi(4))
             timedev(4) = - 0.5d0*(lodi(3)+lodi(4))
          else
             timedev(1) = - 1.0/c**2*(lodi(2)+0.5*(lodi(3)+lodi(4)))
             timedev(2) = - lodi(1)
             timedev(3) = - 0.5/(nodes(1)*c)*(-lodi(3)+lodi(4))
             timedev(4) = - 0.5*(lodi(3)+lodi(4))
          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)

        end function compute_y_timedev_from_LODI_vector_dim2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the x- and y-direction to the time derivatives
        !> of the conservative variables
        !> \f[ \frac{\partial v}{\partial t} = - J^v_p \cdot
        !> \begin{pmatrix}
        !>     \displaystyle{\frac{1}{c^2} \left[ \mathscr{L}_2^1 + \mathscr{L}_2^2 + \frac{1}{2}( \mathscr{L}_3^1 + \mathscr{L}_3^2 + \mathscr{L}_4^1 + \mathscr{L}_4^2 ) \right]} \\\
        !>     \displaystyle{\frac{1}{2 \rho c} ( \mathscr{L}_4^1 - \mathscr{L}_3^1 )} + \mathscr{L}_1^2 \\\
        !>     \mathscr{L}_1^1 + \displaystyle{\frac{1}{2 \rho c} ( \mathscr{L}_4^2 - \mathscr{L}_3^2 )} \\\
        !>     \displaystyle{\frac{1}{2} ( \mathscr{L}_3^1 + \mathscr{L}_3^2 + \mathscr{L}_4^1 + \mathscr{L}_4^2 )}
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 08_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the conservative variables (\f$ \rho, \rho u, \rho v, \rho E \f$)
        !
        !>@param lodi_x
        !> lodi vector in the x-direction, \f$ (\mathscr{L}_1^1, \mathscr{L}_2^1, \mathscr{L}_3^1, \mathscr{L}_4^1)^T \f$
        !
        !>@param lodi_y
        !> lodi vector in the y-direction, \f$ (\mathscr{L}_1^2, \mathscr{L}_2^2, \mathscr{L}_3^2, \mathscr{L}_4^2)^T \f$
        !
        !>@return
        !> contribution of the hyperbolic terms along the x- and
        !> y- directions to the time derivatives of the conservative
        !> variables
        !--------------------------------------------------------------
        function compute_timedev_from_LODI_vectors_dim2d(
     $     nodes, lodi_x, lodi_y) result(timedev)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: lodi_x
          real(rkind), dimension(ne), intent(in) :: lodi_y
          real(rkind), dimension(ne)             :: timedev

          real(rkind)                   :: c
          real(rkind), dimension(ne,ne) :: jacConsPrim

          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then

             timedev(1) = 
     $            - 1.0d0/c**2*(lodi_x(2)+0.5d0*(lodi_x(3)+lodi_x(4))) 
     $            - 1.0d0/c**2*(lodi_y(2)+0.5d0*(lodi_y(3)+lodi_y(4)))

             timedev(2) =
     $            - 0.5d0/(nodes(1)*c)*(-lodi_x(3)+lodi_x(4))
     $            - lodi_y(1)

             timedev(3) =
     $            - lodi_x(1)
     $            - 0.5d0/(nodes(1)*c)*(-lodi_y(3)+lodi_y(4))

             timedev(4) =
     $            - 0.5d0*(lodi_x(3)+lodi_x(4))
     $            - 0.5d0*(lodi_y(3)+lodi_y(4))

          else
             
             timedev(1) = 
     $            - 1.0/c**2*(lodi_x(2)+0.5*(lodi_x(3)+lodi_x(4))) 
     $            - 1.0/c**2*(lodi_y(2)+0.5*(lodi_y(3)+lodi_y(4)))

             timedev(2) =
     $            - 0.5/(nodes(1)*c)*(-lodi_x(3)+lodi_x(4))
     $            - lodi_y(1)

             timedev(3) =
     $            - lodi_x(1)
     $            - 0.5/(nodes(1)*c)*(-lodi_y(3)+lodi_y(4))

             timedev(4) =
     $            - 0.5*(lodi_x(3)+lodi_x(4))
     $            - 0.5*(lodi_y(3)+lodi_y(4))

          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)

        end function compute_timedev_from_LODI_vectors_dim2d

      end module dim2d_prim_module
