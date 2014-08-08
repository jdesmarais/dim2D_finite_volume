      !> @file
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the Navier-Stokes governing
      !> equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the Navier-Stokes governing
      !> equations
      !
      !> @date
      !> 08_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ns2d_prim_module

        use ns2d_parameters  , only : gamma, mach_infty
        use interface_primary, only : get_primary_var
        use parameters_kind  , only : ikind, rkind

        implicit none

        private
        public ::
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       qx_inviscid_x_flux,
     $       qy_inviscid_y_flux,
     $       qxy_transport,
     $       energy_inviscid_x_flux,
     $       energy_inviscid_y_flux


        contains

        
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the mass density \f$ \rho \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        function pressure(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var= (gamma-1.0d0)*(nodes(i,j,4)-0.5d0/nodes(i,j,1)*(
     $            nodes(i,j,2)**2+nodes(i,j,3)**2))
          else
             var= (gamma-1.0)*(nodes(i,j,4)-0.5/nodes(i,j,1)*(
     $            nodes(i,j,2)**2+nodes(i,j,3)**2))
          end if

        end function pressure


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the temperature
        !> \f$ \frac{\gamma(\gamma-1) M_{\infty}^2}{\rho}
        !> \left( \rho E - \frac{1}{2\rho} (q_x^2 + q_y^2) \right)
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> \f$ T \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function temperature(
     $     nodes,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=gamma*mach_infty**2*pressure(nodes,i,j)/nodes(i,j,1)

        end function temperature       


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_x component for the inviscid flux
        !> along the x-axis \f$ \rho u_x^2 + \pressure \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> \f$ \rho u_x^2 + \pressure \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qx_inviscid_x_flux(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var= nodes(i,j,2)**2/nodes(i,j,1) + pressure(nodes,i,j)

        end function qx_inviscid_x_flux


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_y component for the inviscid flux
        !> along the y-axis \f$ \rho u_y^2 + \pressure \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> \f$ \rho u_y^2 + P \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function qy_inviscid_y_flux(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)**2/nodes(i,j,1) + pressure(nodes,i,j)

        end function qy_inviscid_y_flux


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the momentum_x transported along the y-axis
        !> \f$ \rho u_x u_y \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        function qxy_transport(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)*nodes(i,j,3)/nodes(i,j,1)

        end function qxy_transport


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the energy-component of the insvicid flux
        !> along the x-axis \f$ (\rho E + P) u_x \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> \f$ (\rho E + P) u_x \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function energy_inviscid_x_flux(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=(nodes(i,j,4)+pressure(nodes,i,j))*nodes(i,j,2)/nodes(i,j,1)

        end function energy_inviscid_x_flux


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the energy-component of the insvicid flux
        !> along the y-axis \f$ (\rho E + P) u_y \f$
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> \f$ (\rho E + P) u_y \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function energy_inviscid_y_flux(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=(nodes(i,j,4)+pressure(nodes,i,j))*nodes(i,j,3)/nodes(i,j,1)

        end function energy_inviscid_y_flux

      end module ns2d_prim_module
