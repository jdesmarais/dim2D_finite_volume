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
        use parameters_input , only : ne
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
     $       energy_inviscid_y_flux,
     $       speed_of_sound,
     $       compute_jacobian_prim_to_cons,
     $       compute_jacobian_cons_to_prim,
     $       cons_lodi_matrix_x,
     $       cons_lodi_matrix_y,
     $       compute_x_timedev_from_LODI_vector,
     $       compute_y_timedev_from_LODI_vector,
     $       compute_timedev_from_LODI_vectors

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


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the speed of sound given the data at the grid
        !> point location
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data where the speed of sound
        !> is evaluated
        !
        !>@param var
        !> \f$ \sqrt{\frac{\gamma P}{\rho}} = 
        !> \sqrt{\frac{\gamma(\gamma-1)}{\rho} \left[ \rho E - \frac{1}{2 \rho}
        !> \left( q_x^2 + q_y^2 right) \right]} \f$
        !---------------------------------------------------------------
        function speed_of_sound(nodes) result(var)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind)                            :: var
        
          
          if(rkind.eq.8) then
             var = Sqrt(
     $            gamma*(gamma-1.0d0)/nodes(1)*(
     $            nodes(4) - 0.5d0/nodes(1)*(nodes(2)**2+nodes(3)**2)))

          else
             var = Sqrt(
     $            gamma*(gamma-1.0)/nodes(1)*(
     $            nodes(4) - 0.5/nodes(1)*(nodes(2)**2+nodes(3)**2)))

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
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacPrimCons
        !> jacobian matrix for conservative to primitive
        !> variables \f$ \frac{\partial p}{\partial v} \f$
        !--------------------------------------------------------------
        function compute_jacobian_prim_to_cons(nodes)
     $     result(jacPrimCons)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacPrimCons

          real(rkind)                   :: ux
          real(rkind)                   :: uy

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)


          if(rkind.eq.8) then

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

             jacPrimCons(1,4) = 0.5d0*(gamma-1.0d0)*(ux**2+uy**2)
             jacPrimCons(2,4) = -(gamma-1.0d0)*ux
             jacPrimCons(3,4) = -(gamma-1.0d0)*uy
             jacPrimCons(4,4) = gamma-1.0d0             

          else

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

             jacPrimCons(1,4) = 0.5d0*(gamma-1.0d0)*(ux**2+uy**2)
             jacPrimCons(2,4) = -(gamma-1.0d0)*ux
             jacPrimCons(3,4) = -(gamma-1.0d0)*uy
             jacPrimCons(4,4) = gamma-1.0d0           

          end if

        end function compute_jacobian_prim_to_cons


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for conservative
        !> to primitive variables
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacConsPrim
        !> jacobian matrix for primitive to conservative
        !> variables \f$ \frac{\partial v}{\partial p} \f$
        !--------------------------------------------------------------
        function compute_jacobian_cons_to_prim(nodes)
     $     result(jacConsPrim)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacConsPrim

          real(rkind) :: ux
          real(rkind) :: uy

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)          


          if(rkind.eq.8) then
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

             jacConsPrim(1,4) = 0.5d0*(ux**2+uy**2)
             jacConsPrim(2,4) = nodes(1)*ux
             jacConsPrim(3,4) = nodes(1)*uy
             jacConsPrim(4,4) = 1.0d0/(gamma-1.0d0)

          else
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

             jacConsPrim(1,4) = 0.5*(ux**2+uy**2)
             jacConsPrim(2,4) = nodes(1)*ux
             jacConsPrim(3,4) = nodes(1)*uy
             jacConsPrim(4,4) = 1.0/(gamma-1.0)

          end if             
             
        end function compute_jacobian_cons_to_prim


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the conservative LODI matrix in the x-direction
        !
        !> @date
        !> 03_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return var
        !> conservative LODI matrix in the x-direction
        !--------------------------------------------------------------
        function cons_lodi_matrix_x(nodes) result(var)

            implicit none

            real(rkind), dimension(ne)   , intent(in) :: nodes
            real(rkind), dimension(ne,ne)             :: var

            real(rkind) :: u,v,c

            u = nodes(2)/nodes(1)
            v = nodes(3)/nodes(1)
            c = speed_of_sound(nodes)

            if(rkind.eq.8) then

               var(1,1) = -v/nodes(1)
               var(2,1) = 0.0d0
               var(3,1) = 1.0d0/nodes(1)
               var(4,1) = 0.0d0
               
               var(1,2) = c**2-0.5d0*(u**2+v**2)*(gamma-1.0d0)
               var(2,2) = u*(gamma-1.0d0)
               var(3,2) = v*(gamma-1.0d0)
               var(4,2) = 1.0d0-gamma
               
               var(1,3) = c*u+0.5d0*(u**2+v**2)*(gamma-1.0d0)
               var(2,3) = u*(1.0d0-gamma)-c
               var(3,3) = v*(1.0d0-gamma)
               var(4,3) = gamma-1.0d0
               
               var(1,4) = 0.5d0*(u**2+v**2)*(gamma-1.0d0)-c*u
               var(2,4) = u*(1.0d0-gamma)+c
               var(3,4) = v*(1.0d0-gamma)
               var(4,4) = gamma-1.0d0

            else
               
               var(1,1) = -v/nodes(1)
               var(2,1) = 0.0
               var(3,1) = 1.0/nodes(1)
               var(4,1) = 0.0
               
               var(1,2) = c**2-0.5*(u**2+v**2)*(gamma-1.0)
               var(2,2) = u*(gamma-1.0)
               var(3,2) = v*(gamma-1.0)
               var(4,2) = 1.0-gamma
               
               var(1,3) = c*u+0.5*(u**2+v**2)*(gamma-1.0)
               var(2,3) = u*(1.0-gamma)-c
               var(3,3) = v*(1.0-gamma)
               var(4,3) = gamma-1.0
               
               var(1,4) = 0.5*(u**2+v**2)*(gamma-1.0)-c*u
               var(2,4) = u*(1.0-gamma)+c
               var(3,4) = v*(1.0-gamma)
               var(4,4) = gamma-1.0

            end if
            
        end function cons_lodi_matrix_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the conservative LODI matrix in the y-direction
        !
        !> @date
        !> 03_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return var
        !> conservative LODI matrix in the y-direction
        !--------------------------------------------------------------
        function cons_lodi_matrix_y(nodes) result(var)

            implicit none

            real(rkind), dimension(ne)   , intent(in) :: nodes
            real(rkind), dimension(ne,ne)             :: var

            real(rkind) :: u,v,c

            u = nodes(2)/nodes(1)
            v = nodes(3)/nodes(1)
            c = speed_of_sound(nodes(:))


            if(rkind.eq.8) then
               
               var(1,1) = -u/nodes(1)
               var(2,1) = 1.0d0/nodes(1)
               var(3,1) = 0.0d0
               var(4,1) = 0.0d0
               
               var(1,2) = c**2-0.5d0*(u**2+v**2)*(gamma-1.0d0)
               var(2,2) = u*(gamma-1.0d0)
               var(3,2) = v*(gamma-1.0d0)
               var(4,2) = 1.0d0-gamma

               var(1,3) = c*v+0.5d0*(u**2+v**2)*(gamma-1.0d0)
               var(2,3) = u*(1.0d0-gamma)
               var(3,3) = v*(1.0d0-gamma)-c
               var(4,3) = gamma-1.0d0

               var(1,4) = 0.5d0*(u**2+v**2)*(gamma-1.0d0)-c*v
               var(2,4) = u*(1.0d0-gamma)
               var(3,4) = v*(1.0d0-gamma)+c
               var(4,4) = gamma-1.0d0

            else

               var(1,1) = -u/nodes(1)
               var(2,1) = 1.0/nodes(1)
               var(3,1) = 0.0
               var(4,1) = 0.0
               
               var(1,2) = c**2-0.5*(u**2+v**2)*(gamma-1.0)
               var(2,2) = u*(gamma-1.0)
               var(3,2) = v*(gamma-1.0)
               var(4,2) = 1.0-gamma

               var(1,3) = c*v+0.5*(u**2+v**2)*(gamma-1.0)
               var(2,3) = u*(1.0-gamma)
               var(3,3) = v*(1.0-gamma)-c
               var(4,3) = gamma-1.0

               var(1,4) = 0.5*(u**2+v**2)*(gamma-1.0)-c*v
               var(2,4) = u*(1.0-gamma)
               var(3,4) = v*(1.0-gamma)+c
               var(4,4) = gamma-1.0

            end if
            
        end function cons_lodi_matrix_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param lodi
        !> lodi vector in the y-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !--------------------------------------------------------------
        function compute_x_timedev_from_LODI_vector(
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
             timedev(2) = - 0.5d0/(nodes(1)*c)*(lodi(4)-lodi(3))
             timedev(3) = - lodi(1)
             timedev(4) = - 0.5d0*(lodi(3)+lodi(4))
          else
             timedev(1) = - 1.0/c**2*(lodi(2)+0.5*(lodi(3)+lodi(4)))
             timedev(2) = - 0.5/(nodes(1)*c)*(lodi(4)-lodi(3))
             timedev(3) = - lodi(1)
             timedev(4) = - 0.5*(lodi(3)+lodi(4))
          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)          

        end function compute_x_timedev_from_LODI_vector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the y-direction to the time derivatives of
        !> the conservative variables
        !
        !> @date
        !> 05_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param lodi
        !> LODI vector in the y-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives of
        !> the conservative variables
        !--------------------------------------------------------------
        function compute_y_timedev_from_LODI_vector(
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
             timedev(3) = - 0.5d0/(nodes(1)*c)*(lodi(4)-lodi(3))
             timedev(4) = - 0.5d0*(lodi(3)+lodi(4))
          else
             timedev(1) = - 1.0/c**2*(lodi(2)+0.5*(lodi(3)+lodi(4)))
             timedev(2) = - lodi(1)
             timedev(3) = - 0.5/(nodes(1)*c)*(lodi(4)-lodi(3))
             timedev(4) = - 0.5*(lodi(3)+lodi(4))
          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)

        end function compute_y_timedev_from_LODI_vector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms
        !> along the x- and y-direction to the time derivatives
        !> of the conservative variables
        !
        !> @date
        !> 08_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param lodi_x
        !> lodi vector in the x-direction
        !
        !>@param lodi_y
        !> lodi vector in the y-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic terms along the x- and
        !> y- directions to the time derivatives of the conservative
        !> variables
        !--------------------------------------------------------------
        function compute_timedev_from_LODI_vectors(
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
     $            - 0.5d0/(nodes(1)*c)*(lodi_x(4)-lodi_x(3))
     $            - lodi_y(1)

             timedev(3) =
     $            - lodi_x(1)
     $            - 0.5d0/(nodes(1)*c)*(lodi_y(4)-lodi_y(3))

             timedev(4) =
     $            - 0.5d0*(lodi_x(3)+lodi_x(4))
     $            - 0.5d0*(lodi_y(3)+lodi_y(4))

          else
             
             timedev(1) = 
     $            - 1.0/c**2*(lodi_x(2)+0.5*(lodi_x(3)+lodi_x(4))) 
     $            - 1.0/c**2*(lodi_y(2)+0.5*(lodi_y(3)+lodi_y(4)))

             timedev(2) =
     $            - 0.5/(nodes(1)*c)*(lodi_x(4)-lodi_x(3))
     $            - lodi_y(1)

             timedev(3) =
     $            - lodi_x(1)
     $            - 0.5/(nodes(1)*c)*(lodi_y(4)-lodi_y(3))

             timedev(4) =
     $            - 0.5*(lodi_x(3)+lodi_x(4))
     $            - 0.5*(lodi_y(3)+lodi_y(4))

          end if

          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the primitive
          !variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          timedev = MATMUL(timedev,jacConsPrim)

        end function compute_timedev_from_LODI_vectors

      end module ns2d_prim_module
