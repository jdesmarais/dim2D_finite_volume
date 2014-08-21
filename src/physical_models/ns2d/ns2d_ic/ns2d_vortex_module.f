      !> @file
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Navier-Stokes equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Navier-Stokes equations
      !
      !> @date
      !> 12_08_2014 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module ns2d_vortex_module

        use ns2d_parameters , only : gamma, mach_infty
        use parameters_kind , only : ikind, rkind

        implicit none

        private
        public :: apply_vortex_ic

        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !---------------------------------------------------------------
        subroutine apply_vortex_ic(nodes,x_map,y_map,velocity_meanflow)

          implicit none

          real(rkind), dimension(:,:,:)      , intent(inout) :: nodes
          real(rkind), dimension(:)          , intent(in)    :: x_map
          real(rkind), dimension(:)          , intent(in)    :: y_map
          real(rkind), dimension(2), optional, intent(in)    :: velocity_meanflow
          

          real(rkind)    :: mass_infty
          integer(ikind) :: i,j
          real(rkind)    :: l,amp,R
          real(rkind)    :: u0_mean_flow, v0_mean_flow
          real(rkind)    :: p_infty


          !as the scaling for mass density in the NS equations
          !is the mass density in the far field, we have
          !mass_infty=1.0
          mass_infty = 1.0d0

          
          !as the scaling for velocity in the NS equations
          !is the velocity norm in the far field, we should
          !have u0_mean_flow**2+v0_meanflow**2=1.0

          if(present(velocity_meanflow)) then
             u0_mean_flow = velocity_meanflow(1)
             v0_mean_flow = velocity_meanflow(2)
          else
             if(rkind.eq.8) then
                u0_mean_flow = 0.0d0
                v0_mean_flow = 0.0d0
             else
                u0_mean_flow = 0.0
                v0_mean_flow = 0.0
             end if
          end if

          !vortex located at the center of the computational domain
          !the vortex characteristics scales with the size of the
          !computational domain
          l   = 0.25d0*(x_map(size(x_map,1))-x_map(1))
          amp = -0.5d0*l
          R   = 0.15d0*l


          !computation of the pressure in the far field
          if((u0_mean_flow**2+v0_mean_flow**2).le.(1.0e-8)) then
             P_infty = 1.0d0
          else
             P_infty = 1.0d0/(gamma*mach_infty**2)
          end if

          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                
                !mass density: same as in the far field
                nodes(i,j,1) =  mass_infty

                !momentum-x
                nodes(i,j,2) = 
     $               nodes(i,j,1)*(
     $               u0_mean_flow+
     $               get_vortex_velocity_x(
     $                  x_map(i),y_map(j),
     $                  nodes(i,j,1),
     $                  amp,R)
     $               )

                !momentum-y
                nodes(i,j,3) =
     $               nodes(i,j,1)*(
     $               v0_mean_flow+
     $               get_vortex_velocity_y(
     $                  x_map(i),y_map(j),
     $                  nodes(i,j,1),
     $                  amp,R)
     $               )

                !total energy
                nodes(i,j,4) =  0.5d0/nodes(i,j,1)*(
     $               nodes(i,j,2)**2 + nodes(i,j,3)**2)+
     $               get_pressure(
     $               x_map(i),y_map(j),
     $               nodes(i,j,1),amp,R,P_infty)/(gamma-1.0d0)

             end do
          end do

        end subroutine apply_vortex_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x1
        !> coordinate along the x-axis
        !
        !>@param x2
        !> coordinate along the y-axis
        !
        !>@param mass
        !> mass density at (x1,x2)
        !
        !>@param amp
        !> vortex peak amplitude
        !
        !>@param R
        !> radius of the vortex
        !
        !>@return ux
        !> velocity along the x-axis at (x1,x2)
        !---------------------------------------------------------------
        function get_vortex_velocity_x(x1,x2,mass,amp,R)
     $     result(ux)

          implicit none

          real(rkind), intent(in) :: x1
          real(rkind), intent(in) :: x2
          real(rkind), intent(in) :: mass
          real(rkind), intent(in) :: amp
          real(rkind), intent(in) :: R
          real(rkind)             :: ux
          
          if(rkind.eq.8) then
             ux = - 1.0d0/(mass*R**2)*amp*x2*Exp(-(x1**2+x2**2)/(2.0d0*R**2))
          else
             ux = - 1.0/(mass*R**2)*amp*x2*Exp(-(x1**2+x2**2)/(2.0*R**2))
          end if

        end function get_vortex_velocity_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x1
        !> coordinate along the x-axis
        !
        !>@param x2
        !> coordinate along the y-axis
        !
        !>@param mass
        !> mass density at (x1,x2)
        !
        !>@param amp
        !> vortex peak amplitude
        !
        !>@param R
        !> radius of the vortex
        !
        !>@return uy
        !> velocity along the y-axis at (x1,x2)
        !---------------------------------------------------------------
        function get_vortex_velocity_y(x1,x2,mass,amp,R)
     $     result(uy)

          implicit none

          real(rkind), intent(in) :: x1
          real(rkind), intent(in) :: x2
          real(rkind), intent(in) :: mass
          real(rkind), intent(in) :: amp
          real(rkind), intent(in) :: R
          real(rkind)             :: uy
          
          if(rkind.eq.8) then
             uy =   1.0d0/(mass*R**2)*amp*x1*Exp(-(x1**2+x2**2)/(2.0d0*R**2))
          else
             uy =   1.0/(mass*R**2)*amp*x1*Exp(-(x1**2+x2**2)/(2.0*R**2))
          end if          

        end function get_vortex_velocity_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x1
        !> coordinate along the x-axis
        !
        !>@param x2
        !> coordinate along the y-axis
        !
        !>@param mass
        !> mass density at (x1,x2)
        !
        !>@param amp
        !> vortex peak amplitude
        !
        !>@param R
        !> radius of the vortex
        !
        !>@param P_infty
        !> pressure in the far field
        !
        !>@return P
        !> pressure at (x1,x2)
        !---------------------------------------------------------------
        function get_pressure(x1,x2,mass,amp,R,P_infty)
     $     result(P)

          implicit none

          real(rkind), intent(in) :: x1
          real(rkind), intent(in) :: x2
          real(rkind), intent(in) :: mass
          real(rkind), intent(in) :: amp
          real(rkind), intent(in) :: R
          real(rkind), intent(in) :: P_infty
          real(rkind)             :: P
          
          if(rkind.eq.8) then
             P = P_infty + mass*amp**2/R**2*Exp(-(x1**2+x2**2)/(2.0d0*R**2))
          else
             P = P_infty + mass*amp**2/R**2*Exp(-(x1**2+x2**2)/(2.0*R**2))
          end if

        end function get_pressure        

      end module ns2d_vortex_module
