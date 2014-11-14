      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a convected vortex
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a convected vortex
      !
      !> @date
      !> 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use ns2d_parameters, only :
     $     gamma, mach_infty

        use ic_abstract_class, only :
     $     ic_abstract

        use parameters_constant, only :
     $       x_direction,
     $       y_direction,
     $       xy_direction

        use parameters_input, only :
     $       flow_direction

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic


        real(rkind) :: u0_x_flow  = 1.0d0
        real(rkind) :: u0_y_flow  = 0.0d0
        real(rkind) :: u0_xy_flow = 0.5d0*SQRT(2.0d0)
        real(rkind) :: v0_x_flow  = 0.0d0
        real(rkind) :: v0_y_flow  = 1.0d0
        real(rkind) :: v0_xy_flow = 0.5d0*SQRT(2.0d0)


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for a convected vortex
        !
        !> @param apply_initial_conditions
        !> set the initial conditions
        !
        !> @param get_mach_ux_infty
        !> get the mach number along the x-direction in the far field
        !
        !> @param get_mach_uy_infty
        !> get the mach number along the y-direction in the far field
        !
        !> @param get_u_in
        !> get the x-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_v_in
        !> get the y-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_T_in
        !> get the temperature at the edge of the computational
        !> domain
        !
        !> @param get_P_out
        !> get the pressure at the edge of the computational domain
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          contains

          procedure, nopass :: apply_ic
          procedure, nopass :: get_mach_ux_infty
          procedure, nopass :: get_mach_uy_infty
          procedure, nopass :: get_u_in
          procedure, nopass :: get_v_in
          procedure, nopass :: get_T_in
          procedure, nopass :: get_P_out

        end type ic


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
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:)      , intent(inout) :: nodes
          real(rkind), dimension(:)          , intent(in)    :: x_map
          real(rkind), dimension(:)          , intent(in)    :: y_map
          

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

          select case(flow_direction)

            case(x_direction)
               u0_mean_flow = u0_x_flow
               v0_mean_flow = v0_x_flow

            case(y_direction)
               u0_mean_flow = u0_y_flow
               v0_mean_flow = v0_y_flow

            case(xy_direction)
               u0_mean_flow = u0_xy_flow
               v0_mean_flow = v0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select


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

        end subroutine apply_ic        


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          select case(flow_direction)

            case(x_direction)
               var = u0_x_flow*mach_infty/Sqrt(u0_x_flow**2+v0_x_flow**2)
               
            case(y_direction)
               var = u0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)

            case(xy_direction)
               var = u0_xy_flow*mach_infty/Sqrt(u0_xy_flow**2+v0_xy_flow**2)
               
            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_mach_ux_infty'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          select case(flow_direction)

            case(x_direction)
               var = v0_x_flow*mach_infty/Sqrt(u0_x_flow**2+v0_x_flow**2)
               
            case(y_direction)
               var = v0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)

            case(xy_direction)
               var = v0_xy_flow*mach_infty/Sqrt(u0_xy_flow**2+v0_xy_flow**2)

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_mach_uy_infty'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_mach_uy_infty


        !get the x-component of the velocity enforced
        !at the edge of the computational domain
        function get_u_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          select case(flow_direction)

            case(x_direction)
               var = u0_x_flow
               
            case(y_direction)
               var = u0_y_flow

            case(xy_direction)
               var = u0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_u_in'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_u_in


        !get the y-component of the velocity enforced
        !at the edge of the computational domain
        function get_v_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          select case(flow_direction)

            case(x_direction)
               var = v0_x_flow

            case(y_direction)
               var = v0_y_flow

            case(xy_direction)
               var = v0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_v_in'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_v_in


        !get the temperature enforced at the edge of the
        !computational domain
        function get_T_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 1.0d0
          else
             var = 1.0
          end if

        end function get_T_in


        !get the pressure enforced at the edge of the
        !computational domain
        function get_P_out(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 1.0d0/(gamma*mach_infty**2)
          else
             var = 1.0/(gamma*mach_infty**2)
          end if

        end function get_P_out


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

      end module ic_class
