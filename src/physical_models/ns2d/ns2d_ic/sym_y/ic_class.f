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
     $       left,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       flow_direction

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic


        real(rkind) :: u0_x_flow = 1.0d0
        real(rkind) :: u0_y_flow = 0.0d0
        real(rkind) :: v0_x_flow = 0.0d0
        real(rkind) :: v0_y_flow = 1.0d0


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
          

          real(rkind) :: x_s
          real(rkind) :: y_s
          real(rkind) :: nodes_s

          nodes_s = nodes(1,1,1)
          x_s = x_map(1)
          y_s = y_map(1)

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
               if(side.eqv.left) then
                  var = -v0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)
               else
                  var = v0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)
               end if

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
               var =-u0_y_flow

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
               if(y.gt.0) then
                  var =-v0_x_flow
               else
                  var = v0_x_flow
               end if

            case(y_direction)
               if(y.gt.0) then
                  var =-v0_y_flow
               else
                  var = v0_y_flow
               end if

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

      end module ic_class
