      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for steady state
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for steady state
      !
      !> @date
      !> 11_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_parameters, only :
     $       cv_r

        use dim2d_prim_module, only :
     $       speed_of_sound

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       T0

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: ic

        
        !set the initial temperature in the field
        real(rkind), parameter :: u0_flow = 0.5d0*SQRT(2.0d0)
        real(rkind), parameter :: v0_flow = 0.5d0*SQRT(2.0d0)


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for phase separation
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
          procedure,   pass :: get_far_field

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
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param y_map
        !> map of y-coordinates
        !---------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map          

          
          !local variables
          real(rkind), dimension(ne) :: cst_nodes
          real(rkind)                :: t,x,y
          integer(ikind)             :: i,j
          integer                    :: k
          
          x = x_map(1)
          y = y_map(1)

          cst_nodes = get_far_field_cst(t,x,y)

          do k=1,ne
             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
                   nodes(i,j,k) = cst_nodes(k)
                end do
             end do
          end do

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s
          real(rkind) :: c

          side_s = side

          c = speed_of_sound_liquid()

          var = u0_flow/c

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s
          real(rkind) :: c

          side_s = side

          c = speed_of_sound_liquid()

          var = v0_flow/c

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

          var = u0_flow

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

          var = v0_flow

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

          var = T0

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
          real(rkind) :: mass

          t_s = t
          x_s = x
          y_s = y

          mass = get_mass_density_liquid(T0)

          if(rkind.eq.8) then
             var = 8.0d0*mass*T0/(3.0d0-mass) - 3.0d0*mass**2
          else
             var = 8.0*mass*T0/(3.0-mass) - 3.0*mass**2
          end if

        end function get_P_out


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the governing variables imposed in the far field
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@return var
        !> governing variables in the far-field
        !--------------------------------------------------------------
        function get_far_field(this,t,x,y) result(var)

          implicit none

          class(ic)     , intent(in) :: this
          real(rkind)   , intent(in) :: t
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var


          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: d_liq
          
          t_s = t
          x_s = x
          y_s = y


          d_liq = get_mass_density_liquid(T0)

          if(rkind.eq.8) then

             var(1) = d_liq
             var(2) = d_liq*u0_flow
             var(3) = d_liq*v0_flow
             var(4) = 0.5d0*d_liq*((u0_flow)**2 + (v0_flow)**2) + d_liq*(8.0d0/3.0d0*cv_r*T0-3.0d0*d_liq)

          else

             var(1) = d_liq
             var(2) = d_liq*u0_flow
             var(3) = d_liq*v0_flow
             var(4) = 0.5*d_liq*((u0_flow)**2 + (v0_flow)**2) + d_liq*(8.0/3.0*cv_r*T0-3.0*d_liq)

          end if

        end function get_far_field


        function get_far_field_cst(t,x,y) result(var)

          implicit none

          real(rkind)   , intent(in) :: t
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var


          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: d_liq
          
          t_s = t
          x_s = x
          y_s = y


          d_liq = get_mass_density_liquid(T0)

          if(rkind.eq.8) then

             var(1) = d_liq
             var(2) = d_liq*u0_flow
             var(3) = d_liq*v0_flow
             var(4) = 0.5d0*d_liq*((u0_flow)**2 + (v0_flow)**2) + d_liq*(8.0d0/3.0d0*cv_r*T0-3.0d0*d_liq)

          else

             var(1) = d_liq
             var(2) = d_liq*u0_flow
             var(3) = d_liq*v0_flow
             var(4) = 0.5*d_liq*((u0_flow)**2 + (v0_flow)**2) + d_liq*(8.0/3.0*cv_r*T0-3.0*d_liq)

          end if

        end function get_far_field_cst


        function speed_of_sound_liquid() result(c)

          implicit none

          real(rkind), dimension(ne) :: nodes
          real(rkind)                :: c

          real(rkind) :: t,x,y

          nodes = get_far_field_cst(t,x,y)
          c = speed_of_sound(nodes)

        end function speed_of_sound_liquid

      end module ic_class
