      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a steady state
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a steady state
      !
      !> @date
      !> 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use ic_abstract_class, only :
     $     ic_abstract

        use parameters_constant, only :
     $       x_direction,
     $       y_direction

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain at steady state
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
          procedure, nopass :: get_mach_ux_infty => get_far
          procedure, nopass :: get_mach_uy_infty => get_far
          procedure, nopass :: get_u_in  => get_var
          procedure, nopass :: get_v_in  => get_var
          procedure, nopass :: get_T_in  => get_var
          procedure, nopass :: get_P_out => get_var

        end type ic


        contains


        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          
          integer(ikind) :: i,j
          real(rkind)    :: xc,yc

          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                
                nodes(i,j,1) =  1.0
                nodes(i,j,2) =  0.0
                nodes(i,j,3) =  0.0
                nodes(i,j,4) =  2.0

             end do
          end do

          xc = x_map(1)
          yc = y_map(1)

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_far() result(var)

          implicit none

          real(rkind) :: var
          
          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_far


        !get the variable enforced at the edge of the
        !computational domain
        function get_var(t,x,y) result(var)

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
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_var

      end module ic_class
