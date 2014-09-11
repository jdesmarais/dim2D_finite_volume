      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain
      !
      !> @date
      !> 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_abstract_class

        use parameters_kind, only :
     $     rkind

        implicit none

        
        private
        public :: ic_abstract


        !> @class ic_abstract
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain
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
        type, abstract :: ic_abstract

          contains

          procedure(ic_proc) , nopass, deferred :: apply_ic
          procedure(far_proc), nopass, deferred :: get_mach_ux_infty
          procedure(far_proc), nopass, deferred :: get_mach_uy_infty
          procedure(var_proc), nopass, deferred :: get_u_in
          procedure(var_proc), nopass, deferred :: get_v_in
          procedure(var_proc), nopass, deferred :: get_T_in
          procedure(var_proc), nopass, deferred :: get_P_out

        end type ic_abstract


        abstract interface
        
          !apply the initial conditions for the governing variables
          subroutine ic_proc(nodes,x_map,y_map)

            import rkind

            real(rkind), dimension(:,:,:), intent(inout) :: nodes
            real(rkind), dimension(:)    , intent(in)    :: x_map
            real(rkind), dimension(:)    , intent(in)    :: y_map

          end subroutine ic_proc


          !get the variable enforced at the edge of the
          !computational domain
          function far_proc(side) result(var)

            import rkind

            logical    , intent(in) :: side
            real(rkind)             :: var

          end function far_proc

        
          !get the variable enforced at the edge of the
          !computational domain
          function var_proc(t,x,y) result(var)

            import rkind

            real(rkind), intent(in) :: t
            real(rkind), intent(in) :: x
            real(rkind), intent(in) :: y
            real(rkind)             :: var

          end function var_proc

        end interface

      end module ic_abstract_class
