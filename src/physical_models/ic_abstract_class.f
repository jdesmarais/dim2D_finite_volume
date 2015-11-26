      !> @file
      !> abstract class encapsulating interfaces to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> abstract class encapsulating interfaces to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain
      !
      !> @date
      !> 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_abstract_class

        use parameters_input, only :
     $     ne

        use parameters_kind, only :
     $     rkind

        implicit none

        
        private
        public :: ic_abstract


        !> @class ic_abstract
        !> abstract class encapsulating interfaces to apply the initial
        !> conditions and to evaluate the conditions enforced
        !> at the edge of the computational domain
        !---------------------------------------------------------------
        type, abstract :: ic_abstract

          contains

          procedure(ic_proc)       , nopass, deferred :: apply_ic          !<@brief set the initial conditions                                                 
          procedure(far_proc)      , nopass, deferred :: get_mach_ux_infty !<@brief get the Mach number along the x-direction in the far field                 
          procedure(far_proc)      , nopass, deferred :: get_mach_uy_infty !<@brief get the Mach number along the y-direction in the far field                 
          procedure(var_proc)      , nopass, deferred :: get_u_in          !<@brief get the x-component of the velocity at the edge of the computational domain
          procedure(var_proc)      , nopass, deferred :: get_v_in          !<@brief get the y-component of the velocity at the edge of the computational domain
          procedure(var_proc)      , nopass, deferred :: get_T_in          !<@brief get the temperature at the edge of the computational domain                
          procedure(var_proc)      , nopass, deferred :: get_P_out         !<@brief get the pressure at the edge of the computational domain                   
          procedure(far_field_proc),   pass, deferred :: get_far_field     !<@brief get the governing variables imposed at the edge of the computational domain

        end type ic_abstract


        abstract interface
        
          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to apply the initial conditions
          !> to the computational domain
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param this
          !> physical model
          !
          !>@param nodes
          !> array with the grid point data    
          !
          !>@param x_map
          !> array with the x-coordinates
          !
          !>@param y_map
          !> array with the y-coordinates                
          !--------------------------------------------------------------
          subroutine ic_proc(nodes,x_map,y_map)

            import rkind

            real(rkind), dimension(:,:,:), intent(inout) :: nodes
            real(rkind), dimension(:)    , intent(in)    :: x_map
            real(rkind), dimension(:)    , intent(in)    :: y_map

          end subroutine ic_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to obtain the value of the variable
          !> imposed at the edge of the computational domain
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param side
          !> left or right side for the x-direction,
          !> top or bottom for the y-direction
          !
          !>@return
          !> variable imposed at the edge of the computational
          !> domain
          !--------------------------------------------------------------
          function far_proc(side) result(var)

            import rkind

            logical    , intent(in) :: side
            real(rkind)             :: var

          end function far_proc

        
          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to obtain the value of the variable
          !> imposed at the edge of the computational domain
          !> depending on time and coordinates
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
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
          !>@return
          !> variable imposed at the edge of the computational
          !> domain
          !--------------------------------------------------------------
          function var_proc(t,x,y) result(var)

            import rkind

            real(rkind), intent(in) :: t
            real(rkind), intent(in) :: x
            real(rkind), intent(in) :: y
            real(rkind)             :: var

          end function var_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to obtain the value of the variables
          !> imposed at the edge of the computational domain
          !> depending on time and coordinates as well as the
          !> state of the object
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param this
          !> object encapsulating the initial conditions and
          !> the state of the conditions imposed in the far-field
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
          !>@return
          !> variable imposed at the edge of the computational
          !> domain
          !--------------------------------------------------------------
          function far_field_proc(this,t,x,y) result(var)

            import ic_abstract
            import ne
            import rkind

            class(ic_abstract)        , intent(in) :: this
            real(rkind)               , intent(in) :: t
            real(rkind)               , intent(in) :: x
            real(rkind)               , intent(in) :: y
            real(rkind), dimension(ne)             :: var

          end function far_field_proc

        end interface

      end module ic_abstract_class
