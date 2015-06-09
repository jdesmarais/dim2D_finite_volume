      !> @file
      !> abstract class encapsulating subroutines to integrate
      !> the governing equations using the time discretisation
      !> method and the boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> abstract class encapsulating subroutines to integrate
      !> the governing equations using the time discretisation
      !> method and the boundary conditions
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_integrator_abstract_class

        use field_abstract_class, only :
     $       field_abstract

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: td_integrator_abstract


        !> @class td_integrator_abstract
        !> abstract class encapsulating subroutines to integrate
        !> the governing equations using the time discretisation
        !> method and the boundary conditions
        !>
        !> @param integrate
        !> integrate the computational field for dt
        !---------------------------------------------------------------
        type, abstract :: td_integrator_abstract

          contains

          procedure(timeInt_proc)    , nopass, deferred :: integrate
          procedure(timeInt_ext_proc), nopass, deferred :: integrate_ext
          

        end type td_integrator_abstract


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to integrate the governing equations using
          !> space discretisation operators, physical model, 
          !> time discretisation operators and boundary conditions
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param dt
          !> time step integrated
          !--------------------------------------------------------------
          subroutine timeInt_proc(field_used,dt)

            import field_abstract
            import rkind

            class(field_abstract), intent(inout) :: field_used
            real(rkind)          , intent(in)    :: dt

          end subroutine timeInt_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to integrate the governing equations on
          !> the computationaal domain and its extension
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param dt
          !> time step integrated
          !--------------------------------------------------------------
          subroutine timeInt_ext_proc(field_used,dt)

            import field_abstract
            import rkind

            class(field_abstract), intent(inout) :: field_used
            real(rkind)          , intent(in)    :: dt

          end subroutine timeInt_ext_proc

        end interface

      end module td_integrator_abstract_class
