      !> @file
      !> abstract class encapsulating subroutines to integrate
      !> the governing equations using the time discretisation
      !> method and the boundary conditions on a parallel memory
      !> distributed system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> abstract class encapsulating subroutines to integrate
      !> the governing equations using the time discretisation
      !> method and the boundary conditions on a parallel memory
      !> distributed system
      !
      !> @date
      !> 27_08_2013 - initial version - J.L. Desmarais
      !> 25_09_2013 - update for fv_operators_par - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_integrator_abstract_par_class

        use field_abstract_par_class, only : field_abstract_par
        use parameters_kind         , only : rkind

        implicit none

        private
        public :: td_integrator_abstract_par


        !> @class td_integrator
        !> abstract class encapsulating subroutines to integrate
        !> the governing equations using the time discretisation
        !> method and the boundary conditions on a parallel memory
        !> distributed system
        !>
        !> @param integrate
        !> integrate the computational field for dt
        !---------------------------------------------------------------
        type, abstract :: td_integrator_abstract_par

          contains

          procedure(integrate_proc), nopass, deferred :: integrate

        end type td_integrator_abstract_par


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to integrate the governing equations using
          !> space discretisation operators, physical model, 
          !> time discretisation operators and boundary conditions
          !> on a parallel memory distributed system
          !
          !> @date
          !> 27_08_2013 - initial version - J.L. Desmarais
          !> 25_09_2013 - update for fv_operators_par - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables and the
          !> cartesian communicator between the tiles
          !
          !>@param dt
          !> time step integrated
          !--------------------------------------------------------------
          subroutine integrate_proc(field_used, dt)

            import field_abstract_par
            import rkind

            class(field_abstract_par), intent(inout) :: field_used
            real(rkind)              , intent(in)    :: dt

          end subroutine integrate_proc

        end interface

      end module td_integrator_abstract_par_class
