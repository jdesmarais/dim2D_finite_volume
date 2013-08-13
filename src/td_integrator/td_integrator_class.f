      !> @file
      !> abstract class encapsulating subroutines to integrate
      !> the governing equations using the time discretisation
      !> method and the boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> in the governing equations
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_integrator_class

        use field_class , only : field
        use parameters_kind, only : rkind

        implicit none


        type, abstract :: td_integrator

          contains

          procedure(integrate_proc), nopass, deferred :: integrate

        end type td_integrator


        abstract interface

          subroutine integrate_proc(field_bc_used, dt)

            import field_bc
            import rkind

            class(field_bc), intent(inout) :: field_bc_used
            real(rkind)    , intent(in)    :: dt

          end subroutine integrate_proc

        end interface

      end module td_integrator_class
