      !> @file
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_class

        use field_abstract_class , only : field_abstract
        use parameters_kind      , only : rkind
        use td_integrator_class  , only : td_integrator

        implicit none


        private
        public :: field


        !> @class field
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !---------------------------------------------------------------
        type, extends(field_abstract) :: field

          type(td_integrator) :: td_integrator_used

          contains

          procedure, pass :: integrate

        end type field


        contains


        !integrate the field in time from t to t+dt
        subroutine integrate(this, dt)

          implicit none

          class(field), intent(inout) :: this
          real(rkind) , intent(in)    :: dt

          call this%td_integrator_used%integrate(this,dt)

        end subroutine integrate

      end module field_class
