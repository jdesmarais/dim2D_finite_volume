      !> @file
      !> class encapsulating the main tables for the variables and the
      !> coordinates and the operators to integrate them in time
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates and the operators to integrate them in time
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

          type(td_integrator) :: td_integrator_used !<@brief time integration operators

          contains

          procedure, pass :: integrate !<@brief integrate the computational domain in time

        end type field


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> integrate the governing variables in time from t to t+dt
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        subroutine integrate(this, dt)

          implicit none

          class(field), intent(inout) :: this
          real(rkind) , intent(in)    :: dt

          call this%td_integrator_used%integrate(this,dt)

          this%time = this%time + dt

        end subroutine integrate

      end module field_class
