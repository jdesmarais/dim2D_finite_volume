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
      module td_integrator_class

        use cg_operators_class, only : cg_operators
        use field_class       , only : field
        use parameters_kind   , only : rkind
        use phy_model_eq_class, only : phy_model_eq
        use td_operators_class, only : td_operators

        implicit none


        !> @class td_integrator
        !> abstract class encapsulating subroutines to integrate
        !> the governing equations using the time discretisation
        !> method and the boundary conditions
        !>
        !> @param integrate
        !> integrate the computational field for dt
        !---------------------------------------------------------------
        type, abstract :: td_integrator

          contains

          procedure(integrate_proc), nopass, deferred :: integrate

        end type td_integrator


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
          !>@param sd
          !> space discretization operators
          !
          !>@param p
          !> physical model
          !
          !>@param td
          !> time discretisation operators
          !
          !>@param dt
          !> time step integrated
          !--------------------------------------------------------------
          subroutine integrate_proc(field_used, sd, p_model, td, dt)


            import cg_operators
            import field
            import phy_model_eq
            import rkind
            import td_operators

            class(field)       , intent(inout) :: field_used
            type(cg_operators) , intent(in)    :: sd
            class(phy_model_eq), intent(in)    :: p_model
            class(td_operators), intent(in)    :: td
            real(rkind)        , intent(in)    :: dt

          end subroutine integrate_proc

        end interface

      end module td_integrator_class
