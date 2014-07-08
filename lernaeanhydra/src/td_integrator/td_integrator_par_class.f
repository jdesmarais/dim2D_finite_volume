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
      module td_integrator_par_class

        use bc_operators_par_class, only : bc_operators_par
        use cg_operators_class    , only : cg_operators
        use field_par_class       , only : field_par
        use parameters_kind       , only : rkind
        use dim2d_eq_class        , only : dim2d_eq
        use fv_operators_par_class, only : fv_operators_par

        implicit none

        private
        public :: td_integrator_par


        !> @class td_integrator
        !> abstract class encapsulating subroutines to integrate
        !> the governing equations using the time discretisation
        !> method and the boundary conditions on a parallel memory
        !> distributed system
        !>
        !> @param integrate
        !> integrate the computational field for dt
        !---------------------------------------------------------------
        type, abstract :: td_integrator_par

          contains

          procedure(integrate_proc), nopass, deferred :: integrate

        end type td_integrator_par


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
          !>@param sd
          !> space discretization operators
          !
          !>@param p
          !> physical model
          !
          !>@param td_par
          !> time discretisation operators
          !
          !>@param bc_par_used
          !> boundary conditions for a parallel memory distributed
          !> system
          !
          !>@param dt
          !> time step integrated
          !--------------------------------------------------------------
          subroutine integrate_proc(
     $       field_used, sd, p_model, td_par, bc_par_used, dt)

            import bc_operators_par
            import cg_operators
            import field_par
            import dim2d_eq
            import rkind
            import fv_operators_par

            type(field_par)       , intent(inout) :: field_used
            type(cg_operators)    , intent(in)    :: sd
            type(dim2d_eq)        , intent(in)    :: p_model
            type(fv_operators_par), intent(in)    :: td_par
            type(bc_operators_par), intent(in)    :: bc_par_used
            real(rkind)           , intent(in)    :: dt

          end subroutine integrate_proc

        end interface

      end module td_integrator_par_class
