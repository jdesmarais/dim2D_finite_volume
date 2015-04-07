      !> @file
      !> class extending the 'field' class to integrate mpi
      !> attributes that will allows the tile to communicate
      !> with its neighbours
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class extending the 'field' class to integrate mpi
      !> attributes that will allows the tile to communicate
      !> with its neighbours
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_par_class
      
        use field_abstract_par_class, only : field_abstract_par
        use mpi
        use parameters_kind         , only : rkind
        use td_integrator_par_class , only : td_integrator_par

        implicit none

        private
        public :: field_par


        !> @class field
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !>
        !> @param comm_2d
        !> attribute identifying the mpi main communicator between the
        !> tiles
        !>
        !> @param usr_rank
        !> attribute identifying the processor computing the tile
        !---------------------------------------------------------------
        type, extends(field_abstract_par) :: field_par

          type(td_integrator_par) :: td_integrator_used

          contains

          procedure, pass :: integrate

        end type field_par


        contains


        !integrate the field in time from t to t+dt
        subroutine integrate(this, dt)

          implicit none

          class(field_par), intent(inout) :: this
          real(rkind)     , intent(in)    :: dt

          call this%td_integrator_used%integrate(this, dt)
          
        end subroutine integrate

      end module field_par_class
