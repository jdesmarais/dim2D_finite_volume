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

        use parameters_constant, only : periodic_xy_choice,
     $                                  reflection_xy_choice
        use parameters_input   , only : npx, npy, bc_choice
        use field_class        , only : field
        use mpi

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
        type, extends(field) :: field_par

          integer :: comm_2d

          contains

          procedure, pass :: ini_cartesian_communicator

        end type field_par


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the comm_2d and usr_rank
        !> attributes of the tile
        !
        !> @date
        !> 21_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> 'field_par' object initialized
        !
        !>@param nb_tiles
        !> number of tiles in the x and y directions
        !
        !>@param bc_chosen
        !> boundary conditions chosen
        !--------------------------------------------------------------
        subroutine ini_cartesian_communicator(this)

          implicit none

          !< dummy arguments
          class(field_par), intent(inout) :: this

          
          !< local variables
          integer, dimension(2) :: nb_tiles
          integer               :: comm_old
          integer               :: ndims
          logical, dimension(2) :: periods
          logical               :: reorganisation
          integer               :: ierror


          !< set the characteristics of the cartesian grid
          nb_tiles(1)    = npx
          nb_tiles(2)    = npy

          select case(bc_choice)

            case(periodic_xy_choice)
               periods(1)     = .true.
               periods(2)     = .true.

            case(reflection_xy_choice)
               periods(1)     = .false.
               periods(2)     = .false.

            case default
               print *, 'bc_choice not implemented in'
               stop 'splitting the field into tiles'

          end select

          comm_old       = MPI_COMM_WORLD
          ndims          = 2
          reorganisation = .true.


          !< create Cartesian communicator
          call MPI_CART_CREATE(
     $         comm_old,
     $         ndims,
     $         nb_tiles,
     $         periods,
     $         reorganisation,
     $         this%comm_2d,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'field_par_class: MPI_CART_CREATE failed'
          end if

        end subroutine ini_cartesian_communicator

      end module field_par_class
