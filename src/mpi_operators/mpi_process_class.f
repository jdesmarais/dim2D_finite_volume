      !> @file
      !> class encapsulating subroutines for the initilization
      !> and the finalization of mpi processes
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines for the initilization
      !> and the finalization of mpi processes
      !
      !> @date
      ! 01_05_2013 - initial version - J.L. Desmarais
      ! 21_08_2013 - adaptation      - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_process_class

        use mpi

        use parameters_constant, only :
     $       periodic_xy_choice,
     $       reflection_xy_choice,
     $       wall_xy_choice,
     $       wall_x_reflection_y_choice,
     $       hedstrom_xy_choice,
     $       hedstrom_xy_corners_choice,
     $       hedstrom_x_reflection_y_choice,
     $       poinsot_xy_choice,
     $       yoolodato_xy_choice

        use parameters_input, only :
     $       npx,npy,bc_choice

        implicit none


        private
        public :: mpi_process


        !> @class mpi_process
        !> class encapsulating subroutines for the initialization and
        !> finalization of mpi processes
        !
        !> @param ini_mpi
        !> initialize the mpi process
        !
        !> @param ini_cartesian_communicator
        !> initialize a 2D-grid cartesian communicator between the
        !> computational tiles
        !
        !> @param finalize_mpi
        !> finalize the mpi process
        !---------------------------------------------------------------
        type :: mpi_process

          contains

          procedure, nopass :: ini_mpi
          procedure, nopass :: ini_cartesian_communicator
          procedure, nopass :: finalize_mpi

        end type mpi_process


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialization of the mpi process and check if the total
        !> number of processors is equal to the total number of 
        !> tiles
        !
        !> @date
        !> 01_05_2013 - initial version - J.L. Desmarais
        !> 21_08_2013 - adaptation      - J.L. Desmarais
        !--------------------------------------------------------------
        subroutine ini_mpi()

          implicit none

          !< local variables
          integer :: ierror
          integer :: total_nb_tiles
          integer :: total_nb_procs


          !< initialize the mpi process
          call MPI_INIT(ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'mpi_process : MPI_INIT failed'
          end if
          
          
          !< get the total number of tiles
          total_nb_tiles = npx*npy
          
          
          !< get the total number of processors
          call MPI_COMM_SIZE(
     $         MPI_COMM_WORLD,
     $         total_nb_procs,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'mpi_process : MPI_COMM_SIZE failed'
          end if
          
          
          !< compare the total number of processors
          !< with the total number of tiles
          if(total_nb_procs.ne.total_nb_tiles) then
             print '(''nb_tiles: '', I5)', total_nb_tiles
             print '(''nb_procs: '', I5)', total_nb_procs
             print *, 'please use total_nb_tiles=total_nb_procs'
             call finalize_mpi()
             stop
          end if
          
        end subroutine ini_mpi


        !initialize the cartesian communicator between the
        !tiles to have a 2D grid configuration
        subroutine ini_cartesian_communicator(comm_2d, usr_rank)

          implicit none

          integer, intent(out) :: comm_2d
          integer, intent(out) :: usr_rank

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

            case(reflection_xy_choice,
     $           wall_xy_choice,
     $           wall_x_reflection_y_choice,
     $           hedstrom_xy_choice,
     $           hedstrom_xy_corners_choice,
     $           hedstrom_x_reflection_y_choice,
     $           poinsot_xy_choice,
     $           yoolodato_xy_choice)
               periods(1)     = .false.
               periods(2)     = .false.

            case default
               print *, 'mpi_process_class:'
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
     $         comm_2d,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call finalize_mpi()
             stop 'mpi_process_class: MPI_CART_CREATE failed'
          end if


          !< get the rank of the processor computing the tile
          call MPI_COMM_RANK(comm_2d, usr_rank, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call finalize_mpi()
             stop 'mpi_process_comm: MPI_COMM_RANK failed'
          end if

        end subroutine ini_cartesian_communicator



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> finalization of the mpi processes
        !
        !> @date
        !> 01_05_2013 - initial version - J.L. Desmarais
        !> 21_08_2013 - adaptation      - J.L. Desmarais
        !--------------------------------------------------------------
        subroutine finalize_mpi()

          implicit none

          !< local variables
          integer :: ierror

          call MPI_FINALIZE(ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'mpi_process : MPI_FINALIZE failed'
          end if

        end subroutine finalize_mpi

      end module mpi_process_class
