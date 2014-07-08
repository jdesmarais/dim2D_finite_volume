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

        use parameters_input, only : npx,npy
        use mpi

        implicit none


        private
        public :: mpi_process


        !> @class mpi_process
        !> class encapsulating subroutines for the initialization and
        !> finalization of mpi processes
        !>
        !> @param ini_mpi
        !> initialize the mpi process
        !>
        !> @param finalize_mpi
        !> finalize the mpi process
        !---------------------------------------------------------------
        type :: mpi_process

          contains

          procedure, nopass :: ini_mpi
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
