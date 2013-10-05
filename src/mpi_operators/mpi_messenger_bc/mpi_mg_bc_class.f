      !> @file
      !> class encapsulating subroutines to exchange information
      !> between the tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to exchange information
      !> between the tiles
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_mg_bc_class

        use cg_operators_class  , only : cg_operators
        use field_par_class     , only : field_par
        use mpi
        use mpi_mg_neighbours   , only : ini_neighbours_proc_id
        use mpi_mg_derived_types, only : ini_mpi_derived_types
        use parameters_input           , only : nx,ny,ne,npx,npy
        
        implicit none

        private
        public :: mpi_mg_bc


        !> @class mpi_mg_bc
        !> class encapsulating subroutines to exchange information
        !> between the tiles
        !>
        !> @param com_rank
        !> table containing the rank of the processors
        !> computing the neighbouring tiles
        !>
        !> @param com_send
        !> table containing the derived type to send
        !> data to the neighbouring processors
        !>
        !> @param com_recv
        !> table containing the derived type to send
        !> data to the neighbouring processors
        !---------------------------------------------------------------
        type :: mpi_mg_bc

          integer, dimension(4) :: com_rank                 
          integer, dimension(4) :: com_send                 
          integer, dimension(4) :: com_recv

          contains

          procedure, pass :: initialize

        end type mpi_mg_bc


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the attributes of
        !> the 'mpi_mg_bc' object
        !
        !> @date
        !> 21_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> 'mpi_mg_bc' object initialized
        !--------------------------------------------------------------
        subroutine initialize(this,f_used,s)
        
          implicit none

          class(mpi_mg_bc), intent(inout) :: this
          class(field_par)       , intent(in)    :: f_used
          type(cg_operators)     , intent(in)    :: s

          !< local variables
          integer, dimension(2) :: nb_tiles


          !< initialization of local variables
          nb_tiles(1) = npx
          nb_tiles(2) = npy

          !< initialize the 'com_rank' attribute by checking
          !> the rank of the processors computing the neighbours
          call ini_neighbours_proc_id(
     $         f_used%comm_2d, nb_tiles, this%com_rank)


          !< initialize the 'com_send' and 'com_recv' attributes
          !> by creating the MPI derived types needed
          call ini_mpi_derived_types(
     $         nx, ny, ne,
     $         this%com_send, this%com_recv)

        end subroutine initialize

      end module mpi_mg_bc_class
      
