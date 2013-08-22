      !> @file
      !> module encapsulating subroutines to initialize the MPI
      !> derived types used to exchange data between tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to initialize the MPI
      !> derived types used to exchange data between tiles
      !
      !> @date
      ! 22_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_mg_derived_types

        use parameters_constant, only : N,S,E,W
        use parameters_kind    , only : ikind, rkind
        use mpi

        implicit none

        private
        public :: ini_mpi_derived_types

        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the mpi derived types
        !> for the exchange of data between the tiles
        !
        !> @date
        !> 21_08_2012 - initial version - J.L. Desmarais
        !
        !>@param comm_2d
        !> MPI Cartesian communicator betwwen the tiles
        !
        !>@param tile_nx
        !> number of gridpoints along the x-axis for the tile
        !
        !>@param tile_ny
        !> number of gridpoints along the y-axis for the tile
        !
        !> @param com_rank
        !> rank of the processor computing the neighbouring tiles
        !--------------------------------------------------------------
        subroutine ini_mpi_derived_types(
     $       tile_nx, tile_ny, tile_ne, bc_size,
     $       com_send, com_recv)

          implicit none

          integer(ikind)       , intent(in)  :: tile_nx
          integer(ikind)       , intent(in)  :: tile_ny
          integer              , intent(in)  :: tile_ne
          integer              , intent(in)  :: bc_size
          integer, dimension(4), intent(out) :: com_send
          integer, dimension(4), intent(out) :: com_recv


          !< local variables
          integer(ikind) :: i_min,i_max,j_min,j_max
          integer        :: k

          
          !< create the derived types for sending data
          !> by determining the subarray borders of the
          !> data sent. Then the subarray derived type
          !> is committed
          do k=1,4

             !< determine the subarray borders
             select case(k)

               case(N)
                  i_min = 1 
                  i_max = tile_nx
                  j_min = tile_ny-2*bc_size+1
                  j_max = tile_ny-bc_size

               case(S)
                  i_min = 1
                  i_max = tile_nx
                  j_min = bc_size+1
                  j_max = 2*bc_size

               case(E)
                  i_min = tile_nx-2*bc_size+1
                  i_max = tile_nx-bc_size
                  j_min = bc_size+1
                  j_max = tile_ny-bc_size

               case(W)
                  i_min = bc_size+1
                  i_max = 2*bc_size
                  j_min = bc_size+1
                  j_max = tile_ny-bc_size

               case default
                  stop 'mpi_mg_ini_derived_types: wrong case'

             end select


             !< create the derived type
             com_send(k)= create_mpi_subarray(
     $            tile_nx, tile_ny, tile_ne,
     $            i_min, i_max, j_min, j_max)

          end do


          !< create the derived types for receving data
          !> by determining the subarray borders of the
          !> data sent. Then the subarray derived type
          !> is committed
          do k=1,4

             !< determine the subarray borders
             select case(k)

               case(N)
                  i_min = 1
                  i_max = tile_nx
                  j_min = tile_ny-bc_size+1
                  j_max = tile_ny

               case(S)
                  i_min = 1
                  i_max = tile_nx
                  j_min = 1
                  j_max = bc_size

               case(E)
                  i_min = tile_nx-bc_size+1
                  i_max = tile_nx
                  j_min = bc_size+1
                  j_max = tile_ny-bc_size

               case(W)
                  i_min = 1
                  i_max = bc_size
                  j_min = bc_size+1
                  j_max = tile_ny-bc_size

               case default
                  stop 'mpi_mg_ini_derived_types: wrong case'

             end select


             !< create the derived type
             com_recv(k)= create_mpi_subarray(
     $            tile_nx, tile_ny, tile_ne,
     $            i_min, i_max, j_min, j_max)

          end do


        end subroutine ini_mpi_derived_types
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function commiting a subarray mpi derived type
        !> from the borders of the subarray
        !
        !> @date
        !> 22_08_2013 - initial version - J.L. Desmarais
        !
        !>@param tile_nx
        !> number of gridpoints along the x-axis for the tile
        !
        !>@param tile_ny
        !> number of gridpoints along the y-axis for the tile
        !
        !>@param tile_ne
        !> number of governing equations for the tile
        !
        !> @param i_min
        !> indice along the x-axis for the SW and NW borders
        !
        !> @param i_max
        !> indice along the x-axis for the SE and NE borders
        !
        !> @param j_min
        !> indice along the y-axis for the SE and SW borders
        !
        !> @param j_max
        !> indice along the y-axis for the NE and NW borders
        !
        !> @param type_subarray
        !> integer identifying the new MPI derived type for the
        !> subarray
        !--------------------------------------------------------------
        function create_mpi_subarray(
     $     tile_nx, tile_ny, tile_ne,
     $     i_min, i_max, j_min, j_max)
     $     result(type_subarray)

          implicit none

          integer(ikind), intent(in) :: tile_nx, tile_ny
          integer       , intent(in) :: tile_ne
          integer(ikind), intent(in) :: i_min, i_max, j_min, j_max
          integer                    :: type_subarray

          !local variable
          integer(ikind)               :: sub_lines_nb
          integer(ikind)               :: sub_columns_nb
          integer(ikind), dimension(3) :: profile_array
          integer(ikind), dimension(3) :: profile_subarray
          integer(ikind), dimension(3) :: beginning_point
          
          integer(ikind) :: ierror


          sub_lines_nb       = i_max-i_min+1
          sub_columns_nb     = j_max-j_min+1

          
          profile_array(1)   = tile_nx !column number
          profile_array(2)   = tile_ny !line number
          profile_array(3)   = tile_ne !var nb


          profile_subarray(1)= sub_lines_nb
          profile_subarray(2)= sub_columns_nb
          profile_subarray(3)= tile_ne


          beginning_point(1) = i_min-1
          beginning_point(2) = j_min-1
          beginning_point(3) = 0

          select case(rkind)
            case(4)

               call MPI_TYPE_CREATE_SUBARRAY(
     $              3,
     $              profile_array,
     $              profile_subarray,
     $              beginning_point, MPI_ORDER_FORTRAN,
     $              MPI_REAL, type_subarray,
     $              ierror)

            case(8)

               call MPI_TYPE_CREATE_SUBARRAY(
     $              3,
     $              profile_array,
     $              profile_subarray,
     $              beginning_point, MPI_ORDER_FORTRAN,
     $              MPI_DOUBLE_PRECISION, type_subarray,
     $              ierror)
               
            case default
               print *, 'mpi_mg_bc_ini_derived_types'
               print *, 'create_derived_types'
               stop 'rkind not compatible with MPI kinds'
                  
          end select


          if(ierror.ne.MPI_SUCCESS) then
             stop 'create subarray: error in creating the subarray'
          end if


          call MPI_TYPE_COMMIT(type_subarray, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'create subarray: error in submitting
     $            the derived type'
          end if

        end function create_mpi_subarray

      end module mpi_mg_derived_types



