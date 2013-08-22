      !> @file
      !> module encapsulating subroutine to get the rank
      !> of the processors computing the neighboring tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutine to get the rank
      !> of the processors computing the neighboring tiles
      !
      !> @date
      ! 23_08_2012 - initial version - J.L. Desmarais
      ! 21_08_2013 - adaptation      - J.L. Desmarais
      !--------------------------------------------------------------
      module mpi_mg_neighbours

        use mpi
        use parameters_constant, only : N,S,E,W

        implicit none

        private
        public :: ini_neighbours_proc_id


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine identifying the processors computing
        !> the neighbouring tiles
        !
        !> @date
        !> 23_08_2012 - initial version - J.L. Desmarais
        !> 21_08_2012 - adaptation      - J.L. Desmarais
        !
        !>@param comm_2d
        !> MPI Cartesian communicator betwwen the tiles
        !
        !>@param tiles_nb
        !> tile number in the x and y directions
        !
        !> @param com_rank
        !> rank of the processor computing the neighbouring tiles
        !--------------------------------------------------------------
        subroutine ini_neighbours_proc_id(
     $       comm_2D,
     $       nb_tiles,
     $       com_rank)

          implicit none

          !< dummy arguments
          integer              , intent(in) :: comm_2D
          integer, dimension(2), intent(in) :: nb_tiles
          integer, dimension(4), intent(out):: com_rank


          !< local variables
          !
          !> dims_nb:              dimensions nb of the Cartesian grid
          !>                       for the communications
          !
          !> tile_proc_rank:       rank of the processor computing
          !>                       the tile
          !
          !> tile_cart_coord:      Cartesian coordinates of the tile
          !
          !> neighbour_cart_coord: Cartesian coordinates of the
          !>                       neighbouring tile
          !
          !> ierror:               success indicator for MPI calls
          !-----------------------------------------------------------
          integer               :: dims_nb             
          integer               :: tile_proc_rank      
          integer, dimension(2) :: tile_cart_coord     
          integer, dimension(2) :: neighbour_cart_coord
          integer               :: ierror              
          integer               :: k

          dims_nb =2


          !< compute the rank of the processor computing the tile
          call MPI_COMM_RANK(comm_2D, tile_proc_rank, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'ini_neighbouring_tiles_processor_id:
     $            get processor rank for tile'
          end if


          !< compute the Cartesian coordinates of the tile
          call MPI_CART_COORDS(comm_2D, tile_proc_rank,
     $         dims_nb, tile_cart_coord, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'ini_neighbouring_tiles_processor_id:
     $            MPI_CART_COORDS failed'
          end if


          !< compute the Cartesian coordinates of the neighbouring tiles
          !> and get the rank of the processors computing them
          do k=1,4
             
             neighbour_cart_coord = tile_cart_coord

             !< get the Cartesian coordinates of the neighbouring tile
             select case(k)

               case(N)
                  neighbour_cart_coord(2)=neighbour_cart_coord(2)+1
                  neighbour_cart_coord(2)=mod(
     $                 neighbour_cart_coord(2), nb_tiles(2))

               case(S)
                  neighbour_cart_coord(2)=neighbour_cart_coord(2)-1
                  neighbour_cart_coord(2)=mod(
     $                 neighbour_cart_coord(2), nb_tiles(2))

                  if(neighbour_cart_coord(2)<0) then
                     neighbour_cart_coord(2)=nb_tiles(2)+
     $                    neighbour_cart_coord(2)
                  end if

               case(E)
                  neighbour_cart_coord(1)=neighbour_cart_coord(1)+1
                  neighbour_cart_coord(1)=mod(
     $                 neighbour_cart_coord(1), nb_tiles(1))

               case(W)
                  neighbour_cart_coord(1)=neighbour_cart_coord(1)-1
                  neighbour_cart_coord(1)=mod(
     $                 neighbour_cart_coord(1),nb_tiles(1))

                  if(neighbour_cart_coord(1)<0) then
                     neighbour_cart_coord(1)=nb_tiles(1)+
     $                    neighbour_cart_coord(1)
                  end if

             end select

            
             !< get the rank of the processor computing the
             !> neighbouring tile
             call MPI_CART_RANK(
     $            comm_2D,
     $            neighbour_cart_coord,
     $            com_rank(k),
     $            ierror)
             if(ierror.ne.MPI_SUCCESS) then
                stop 'ini_neighbouring_tiles_processor_id:
     $               MPI_CART_RANK failed'
             end if
             
          end do

        end subroutine ini_neighbours_proc_id

      end module mpi_mg_neighbours
