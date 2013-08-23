      !> @file
      !> module encapsulating subroutines to choose the type
      !> of procedures used by each tile when computing the
      !> boundary layers
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to choose the type
      !> of procedures used by each tile when computing the
      !> boundary layers
      !
      !> @date
      ! 22_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_mg_ini_bc_proc

        use mpi
        use parameters_constant, only : N,S,E,W,
     $                                  periodic_xy_choice,
     $                                  only_compute_proc,
     $                                  compute_and_exchange_proc,
     $                                  only_exchange_proc,
     $                                  x_direction, y_direction
        use parameters_input   , only : npx,npy,bc_choice
        
        implicit none

        private
        public :: ini_bc_procedures


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the procedures computing
        !> the boundary layers
        !
        !> @date
        !> 22_08_2012 - initial version - J.L. Desmarais
        !
        !>@param comm_2d
        !> MPI Cartesian communicator betwwen the tiles
        !
        !>@param proc_x_choice
        !> type of procedure for the boundary layers in the 
        !> x-direction (E and W)
        !
        !>@param proc_y_choice
        !> type of procedure for the boundary layers in the 
        !> y-direction (N and S)
        !
        !> @param exchange_id
        !> cardinal direction in which the data are exchanged
        !> if the procedure for the computation of the boundary
        !> is both computing and exchanging grid points layers
        !--------------------------------------------------------------
        subroutine ini_bc_procedures(
     $       comm_2d, proc_x_choice, proc_y_choice, exchange_id)
        
          implicit none

          integer              , intent(in)  :: comm_2d
          integer              , intent(out) :: proc_x_choice
          integer              , intent(out) :: proc_y_choice
          integer, dimension(2), intent(out) :: exchange_id


          integer               :: proc_rank
          integer               :: dims_nb
          integer, dimension(2) :: cart_coord
          logical               :: tile_is_on_the_border
          logical               :: periodic_exchange
          integer               :: ierror


          !< get the rank of the processor computing the tile
          call MPI_COMM_RANK(comm_2d, proc_rank, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'ini_bc_procedures: MPI_COMM_RANK failed'
          end if


          !< get the cartesian coordinates of the processor
          !> computing the tile
          dims_nb=2
          call MPI_CART_COORDS(comm_2d, proc_rank,
     $         dims_nb, cart_coord, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             stop 'ini_bc_procedures:MPI_CART_COORDS failed'
          end if


          !< get the proc_x_choice and the exchange_id
          tile_is_on_the_border = 
     $         (cart_coord(1).eq.0).or.(cart_coord(1).eq.(npx-1))

          periodic_exchange =
     $         (bc_choice.eq.periodic_xy_choice).and.(npx.gt.1)

          if(tile_is_on_the_border.and.(.not.periodic_exchange)) then
             if(npx.eq.1) then
                proc_x_choice = only_compute_proc
             else
                proc_x_choice = compute_and_exchange_proc
                if(cart_coord(1).eq.0) then
                   exchange_id(x_direction)=E
                else
                   exchange_id(x_direction)=W
                end if
             end if

          else
             proc_x_choice = only_exchange_proc
          end if
          

          !< get the proc_y_choice and the exchange_id
          tile_is_on_the_border = 
     $         (cart_coord(2).eq.0).or.(cart_coord(2).eq.(npy-1))

          periodic_exchange =
     $         (bc_choice.eq.periodic_xy_choice).and.(npy.gt.1)

          if(tile_is_on_the_border.and.(.not.periodic_exchange)) then
             if(npy.eq.1) then
                proc_y_choice = only_compute_proc
             else
                proc_y_choice = compute_and_exchange_proc
                if(cart_coord(2).eq.0) then
                   exchange_id(y_direction)=N
                else
                   exchange_id(y_direction)=S
                end if
             end if

          else
             proc_y_choice = only_exchange_proc
          end if

        end subroutine ini_bc_procedures

      end module mpi_mg_ini_bc_proc
