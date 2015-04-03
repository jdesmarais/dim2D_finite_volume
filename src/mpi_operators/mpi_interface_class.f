      !> @file
      !> class encapsulating procedures to identify which processors
      !> are computing the neighboring tiles and to exchange nodes
      !> between the tiles using non-blocking communications
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating procedures to identify which processors
      !> are computing the neighboring tiles and to exchange nodes
      !> between the tiles using non-blocking communications
      !
      !> @date
      ! 03_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_interface_class

        use mpi

        use mpi_mg_bc_ext_class, only :
     $       mpi_mg_bc_ext

        use mpi_process_class, only :
     $       mpi_process

        use mpi_requests_module, only :
     $       create_requests_for_one_direction,
     $       only_exchange_twice

        use parameters_bf_layer, only :
     $       no_overlap,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SW_corner_type,
     $       SE_corner_type

        use parameters_constant, only :
     $       N,S,E,W,
     $       only_compute_proc,
     $       compute_and_exchange_proc,
     $       only_exchange_proc

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       npx,npy

        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none


        private
        public :: mpi_interface


        !> @class mpi_interface
        !> class encapsulating procedures to identify which processors
        !> are computing the neighboring tiles and to exchange nodes
        !> between the tiles using non-blocking communications
        !>
        !> @param nb_mpi_requests_x
        !> number of mpi requests emitted along the x-direction
        !>
        !> @param mpi_requests_x
        !> ID of the mpi requests emitted along the x-direction
        !>
        !> @param nb_mpi_requests_y
        !> number of mpi requests emitted along the y-direction
        !>
        !> @param mpi_requests_y
        !> ID of the mpi requests emitted along the y-direction
        !>
        !> @param ini
        !> initialize the number of requests for the x- and y-
        !> directions
        !
        !> @param ini_for_timeInt
        !> initialize the bc_sections to know which grid-points 
        !> of the boundary sections should be computed by the tile
        !> as well as the time integration borders
        !
        !> @param MPI_ISENDRECV_XDIR
        !> create the mpi requests for the exchanges along the
        !> x-direction (with the neighboring tiles E and W)
        !
        !> @param MPI_ISENDRECV_YDIR
        !> create the mpi requests for the exchanges along the
        !> y-direction (with the neighboring tiles N and S)
        !
        !> @param MPI_WAITALL_XDIR
        !> wait for the mpi requests to complete for the exchanges
        !> along the x-direction (with the neighboring tiles E and W)
        !
        !> @param MPI_WAITALL_YDIR
        !> wait for the mpi requests to complete for the exchanges
        !> along the y-direction (with the neighboring tiles N and S)
        !---------------------------------------------------------------
        type, extends(mpi_mg_bc_ext) :: mpi_interface

          integer               :: nb_mpi_requests_x
          integer, dimension(2) :: mpi_requests_x
          integer               :: nb_mpi_requests_y
          integer, dimension(2) :: mpi_requests_y

          contains

          procedure, pass :: ini
          procedure, pass :: ini_for_timeInt
          procedure, pass :: MPI_ISENDRECV_XDIR
          procedure, pass :: MPI_ISENDRECV_YDIR
          procedure, pass :: MPI_WAITALL_XDIR
          procedure, pass :: MPI_WAITALL_YDIR

        end type mpi_interface

        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the ID of the processors of the
        !> neighboring tiles as well as the MPI derived types used
        !> to exchange data and the number of MPI requests needed
        !> to exchange data in the x- and y- directions
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !--------------------------------------------------------------
        subroutine ini(this,comm_2d)

          implicit none

          class(mpi_interface), intent(inout) :: this
          integer             , intent(in)    :: comm_2d

          ! initialize the parent objects
          !============================================================
          call this%mpi_mg_bc_ext%ini(comm_2d)


          ! initialize the number of requests for the x- and
          ! y- directions
          !============================================================
          ! only_compute_proc :
          !      no information is exchanged, so there is no need
          !      to wait for exchanges to complete
          !                           
          ! only_exchange_proc :
          !      if the processors for the exchanges ((N and S) or (E and W))
          !      are the same, the data are combined into one data package
          !      and 2 MPI requests are needed: one to send, one to receive
          !      if the processors of the exchanges ((N and S) or (E and W))
          !      are not the same, 4 MPI requests are needed. However,
          !      the subroutine used to exchange the data has already a
          !      method built inside to wait for the requests to complete
          !      so we set nb_mpi_requests = 0 as we do not need to handle
          !      the waiting process
          !
          ! compute_and_exchange_proc :
          !      one of the two cardinal coordinates in the direction
          !      is exchanging and the other one is computed so there
          !      are two requests: one to send, one to receive
          !------------------------------------------------------------
          if(this%proc_x_choice.eq.compute_and_exchange_proc) then
             this%nb_mpi_requests_x = 2
          else
             if(this%proc_x_choice.eq.only_exchange_proc) then
                if(this%com_rank(E).eq.this%com_rank(W)) then
                   this%nb_mpi_requests_x = 2
                else
                   this%nb_mpi_requests_x = 0
                end if
             else
                this%nb_mpi_requests_x = 0
             end if
          end if

          if(this%proc_y_choice.eq.compute_and_exchange_proc) then
             this%nb_mpi_requests_y = 2
          else
             if(this%proc_y_choice.eq.only_exchange_proc) then
                if(this%com_rank(N).eq.this%com_rank(S)) then
                   this%nb_mpi_requests_y = 2
                else
                   this%nb_mpi_requests_y = 0
                end if
             else
                this%nb_mpi_requests_y = 0
             end if
          end if

        end subroutine ini
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialization of the bc_sections for the boundary sections
        !> whose nodes are computed and not exchanged as well as the
        !> bounds for the integration
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> MPI 2D cartesian communicator
        !
        !> @param x_borders
        !> x-borders for the time integration of the tile
        !
        !> @param y_borders
        !> y-borders for the time integration of the tile
        !
        !> @param bc_sections_x
        !> boundary sections along the x-direction
        !
        !> @param bc_sections_y
        !> boundary sections along the y-direction
        !--------------------------------------------------------------
        subroutine ini_for_timeInt(
     $     this,
     $     comm2d,
     $     x_borders,
     $     y_borders,
     $     bc_sections_x,
     $     bc_sections_y)

          implicit none

          class(mpi_interface)                       , intent(inout) :: this
          integer                                    , intent(in)    :: comm2d
          integer(ikind), dimension(2)               , intent(out)   :: x_borders
          integer(ikind), dimension(2)               , intent(out)   :: y_borders
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: bc_sections_x
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: bc_sections_y

          type(mpi_process) :: mpi_process_used

          integer               :: rank
          integer, dimension(2) :: cart_coords
          integer               :: ierror

          integer(ikind) :: edge_i_min
          integer(ikind) :: edge_i_max

          integer :: nb_bc_sections_y
          logical :: left_corner
          logical :: right_corner
          integer :: i


          !============================================================
          ! a boundary section is identified by 5 integer:
          ! [ type of bc_sections, i_min, j_min, extent or corner_overlap, overlap]
          !
          !   - type of bc_sections: N_edge_type   , S_edge_type,
          !                          E_edge_type   , W_edge_type,
          !                          NE_corner_type, NW_corner_type,
          !                          SE_corner_type, SW_corner_type
          !
          !   - i_min : minimum index along the x-direction
          !   - j_min : minimum index along the y-direction
          !   - extent: for an edge-type 
          !   - overlap: whether all grid-points are computed
          !============================================================

          ! determine the bc_sections along the x-direction
          !------------------------------------------------------------
          if(this%proc_x_choice.eq.only_compute_proc) then

             x_borders = [1,nx]

             allocate(bc_sections_x(5,2))

             bc_sections_x(:,1) = [W_edge_type, 1           , bc_size+1, ny-bc_size, no_overlap]
             bc_sections_x(:,2) = [E_edge_type, nx-bc_size+1, bc_size+1, ny-bc_size, no_overlap]

          else

             if(this%proc_x_choice.eq.compute_and_exchange_proc) then

                allocate(bc_sections_x(5,1))

                if(this%exchange_id(1).eq.E) then
                   x_borders = [1,nx-bc_size]
                   bc_sections_x(:,1) = [W_edge_type, 1           , bc_size+1, ny-bc_size, no_overlap]
                else
                   x_borders = [bc_size+1,nx]
                   bc_sections_x(:,1) = [E_edge_type, nx-bc_size+1, bc_size+1, ny-bc_size, no_overlap]
                end if

             else
                x_borders = [bc_size+1,nx-bc_size]
             end if

          end if

          ! determine the bc_sections along the y-direction
          !------------------------------------------------------------
          if(this%proc_y_choice.ne.only_exchange_proc) then

             ! determine the number of bc_sections and whether corners
             ! are computed
             !............................................................
             call MPI_COMM_RANK(comm2d,rank,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                print '(''MPI_COMM_RANK failed'')'
                call mpi_process_used%finalize_mpi()
             end if
             
             call MPI_CART_COORDS(comm2d, rank, 2, cart_coords, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                print '(''MPI_CART_COORDS failed'')'
                call mpi_process_used%finalize_mpi()
             end if
             
             edge_i_min = 1
             edge_i_max = nx
             
             
             nb_bc_sections_y = 1 !for the edge
             
             if(cart_coords(1).eq.0) then
                
                nb_bc_sections_y = nb_bc_sections_y+1 !for the left corner
                left_corner      = .true.
                edge_i_min       = bc_size+1
                
             end if

             if(cart_coords(1).eq.(npx-1)) then

                nb_bc_sections_y = nb_bc_sections_y+1 !for the right corner
                right_corner     = .true.
                edge_i_max       = nx-bc_size

             end if

             
             i=1

             ! determine the bc_sections for an only_compute case
             !............................................................
             if(this%proc_y_choice.eq.only_compute_proc) then

                y_borders = [1,ny]

                allocate(bc_sections_y(5,2*nb_bc_sections_y))

                ! south layer
                call add_south_bc_sections(
     $               left_corner,
     $               right_corner,
     $               edge_i_min,
     $               edge_i_max,
     $               i,
     $               bc_sections_y)

                ! north layer
                call add_north_bc_sections(
     $               left_corner,
     $               right_corner,
     $               edge_i_min,
     $               edge_i_max,
     $               i,
     $               bc_sections_y)


             ! determine the bc_sections for an compute_and_exchange case
             !............................................................
             else

                allocate(bc_sections_y(5,nb_bc_sections_y))

                if(this%exchange_id(2).eq.N) then

                   y_borders = [1,ny-bc_size]
                   
                   ! south layer
                   call add_south_bc_sections(
     $                  left_corner,
     $                  right_corner,
     $                  edge_i_min,
     $                  edge_i_max,
     $                  i,
     $                  bc_sections_y)

                else

                   y_borders = [bc_size+1,ny]

                   ! north layer
                   call add_north_bc_sections(
     $                  left_corner,
     $                  right_corner,
     $                  edge_i_min,
     $                  edge_i_max,
     $                  i,
     $                  bc_sections_y)

                end if

             end if

          else

             y_borders = [bc_size+1,ny-bc_size]

          end if

        end subroutine ini_for_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the south boundary sections to the boundary sections
        !> of the tile
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param left_corner
        !> add the left corner or not
        !
        !> @param right_corner
        !> add the right corner or not
        !
        !> @param edge_i_min
        !> minimum index or the edge along the x-direction
        !
        !> @param edge_i_max
        !> maximum index or the edge along the x-direction
        !
        !> @param i
        !> index where the new bc_section should be added
        !
        !> @param bc_sections_y
        !> boundary sections along the y-direction
        !--------------------------------------------------------------
        subroutine add_south_bc_sections(
     $     left_corner,
     $     right_corner,
     $     edge_i_min,
     $     edge_i_max,
     $     i,
     $     bc_sections_y)

          implicit none

          logical                       , intent(in)    :: left_corner
          logical                       , intent(in)    :: right_corner
          integer(ikind)                , intent(in)    :: edge_i_min
          integer(ikind)                , intent(in)    :: edge_i_max
          integer                       , intent(inout) :: i
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections_y
          

          if(left_corner) then
             bc_sections_y(:,i) = [SW_corner_type, 1, 1, no_overlap, no_overlap]
             i = i+1
          end if

          bc_sections_y(:,i) = [S_edge_type, edge_i_min, 1, edge_i_max, no_overlap]
          i = i+1
          
          if(right_corner) then
             bc_sections_y(:,i) = [SE_corner_type, nx-bc_size+1, 1, no_overlap, no_overlap]
             i = i+1
          end if

        end subroutine add_south_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the north boundary sections to the boundary sections
        !> of the tile
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param left_corner
        !> add the left corner or not
        !
        !> @param right_corner
        !> add the right corner or not
        !
        !> @param edge_i_min
        !> minimum index or the edge along the x-direction
        !
        !> @param edge_i_max
        !> maximum index or the edge along the x-direction
        !
        !> @param i
        !> index where the new bc_section should be added
        !
        !> @param bc_sections_y
        !> boundary sections along the y-direction
        !--------------------------------------------------------------
        subroutine add_north_bc_sections(
     $     left_corner,
     $     right_corner,
     $     edge_i_min,
     $     edge_i_max,
     $     i,
     $     bc_sections_y)

          implicit none

          logical                       , intent(in)    :: left_corner
          logical                       , intent(in)    :: right_corner
          integer(ikind)                , intent(in)    :: edge_i_min
          integer(ikind)                , intent(in)    :: edge_i_max
          integer                       , intent(inout) :: i
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections_y
          

          if(left_corner) then
             bc_sections_y(:,i) = [NW_corner_type, 1, ny-bc_size+1, no_overlap, no_overlap]
             i = i+1
          end if

          bc_sections_y(:,i) = [N_edge_type, edge_i_min, ny-bc_size+1, edge_i_max, no_overlap]
          i = i+1
          
          if(right_corner) then
             bc_sections_y(:,i) = [NE_corner_type, nx-bc_size+1, ny-bc_size+1, no_overlap, no_overlap]
             i = i+1
          end if

        end subroutine add_north_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create requests to send and receive data in the x-direction
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !
        !> @param nodes
        !> grid-points of the tile
        !--------------------------------------------------------------
        subroutine MPI_ISENDRECV_XDIR(
     $     this,
     $     comm_2d,
     $     nodes)

          implicit none

          class(mpi_interface)            , intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          call MPI_ISENDRECV_DIR(
     $         this,
     $         comm_2d,
     $         E,W,
     $         this%exchange_id(1),
     $         this%proc_x_choice,
     $         this%mpi_requests_x,
     $         nodes)

        end subroutine MPI_ISENDRECV_XDIR


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create requests to send and receive data in the y-direction
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !
        !> @param nodes
        !> grid-points of the tile
        !--------------------------------------------------------------
        subroutine MPI_ISENDRECV_YDIR(
     $     this,
     $     comm_2d,
     $     nodes)

          implicit none

          class(mpi_interface)            , intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          call MPI_ISENDRECV_DIR(
     $         this,
     $         comm_2d,
     $         N,S,
     $         this%exchange_id(2),
     $         this%proc_y_choice,
     $         this%mpi_requests_y,
     $         nodes)

        end subroutine MPI_ISENDRECV_YDIR


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create requests to send and received data
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !
        !> @param card1
        !> first cardinal coordinate identifying the exchange in
        !> the direction (E for x-dir, N for y-dir)
        !
        !> @param card2
        !> second cardinal coordinate identifying the exchange in
        !> the direction (W for x-dir, S for y-dir)
        !
        !> @param card_exchange
        !> cardinal coordinate identifying the exchange if one
        !> is computed and one is exchanged (exchange_id(1) for x-dir,
        !> exchange_id(2) for y-dir)
        !
        !> @param proc_n_choice
        !> identification of the process in the direction:
        !> proc_x_choice for x-dir
        !> proc_y_choice for y-dir
        !
        !> @param mpi_requests
        !> mpi requests that should be tracked for waiting when the
        !> exchange is overlap with computations
        !
        !> @param nodes
        !> grid-points of the tile
        !--------------------------------------------------------------
        subroutine MPI_ISENDRECV_DIR(
     $     this,
     $     comm_2d,
     $     card1, card2,
     $     card_exchange,
     $     proc_n_choice,
     $     mpi_requests,
     $     nodes)

          implicit none

          class(mpi_interface)            , intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: card1
          integer                         , intent(in)    :: card2
          integer                         , intent(in)    :: card_exchange
          integer                         , intent(in)    :: proc_n_choice
          integer    , dimension(2)       , intent(out)   :: mpi_requests
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          
          

          type(mpi_process) :: mpi_process_used

          integer :: ierror
          integer :: rank


          ! to prevent non-initialization warnings
          ! if the mi_requests are not initialized
          ! if it is not needed
          mpi_requests(1) = 1


          ! get the rank of the processor
          call MPI_COMM_RANK(comm_2d,rank,ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print '(''MPI_ISENDRECV_DIR'')'
             print '(''MPI_COMM_RANK failed'')'
             call mpi_process_used%finalize_mpi()
             stop
          end if


          ! if we are only exchanging data
          ! we either need 4 or 2 requests
          if(proc_n_choice.eq.only_exchange_proc) then

             ! 2 requests if we are exchanging with the same processor
             if(this%com_rank(card1).eq.this%com_rank(card2)) then
                
                mpi_requests = create_requests_for_one_direction(
     $               this, comm_2d, rank, nodes, card1)

             ! 4 requests if we are exchanging with two different
             ! processors
             ! the waiting process is built in the subroutine
             ! only_exchange_twice such that we do not need to
             ! track the requests for waiting
             else

                call only_exchange_twice(
     $               this, comm_2d, rank, nodes, npx*npy, [card1,card2])

             end if

          else

             ! if we are exchanging data only once, we need
             ! 2 MPI requests
             if(proc_n_choice.eq.compute_and_exchange_proc) then

                mpi_requests = create_requests_for_one_direction(
     $               this, comm_2d, rank, nodes, card_exchange)

             end if

          end if

        end subroutine MPI_ISENDRECV_DIR


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> wait for the requests in the x-direction to complete
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !--------------------------------------------------------------
        subroutine MPI_WAITALL_XDIR(this)

          implicit none

          class(mpi_interface), intent(inout) :: this


          type(mpi_process) :: mpi_process_used

          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer :: ierror


          if(this%nb_mpi_requests_x.ne.0) then

             call MPI_WAITALL(2, this%mpi_requests_x, status, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_process_used%finalize_mpi()
                print *, 'mpi_interface_class'
                print *, 'MPI_WAITALL_XDIR'
                stop 'MPI_WAITALL failed'
             end if
             
          end if

        end subroutine MPI_WAITALL_XDIR


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> wait for the requests in the y-direction to complete
        !
        !> @date
        !> 03_04_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> interface for exchanging the boundary sections with the
        !> neighboring tiles
        !
        !> @param comm2d
        !> cartesian 2D MPI communicator
        !--------------------------------------------------------------
        subroutine MPI_WAITALL_YDIR(this)

          implicit none

          class(mpi_interface), intent(inout) :: this


          type(mpi_process) :: mpi_process_used

          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer :: ierror


          if(this%nb_mpi_requests_y.ne.0) then          

             call MPI_WAITALL(2, this%mpi_requests_y, status, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_process_used%finalize_mpi()
                print *, 'mpi_interface_class'
                print *, 'MPI_WAITALL_YDIR'
                stop 'MPI_WAITALL failed'
             end if
             
          end if

        end subroutine MPI_WAITALL_YDIR

      end module mpi_interface_class
