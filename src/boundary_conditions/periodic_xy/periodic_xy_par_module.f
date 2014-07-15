      !> @file
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using periodic xy boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using periodic xy boundary conditions
      !> only compute and exchange subroutines are needed for the
      !> periodic xy boundary conditions
      !
      !> @date
      ! 27_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module periodic_xy_par_module

        use mpi
        use mpi_mg_bc_class    , only : mpi_mg_bc
        use mpi_process_class  , only : mpi_process
        use mpi_requests_module, only : only_exchange_twice
        use mpi_tag_module     , only : compute_mpi_tag
        use parameters_constant, only : N,S,E,W,x_direction
        use parameters_input   , only : nx,ny,ne,npx,npy,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none
        
        private
        public :: only_compute_along_x,
     $            only_compute_along_y,
     $            only_exchange

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the periodic boundary conditions
        !> along the x-direction
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> main variables on the 2d computational field
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine only_compute_along_x(nodes)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          integer(ikind) :: i,j
          integer        :: k,period_x


          period_x = nx-2*bc_size


          !<compute the east and west boundary layers
          !>without the north and south corners
          do k=1, ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1, bc_size

                   nodes(i,j,k)=nodes(i+period_x,j,k)
                   nodes(i+period_x+bc_size,j,k)=nodes(i+bc_size,j,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the periodic boundary conditions
        !> along the y-direction
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> main variables on the 2d computational field
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine only_compute_along_y(nodes)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          integer(ikind) :: i,j
          integer        :: k,period_y


          period_y = ny-2*bc_size


          !<compute the south and north layers
          !>with the east and west corners
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx

                   nodes(i,j,k)=nodes(i,j+period_y,k)
                   nodes(i,j+period_y+bc_size,k)=nodes(i,j+bc_size,k)

                end do
             end do
          end do

        end subroutine only_compute_along_y
      

c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine computing the boundary layers along the
c$$$        !> x-direction (E and W) using periodic b.c.
c$$$        !
c$$$        !> @date
c$$$        !> 26_08_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> object encapsulating the rank of the processors computing
c$$$        !> the neighbouring tiles as well as the MPI derived types to
c$$$        !> identify the location of the data sent and received
c$$$        !
c$$$        !> @param f_used
c$$$        !> object encapsulating the main variables
c$$$        !
c$$$        !>@param nodes
c$$$        !> table containing the gridpoint data
c$$$        !
c$$$        !>@param bc_size
c$$$        !> size of the boundary layer
c$$$        !
c$$$        !>@param p_model
c$$$        !> physical model
c$$$        !
c$$$        !>@param card_pt
c$$$        !> cardinal point identifying the direction in which the data
c$$$        !> are sent
c$$$        !--------------------------------------------------------------
c$$$        subroutine compute_and_exchange_along_x(
c$$$     $     this, f_used, nodes, bc_size, card_pt)
c$$$
c$$$          implicit none
c$$$
c$$$          class(mpi_mg_bc)                , intent(in)    :: this
c$$$          class(field_par)                , intent(inout) :: f_used
c$$$          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
c$$$          integer                         , intent(in)    :: bc_size
c$$$          integer                         , intent(in)    :: card_pt
c$$$
c$$$
c$$$          !< mpi_op      : mpi process to finalize in case of error
c$$$          !
c$$$          !< mpi_requests: integer identifying the mpi requests
c$$$          !>               to send and receive information in the
c$$$          !>               direction asked by the user
c$$$          !
c$$$          !< cart_pt     : table corresponding to the two cardinal
c$$$          !>               pts of the direction asked by the user
c$$$          !
c$$$          !< nb_procs    : total number of processors in the 
c$$$          !>               communicator
c$$$          !
c$$$          !< status      : table identifying the status of the mpi
c$$$          !>               requests for sending and receiving data
c$$$          !-------------------------------------------------------
c$$$          type(mpi_process)                     :: mpi_op
c$$$          integer, dimension(2)                 :: mpi_requests
c$$$          integer                               :: nb_procs
c$$$          integer, dimension(MPI_STATUS_SIZE,2) :: status
c$$$          integer                               :: ierror,tag,k
c$$$          integer(ikind)                        :: i,j
c$$$          integer                               :: period_x
c$$$
c$$$
c$$$          nb_procs = npx*npy
c$$$
c$$$
c$$$          !< compute the tag identifying the sending MPI request
c$$$          tag = compute_mpi_tag(
c$$$     $         usr_rank, this%com_rank(card_pt),
c$$$     $         comm_2d, nb_procs)
c$$$   
c$$$
c$$$          !< create a send request
c$$$          call MPI_ISSEND(
c$$$     $         nodes, 1, this%com_send(card_pt),
c$$$     $         this%com_rank(card_pt), tag,
c$$$     $         comm_2d, mpi_requests(1),ierror)
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             call mpi_op%finalize_mpi()
c$$$             stop 'reflection_xy_par_module: MPI_ISSEND failed'
c$$$          end if
c$$$             
c$$$
c$$$          
c$$$          !< compute the tag identifying the receving MPI request
c$$$          tag = compute_mpi_tag(
c$$$     $         this%com_rank(card_pt), usr_rank,
c$$$     $         comm_2d, nb_procs)
c$$$
c$$$
c$$$          !< create a receive request
c$$$          call MPI_IRECV(
c$$$     $         nodes, 1, this%com_recv(card_pt),
c$$$     $         this%com_rank(card_pt), tag,
c$$$     $         comm_2d, mpi_requests(2),ierror) !<compute the tag
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             call mpi_op%finalize_mpi()
c$$$             stop 'reflection_xy_par_module: MPI_IRECV failed'
c$$$          end if
c$$$
c$$$
c$$$          !< overlap some communications with computations
c$$$
c$$$          !< compute the period_x
c$$$          period_x = nx-2*bc_size
c$$$
c$$$          select case(card_pt)
c$$$
c$$$            !< if we send along E, we compute W
c$$$            case(E)
c$$$                do k=1,ne
c$$$                   do j=1+bc_size, ny-bc_size
c$$$                      do i=1,bc_size
c$$$                   
c$$$                         nodes(i,j,k) = nodes(i+period_x,j,k)
c$$$                   
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$
c$$$            case(W)
c$$$               do k=1,ne
c$$$                   do j=1+bc_size, ny-bc_size
c$$$                      do i=1,bc_size
c$$$                   
c$$$                         nodes(i+period_x+bc_size,j,k)=
c$$$     $                        nodes(i+bc_size,j,k)
c$$$                   
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$
c$$$            case default
c$$$               call mpi_op%finalize_mpi()
c$$$               print '(''compute_and_exchange_along_x:'')'
c$$$               stop 'cart_pt not recognized'
c$$$          end select
c$$$
c$$$
c$$$          !< wait for all requests to be finished
c$$$          call MPI_WAITALL(2, mpi_requests, status, ierror)
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             print *, 'periodic_xy_par_module'
c$$$             print *, 'compute_and_exchange_along_x'
c$$$             stop 'MPI_WAITALL failed'
c$$$          end if
c$$$
c$$$        end subroutine compute_and_exchange_along_x
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine computing the boundary layers along the
c$$$        !> y-direction (N and S) using reflection b.c.
c$$$        !
c$$$        !> @date
c$$$        !> 26_08_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> object encapsulating the rank of the processors computing
c$$$        !> the neighbouring tiles as well as the MPI derived types to
c$$$        !> identify the location of the data sent and received
c$$$        !
c$$$        !> @param f_used
c$$$        !> object encapsulating the main variables
c$$$        !
c$$$        !>@param nodes
c$$$        !> table containing the gridpoint data
c$$$        !
c$$$        !>@param bc_size
c$$$        !> size of the boundary layer
c$$$        !
c$$$        !>@param p_model
c$$$        !> physical model
c$$$        !
c$$$        !>@param card_pt
c$$$        !> cardinal point identifying the direction in which the data
c$$$        !> are sent
c$$$        !--------------------------------------------------------------
c$$$        subroutine compute_and_exchange_along_y(
c$$$     $     this, f_used, nodes, bc_size, card_pt)
c$$$
c$$$          implicit none
c$$$
c$$$          class(mpi_mg_bc)                , intent(in)    :: this
c$$$          class(field_par)                , intent(inout) :: f_used
c$$$          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
c$$$          integer                         , intent(in)    :: bc_size
c$$$          integer                         , intent(in)    :: card_pt
c$$$
c$$$
c$$$          !< mpi_op      : mpi process to finalize in case of error
c$$$          !
c$$$          !< mpi_requests: integer identifying the mpi requests
c$$$          !>               to send and receive information in the
c$$$          !>               direction asked by the user
c$$$          !
c$$$          !< cart_pt     : table corresponding to the two cardinal
c$$$          !>               pts of the direction asked by the user
c$$$          !
c$$$          !< nb_procs    : total number of processors in the 
c$$$          !>               communicator
c$$$          !
c$$$          !< status      : table identifying the status of the mpi
c$$$          !>               requests for sending and receiving data
c$$$          !-------------------------------------------------------
c$$$          type(mpi_process)                     :: mpi_op
c$$$          integer, dimension(2)                 :: mpi_requests
c$$$          integer                               :: nb_procs
c$$$          integer, dimension(MPI_STATUS_SIZE,2) :: status
c$$$          integer                               :: ierror,tag,k
c$$$          integer(ikind)                        :: i,j
c$$$          integer                               :: period_y
c$$$
c$$$
c$$$          nb_procs = npx*npy
c$$$
c$$$
c$$$          !< compute the tag identifying the sending MPI request
c$$$          tag = compute_mpi_tag(
c$$$     $         usr_rank, this%com_rank(card_pt),
c$$$     $         comm_2d, nb_procs)
c$$$   
c$$$
c$$$          !< create a send request
c$$$          call MPI_ISSEND(
c$$$     $         nodes, 1, this%com_send(card_pt),
c$$$     $         this%com_rank(card_pt), tag,
c$$$     $         comm_2d, mpi_requests(1),ierror)
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             call mpi_op%finalize_mpi()
c$$$             stop 'reflection_xy_par_module: MPI_ISSEND failed'
c$$$          end if
c$$$
c$$$          
c$$$          !< compute the tag identifying the receving MPI request
c$$$          tag = compute_mpi_tag(
c$$$     $         this%com_rank(card_pt), usr_rank,
c$$$     $         comm_2d, nb_procs)
c$$$
c$$$
c$$$          !< create a receive request
c$$$          call MPI_IRECV(
c$$$     $         nodes, 1, this%com_recv(card_pt),
c$$$     $         this%com_rank(card_pt), tag,
c$$$     $         comm_2d, mpi_requests(2),ierror)
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             call mpi_op%finalize_mpi()
c$$$             stop 'reflection_xy_par_module: MPI_IRECV failed'
c$$$          end if
c$$$
c$$$
c$$$          !< overlap some communications with computations
c$$$
c$$$          !< compute the period_y
c$$$          period_y = ny-2*bc_size
c$$$
c$$$          select case(card_pt)
c$$$
c$$$            !< if we send along N, we compute S
c$$$            case(N)
c$$$                do k=1, ne
c$$$                   do j=1, bc_size
c$$$                      do i=1, nx
c$$$                   
c$$$                         nodes(i,j,k)=nodes(i,j+period_y,k)
c$$$                         
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$
c$$$            !< if we send along S, we compute N
c$$$            case(S)
c$$$               do k=1, ne
c$$$                   do j=1, bc_size
c$$$                      do i=1, nx
c$$$                   
c$$$                         nodes(i,j+period_y+bc_size,k)=
c$$$     $                        nodes(i,j+bc_size,k)
c$$$                         
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$
c$$$            case default
c$$$               call mpi_op%finalize_mpi()
c$$$               print '(''compute_and_exchange_along_y:'')'
c$$$               stop 'cart_pt not recognized'
c$$$          end select
c$$$
c$$$
c$$$          !< wait for all requests to be finished
c$$$          call MPI_WAITALL(2, mpi_requests, status, ierror)
c$$$          if(ierror.ne.MPI_SUCCESS) then
c$$$             print *, 'periodic_xy_par_module'
c$$$             print *, 'compute_and_exchange_along_x'
c$$$             stop 'MPI_WAITALL failed'
c$$$          end if
c$$$
c$$$        end subroutine compute_and_exchange_along_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers in one direction
        !> by only exchanging data with its neighbours
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the rank of the processors computing
        !> the neighbouring tiles as well as the MPI derived types to
        !> identify the locatio of the data sent and received
        !
        !> @param f_used
        !> object encapsulating the main variables
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param direction
        !> direction in which the data are sent (x or y axis)
        !--------------------------------------------------------------
        subroutine only_exchange(
     $     this, comm_2d, usr_rank, nodes, direction)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: direction
          

          !< cart_pt     : table corresponding to the two cardinal
          !>               pts of the direction asked by the user
          !
          !< nb_procs    : total number of processors in the 
          !>               communicator
          !-------------------------------------------------------
          integer, dimension(2) :: card_pt
          integer               :: nb_procs

          
          !< choose the direction in which data are send
          if(direction.eq.x_direction) then
             card_pt = [E,W]
          else
             card_pt = [N,S]
          end if               

          nb_procs = npx*npy

          
          !< if the neighbouring tiles to which the processor
          !> is sending and receiving information are computed
          !> by the same processor, the data are gathered in
          !> one MPI structure. So only one sending and one 
          !> receiving requests are required. Otherwise, we need
          !> to divide the sending and receiving requests into
          !> two : the cardinal directions are treated separately
          !-------------------------------------------------------

          !< two sending and two receiving requests are needed
          if(this%com_rank(card_pt(1)).ne.this%com_rank(card_pt(2))) then
             
             call only_exchange_twice(
     $            this, comm_2d, usr_rank, nodes, nb_procs, card_pt)

          !< only one sending and one receiving requests are needed
          else

             call only_exchange_once(
     $            this, comm_2d, usr_rank, nodes, nb_procs, card_pt(1))

          end if

        end subroutine only_exchange


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers in one direction
        !> by only exchanging data with its neighbours using four
        !> requests
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the rank of the processors computing
        !> the neighbouring tiles as well as the MPI derived types to
        !> identify the location of the data sent and received
        !
        !> @param f_used
        !> object encapsulating the main variables
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param nb_procs
        !> total number of processors
        !
        !>@param card_pt
        !> cardinal directions in which the data are sent (x or y axis)
        !--------------------------------------------------------------
        subroutine only_exchange_once(
     $     this, comm_2d, usr_rank, nodes, nb_procs, card_pt)
        
          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: nb_procs
          integer                         , intent(in)    :: card_pt


          !< mpi_op      : mpi process to finalize in case of error
          !
          !< mpi_requests: integer identifying the mpi requests
          !>               to send and receive information in the
          !>               direction asked by the user
          !
          !< status      : table identifying the status of the mpi
          !>               requests for sending and receiving data
          !-------------------------------------------------------
          type(mpi_process)                     :: mpi_op
          integer, dimension(2)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer                               :: ierror,tag


          !< compute the tag identifying the sending MPI request
          tag = compute_mpi_tag(
     $         usr_rank, this%com_rank(card_pt), nb_procs)

          
          !< create a send request
          call MPI_ISSEND(
     $         nodes, 1, this%com_send(card_pt),
     $         this%com_rank(card_pt), tag,
     $         comm_2d, mpi_requests(1),ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'periodic_xy_par_module: MPI_ISSEND failed'
          end if

          
          !< compute the tag identifying the receving MPI request
          tag = compute_mpi_tag(
     $         this%com_rank(card_pt), usr_rank, nb_procs)

          
          !< create a receive request
          call MPI_IRECV(
     $         nodes, 1, this%com_recv(card_pt),
     $         this%com_rank(card_pt), tag,
     $         comm_2d, mpi_requests(2),ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'periodic_xy_par_module: MPI_IRECV failed'
          end if
           
           
          !< wait for all requests to be finished
          call MPI_WAITALL(2, mpi_requests, status, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             print *, 'periodic_xy_par_model'
             print *, 'only_exchange'
             stop 'MPI_WAITALL failed'
          end if

        end subroutine only_exchange_once     

      end module periodic_xy_par_module
