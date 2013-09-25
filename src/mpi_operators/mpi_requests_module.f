      !> @file
      !> module encapsulating subroutines for the computation of 
      !> mpi requests
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of 
      !> mpi requests
      !
      !> @date
      ! 25_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_requests_module

        use field_par_class  , only : field_par
        use mpi
        use mpi_mg_bc_class  , only : mpi_mg_bc
        use mpi_process_class, only : mpi_process
        use mpi_tag_module   , only : compute_mpi_tag
        use parameters_input , only : nx,ny,ne,npx,npy
        use parameters_kind  , only : rkind

        implicit none

        private
        public :: create_requests_for_one_direction,
     $            only_exchange_twice
        
        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the mpi requests for non-blocking 
        !> communications: sending and receiving data in one direction
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object containing informations on the mpi derived types
        !> exchanged as well as the processor ids computing the
        !> neighbouring tiles
        !
        !>@param f_used
        !> object containing information about the communicator used
        !> between the tiles as well as the rank of the processor
        !> computing the current tile
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param cart_pt
        !> cardinal point identifying the direction in which data are
        !> sent
        !--------------------------------------------------------------
        function create_requests_for_one_direction(
     $       this, f_used, nodes, bc_size, card_pt)
     $       result(mpi_requests)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          class(field_par)                , intent(inout) :: f_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          integer                         , intent(in)    :: card_pt
          integer, dimension(2)                          :: mpi_requests

          type(mpi_process) :: mpi_op
          integer           :: nb_procs
          integer           :: ierror
          integer           :: tag


          !total number of processor
          nb_procs = npx*npy


          !< compute the tag identifying the sending MPI request
          tag = compute_mpi_tag(
     $         f_used%usr_rank, this%com_rank(card_pt), nb_procs)
   

          !< create a sending request
          call MPI_ISSEND(
     $         nodes, 1, this%com_send(card_pt),
     $         this%com_rank(card_pt), tag,
     $         f_used%comm_2d, mpi_requests(1),ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'mpi_requests_module: MPI_ISSEND failed'
          end if
             

          
          !< compute the tag identifying the receving MPI request
          tag = compute_mpi_tag(
     $         this%com_rank(card_pt), f_used%usr_rank, nb_procs)


          !< create a receiving request
          call MPI_IRECV(
     $         nodes, 1, this%com_recv(card_pt),
     $         this%com_rank(card_pt), tag,
     $         f_used%comm_2d, mpi_requests(2),ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'mpi_requests_module: MPI_IRECV failed'
          end if

        end function create_requests_for_one_direction


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
        subroutine only_exchange_twice(
     $     this, f_used, nodes, nb_procs, card_pt)
        
          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          class(field_par)                , intent(inout) :: f_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: nb_procs
          integer, dimension(2)           , intent(in)    :: card_pt


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
          integer, dimension(4)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,4) :: status
          integer                               :: ierror,k,tag


           !< create the MPI requests to send and receive in
           !> the direction asked by the user
           do k=1,2
           
              !< compute the tag identifying the sending MPI request
              tag = compute_mpi_tag(
     $             f_used%usr_rank, this%com_rank(card_pt(k)), nb_procs)
           
              !< create a send request
              call MPI_ISSEND(
     $             nodes, 1, this%com_send(card_pt(k)),
     $             this%com_rank(card_pt(k)), tag,
     $             f_used%comm_2d, mpi_requests(2*k-1),ierror)
              if(ierror.ne.MPI_SUCCESS) then
                 call mpi_op%finalize_mpi()
                 stop 'reflection_xy_par_module: MPI_ISSEND failed'
              end if
              
              !< compute the tag identifying the receving MPI request
              tag = compute_mpi_tag(
     $             this%com_rank(card_pt(k)), f_used%usr_rank, nb_procs)
           
              !< create a receive request
              call MPI_IRECV(
     $             nodes, 1, this%com_recv(card_pt(k)),
     $             this%com_rank(card_pt(k)), tag,
     $             f_used%comm_2d, mpi_requests(2*k),ierror)
              if(ierror.ne.MPI_SUCCESS) then
                 call mpi_op%finalize_mpi()
                 stop 'reflection_xy_par_module: MPI_IRECV failed'
              end if
           
           end do
           
           
           !< wait for all requests to be finished
           call MPI_WAITALL(4, mpi_requests, status, ierror)
           if(ierror.ne.MPI_SUCCESS) then
              call mpi_op%finalize_mpi()
              print *, 'reflection_xy_par_model'
              print *, 'only_exchange'
              stop 'MPI_WAITALL failed'
           end if

        end subroutine only_exchange_twice

      end module mpi_requests_module
