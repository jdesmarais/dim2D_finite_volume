      !> @file
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using reflection xy boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using reflection xy boundary conditions
      !> only compute and exchange subroutines are needed for the
      !> reflection xy boundary conditions
      !
      !> @date
      ! 26_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module reflection_xy_par_module

        use dim2d_eq_class      , only : dim2d_eq
        use field_par_class     , only : field_par
        use mpi                 
        use mpi_mg_bc_class     , only : mpi_mg_bc
        use mpi_process_class   , only : mpi_process
        use mpi_requests_module , only : create_requests_for_one_direction,
     $                                   only_exchange_twice
        use mpi_tag_module      , only : compute_mpi_tag
        use parameters_constant , only : x_direction, y_direction,N,S,E,W
        use parameters_input    , only : nx,ny,ne,npx,npy
        use parameters_kind     , only : ikind, rkind
        use reflection_xy_module, only : reflection_x_prefactor,
     $                                   reflection_y_prefactor

        implicit none

        private
        public :: only_compute_along_x, only_compute_along_y,
     $            only_exchange,
     $            compute_and_exchange_along_x, compute_and_exchange_along_y

        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> x-direction (E and W) using reflection b.c.
        !
        !> @date
        !> 23_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param p_model
        !> physical model
        !--------------------------------------------------------------
        subroutine only_compute_along_x(nodes, bc_size, p_model)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          type(dim2d_eq)                  , intent(in)    :: p_model

          
          !< x_prefactor : equal to -1 or +1 depending on the variable
          !>               type: vector_x or not
          integer, dimension(ne) :: x_prefactor
          integer(ikind)         :: i,j
          integer                :: k


          !< compute the prefactor
          x_prefactor = reflection_x_prefactor(p_model)


          !< compute the reflection b.c. in E and W boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                  x_prefactor(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                  x_prefactor(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_x



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> y-direction (N and S) using reflection b.c.
        !
        !> @date
        !> 23_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param p_model
        !> physical model
        !--------------------------------------------------------------
        subroutine only_compute_along_y(nodes, bc_size, p_model)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          type(dim2d_eq)                  , intent(in)    :: p_model

          
          !< y_prefactor : equal to -1 or +1 depending on the variable
          !>               type: vector_y or not
          integer, dimension(ne) :: y_prefactor
          integer(ikind)         :: i,j
          integer                :: k


          !< compute the prefactor
          y_prefactor = reflection_y_prefactor(p_model)


          !< compute the N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  y_prefactor(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  y_prefactor(k)*nodes(i,ny-bc_size-j+1,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> x-direction (E and W) using reflection b.c.
        !
        !> @date
        !> 23_08_2013 - initial version - J.L. Desmarais
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
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param p_model
        !> physical model
        !
        !>@param card_pt
        !> cardinal point identifying the direction in which the data
        !> are sent
        !--------------------------------------------------------------
        subroutine compute_and_exchange_along_x(
     $     this, f_used, nodes, bc_size, p_model, card_pt)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          class(field_par)                , intent(inout) :: f_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          type(dim2d_eq)                  , intent(in)    :: p_model
          integer                         , intent(in)    :: card_pt


          !< mpi_op      : mpi process to finalize in case of error
          !
          !< mpi_requests: integer identifying the mpi requests
          !>               to send and receive information in the
          !>               direction asked by the user
          !
          !< cart_pt     : table corresponding to the two cardinal
          !>               pts of the direction asked by the user
          !
          !< nb_procs    : total number of processors in the 
          !>               communicator
          !
          !< status      : table identifying the status of the mpi
          !>               requests for sending and receiving data
          !-------------------------------------------------------
          type(mpi_process)                     :: mpi_op
          integer, dimension(2)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer                               :: ierror,k
          integer(ikind)                        :: i,j
          integer, dimension(ne)                :: x_prefactor


          !< create two requests in one direction: one for sending
          !> and the other one for receiving data from the same
          !> processor. This processor is identified by the cardinal
          !> direction 'cart_pt'
          !DEC$ FORCEINLINE RECURSIVE
          mpi_requests = create_requests_for_one_direction(
     $         this, f_used, nodes, card_pt)


          !< overlap some communications with computations

          !< compute the reflection prefactor
          x_prefactor = reflection_x_prefactor(p_model)

          select case(card_pt)

            !< if we send along E, we compute W
            case(E)
                do k=1,ne
                   do j=1+bc_size, ny-bc_size
                      !DEC$ IVDEP
                      do i=1,bc_size
                   
                         nodes(i,j,k) = 
     $                        x_prefactor(k)*nodes(2*bc_size+1-i,j,k)
                   
                      end do
                   end do
                end do

            case(W)
               do k=1,ne
                   do j=1+bc_size, ny-bc_size
                      !DEC$ IVDEP
                      do i=1,bc_size
                   
                         nodes(nx-bc_size+i,j,k) = 
     $                        x_prefactor(k)*nodes(nx-bc_size-i+1,j,k)
                   
                      end do
                   end do
                end do

            case default
               call mpi_op%finalize_mpi()
               print '(''compute_and_exchange_along_x:'')'
               stop 'card_pt not recognized'
          end select


          !< wait for all requests to be finished
          call MPI_WAITALL(2, mpi_requests, status, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print *, 'wall_xy_par_module'
             print *, 'compute_and_exchange_along_x'
             stop 'MPI_WAITALL failed'
          end if


        end subroutine compute_and_exchange_along_x
        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> y-direction (N and S) using reflection b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
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
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param p_model
        !> physical model
        !
        !>@param card_pt
        !> cardinal point identifying the direction in which the data
        !> are sent
        !--------------------------------------------------------------
        subroutine compute_and_exchange_along_y(
     $     this, f_used, nodes, bc_size, p_model, card_pt)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          class(field_par)                , intent(inout) :: f_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          type(dim2d_eq)                  , intent(in)    :: p_model
          integer                         , intent(in)    :: card_pt


          !< mpi_op      : mpi process to finalize in case of error
          !
          !< mpi_requests: integer identifying the mpi requests
          !>               to send and receive information in the
          !>               direction asked by the user
          !
          !< cart_pt     : table corresponding to the two cardinal
          !>               pts of the direction asked by the user
          !
          !< nb_procs    : total number of processors in the 
          !>               communicator
          !
          !< status      : table identifying the status of the mpi
          !>               requests for sending and receiving data
          !-------------------------------------------------------
          type(mpi_process)                     :: mpi_op
          integer, dimension(2)                 :: mpi_requests
          integer                               :: nb_procs
          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer                               :: ierror,tag,k
          integer(ikind)                        :: i,j
          integer, dimension(ne)                :: y_prefactor


          !< create two requests in one direction: one for sending
          !> and the other one for receiving data from the same
          !> processor. This processor is identified by the cardinal
          !> direction 'cart_pt'
          !DEC$ FORCEINLINE RECURSIVE
          mpi_requests = create_requests_for_one_direction(
     $       this, f_used, nodes, card_pt)

          !< overlap some communications with computations

          !< compute the reflection prefactor
          y_prefactor = reflection_y_prefactor(p_model)

          select case(card_pt)

            !< if we send along N, we compute S
            case(N)
                do k=1, ne
                   do j=1, bc_size
                      !DEC$ IVDEP
                      do i=1, nx
                   
                         nodes(i,j,k) = 
     $                        y_prefactor(k)*nodes(i,2*bc_size+1-j,k)
                         
                      end do
                   end do
                end do

            !< if we send along S, we compute N
            case(S)
               do k=1, ne
                   do j=1, bc_size
                      !DEC$ IVDEP
                      do i=1, nx
                   
                         nodes(i,ny-bc_size+j,k) = 
     $                        y_prefactor(k)*nodes(i,ny-bc_size-j+1,k)
                         
                      end do
                   end do
                end do

            case default
               call mpi_op%finalize_mpi()
               print '(''compute_and_exchange_along_y:'')'
               stop 'cart_pt not recognized'
          end select


          !< wait for all requests to be finished
          call MPI_WAITALL(2, mpi_requests, status, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print *, 'reflection_xy_par_model'
             print *, 'compute_and_exchange_along_x'
             stop 'MPI_WAITALL failed'
          end if

        end subroutine compute_and_exchange_along_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers in one direction
        !> by only exchanging data with its neighbours
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
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
        subroutine only_exchange(this, f_used, nodes, direction)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          class(field_par)                , intent(inout) :: f_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: direction
          

          integer, dimension(2)                 :: card_pt
          integer                               :: nb_procs
          
          !< choose the direction in which data are send
          if(direction.eq.x_direction) then
             card_pt = [E,W]
          else
             card_pt = [N,S]
          end if               

          nb_procs = npx*npy


          call only_exchange_twice(
     $            this, f_used, nodes, nb_procs, card_pt)

        end subroutine only_exchange

      end module reflection_xy_par_module
