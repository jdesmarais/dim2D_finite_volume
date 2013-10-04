      !> @file
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using wall xy boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using wall xy boundary conditions
      !> only compute and exchange subroutines are needed for the
      !> wall xy boundary conditions
      !
      !> @date
      ! 25_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_xy_par_module

        use cg_operators_class  , only : cg_operators
        use dim2d_eq_class      , only : dim2d_eq
        use field_par_class     , only : field_par
        use mpi                 
        use mpi_mg_bc_class     , only : mpi_mg_bc
        use mpi_process_class   , only : mpi_process
        use mpi_requests_module , only : create_requests_for_one_direction,
     $                                   only_exchange_twice
        use mpi_tag_module      , only : compute_mpi_tag
        use parameters_constant , only : x_direction, y_direction,
     $                                   N,S,E,W
        use parameters_input    , only : nx,ny,ne,npx,npy
        use parameters_kind     , only : ikind, rkind
        use wall_xy_module      , only : wall_prefactor,
     $                                   compute_wall_flux_x,
     $                                   compute_wall_flux_y

        implicit none

        private
        public :: only_compute_along_x,
     $            only_compute_along_y,
     $            only_exchange,
     $            compute_and_exchange_along_x,
     $            compute_and_exchange_along_y,
     $            fluxes_only_compute_along_x,
     $            fluxes_only_compute_along_y,
     $            fluxes_compute_and_exchange_along_x,
     $            fluxes_compute_and_exchange_along_y

        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> x-direction (E and W) using wall b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
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
        subroutine only_compute_along_x(nodes, bc_size, p_model)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer                         , intent(in)    :: bc_size
          type(dim2d_eq)                  , intent(in)    :: p_model

          
          !< prefactor : equal to -1 or +1 depending on the variable
          !>               type: vector_x/vector_y or not
          integer, dimension(ne) :: prefactor
          integer(ikind)         :: i,j
          integer                :: k


          !< compute the prefactor
          prefactor = wall_prefactor(p_model)


          !< compute the wall b.c. in E and W boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                  prefactor(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                  prefactor(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_x



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> y-direction (N and S) using wall b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
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

          
          !< prefactor : equal to -1 or +1 depending on the variable
          !>             type: vector_x/vector_y or not
          integer, dimension(ne) :: prefactor
          integer(ikind)         :: i,j
          integer                :: k


          !< compute the prefactor
          prefactor = wall_prefactor(p_model)


          !< compute the N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  prefactor(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  prefactor(k)*nodes(i,ny-bc_size-j+1,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing one boundary layers along the
        !> x-direction (E or W) using wall b.c. and exchanging
        !> the other one
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
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
          !< card_pt     : table corresponding to the two cardinal
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
          integer, dimension(ne)                :: prefactor


          !< create two requests in one direction: one for sending
          !> and the other one for receiving data from the same
          !> processor. This processor is identified by the cardinal
          !> direction 'card_pt'
          !DEC$ FORCEINLINE RECURSIVE
          mpi_requests = create_requests_for_one_direction(
     $         this, f_used, nodes, card_pt)

          !< overlap some communications with computations

          !< compute the wall prefactor
          prefactor = wall_prefactor(p_model)

          select case(card_pt)

            !< if we send along E, we compute W
            case(E)
                do k=1,ne
                   do j=1+bc_size, ny-bc_size
                      !DEC$ IVDEP
                      do i=1,bc_size
                   
                         nodes(i,j,k) = 
     $                        prefactor(k)*nodes(2*bc_size+1-i,j,k)
                   
                      end do
                   end do
                end do

            case(W)
               do k=1,ne
                   do j=1+bc_size, ny-bc_size
                      !DEC$ IVDEP
                      do i=1,bc_size
                   
                         nodes(nx-bc_size+i,j,k) = 
     $                        prefactor(k)*nodes(nx-bc_size-i+1,j,k)
                   
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
        !> y-direction (N and S) using wall b.c.
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
          !< status      : table identifying the status of the mpi
          !>               requests for sending and receiving data
          !-------------------------------------------------------
          type(mpi_process)                     :: mpi_op
          integer, dimension(2)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,2) :: status
          integer                               :: ierror,k
          integer(ikind)                        :: i,j
          integer, dimension(ne)                :: prefactor


          !< create two requests in one direction: one for sending
          !> and the other one for receiving data from the same
          !> processor. This processor is identified by the cardinal
          !> direction 'card_pt'
          !DEC$ FORCEINLINE RECURSIVE
          mpi_requests = create_requests_for_one_direction(
     $         this, f_used, nodes, card_pt)


          !< overlap some communications with computations

          !< compute the wall prefactor
          prefactor = wall_prefactor(p_model)

          select case(card_pt)

            !< if we send along N, we compute S
            case(N)
                do k=1, ne
                   do j=1, bc_size
                      !DEC$ IVDEP
                      do i=1, nx
                   
                         nodes(i,j,k) = 
     $                        prefactor(k)*nodes(i,2*bc_size+1-j,k)
                         
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
     $                        prefactor(k)*nodes(i,ny-bc_size-j+1,k)
                         
                      end do
                   end do
                end do

            case default
               call mpi_op%finalize_mpi()
               print '(''wall compute_and_exchange_along_y:'')'
               stop 'card_pt not recognized'
          end select


          !< wait for all requests to be finished
          call MPI_WAITALL(2, mpi_requests, status, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print *, 'wall_xy_par_model'
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
          


          !< card_pt     : table corresponding to the two cardinal
          !>               pts of the direction asked by the user
          !
          !< nb_procs    : total number of processors in the 
          !>               communicator
          !-------------------------------------------------------
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
     $         this, f_used, nodes, nb_procs, card_pt)

        end subroutine only_exchange

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> x-direction (E and W) using wall b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object containing the variables of the governing euqations
        !
        !>@param s_op
        !> spatial discretization operator
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !--------------------------------------------------------------
        subroutine fluxes_only_compute_along_x(f_used,s_op,flux_x)

          implicit none

          type(field_par)                   , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s_op
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          integer        :: k, bc_size
          integer(ikind), dimension(2) :: id


          !< get the size of the boundary layer
          bc_size = s_op%get_bc_size()

          !< provide the x-indices modified
          id(1)=bc_size+1
          id(2)=nx+1-bc_size

          !< modify the fluxes along the x-direction
          !> at the E and W borders
          !> W border: i= bc_size+1
          !> E border: i= nx-bc_size+1
          do k=1,2
             !DEC$ FORCEINLINE RECURSIVE
             call compute_wall_flux_x(f_used,s_op,id(k),flux_x)
          end do

        end subroutine fluxes_only_compute_along_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary fluxes along the
        !> y-direction (N and S) using wall b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object containing the variables of the governing euqations
        !
        !>@param s_op
        !> spatial discretization operator
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !--------------------------------------------------------------
        subroutine fluxes_only_compute_along_y(f_used,s_op,flux_y)

          implicit none

          type(field_par)                   , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s_op
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer        :: k, bc_size
          integer(ikind), dimension(2) :: id

          !< get the size of the boundary layer
          bc_size = s_op%get_bc_size()

          !< provide the y-indices modified
          id(1)=bc_size+1
          id(2)=ny-bc_size+1

          !< modify the fluxes along the y-direction
          !> at the N and S borders
          !> S border: j= bc_size+1
          !> N border: j= ny-bc_size+1
          do k=1,2
             !DEC FORCEINLINE RECURSIVE
             call compute_wall_flux_y(f_used,s_op,id(k),flux_y)
          end do          

        end subroutine fluxes_only_compute_along_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine modifying the fluxes along the x-direction
        !> for the E or the W layer
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !> @param f_used
        !> object encapsulating the main variables
        !
        !>@param s_op
        !> space discretization operators
        !
        !>@param card_pt
        !> cardinal point defining the direction in which data are
        !> exchanged with the neighbouring tile
        !
        !>@param flux_x
        !> flux along the x-direction
        !--------------------------------------------------------------
        subroutine fluxes_compute_and_exchange_along_x(
     $     f_used, s_op, card_pt, flux_x)

          implicit none

          type(field_par)                   , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s_op
          integer                           , intent(in)    :: card_pt
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          type(mpi_process) :: mpi_op
          integer           :: bc_size,i

          !get the size of the boundary layers
          bc_size = s_op%get_bc_size()

          !select the x-indice of the flux modified: either E or W
          select case(card_pt)
            case(E)
               i=bc_size+1
            case(W)
               i=nx-bc_size+1
            case default
               call mpi_op%finalize_mpi()
               print '(''fluxes_compute_and_exchange_along_x:'')'
               stop 'card_pt not recognized'
          end select

          !modify the fluxes
          call compute_wall_flux_x(f_used,s_op,i,flux_x)

        end subroutine fluxes_compute_and_exchange_along_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine modifying the fluxes along the y-direction
        !> for the N or the S layer
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !> @param f_used
        !> object encapsulating the main variables
        !
        !>@param s_op
        !> space discretization operators
        !
        !>@param card_pt
        !> cardinal point defining the direction in which data are
        !> exchanged with the neighbouring tile
        !
        !>@param flux_y
        !> flux along the y-direction
        !--------------------------------------------------------------
        subroutine fluxes_compute_and_exchange_along_y(
     $     f_used, s_op, card_pt, flux_y)

          implicit none

          type(field_par)                   , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s_op
          integer                           , intent(in)    :: card_pt
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          type(mpi_process) :: mpi_op
          integer           :: bc_size,j

          !get the size of the boundary layers
          bc_size = s_op%get_bc_size()

          !select the x-indice of the flux modified: either E or W
          select case(card_pt)
            case(N)
               j=bc_size+1
            case(S)
               j=ny-bc_size+1
            case default
               call mpi_op%finalize_mpi()
               print '(''fluxes_compute_and_exchange_along_x:'')'
               stop 'card_pt not recognized'
          end select

          !modify the fluxes
          call compute_wall_flux_y(f_used,s_op,j,flux_y)

        end subroutine fluxes_compute_and_exchange_along_y

      end module wall_xy_par_module
