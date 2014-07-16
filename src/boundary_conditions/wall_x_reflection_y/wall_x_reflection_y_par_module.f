      !> @file
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using wall and reflection boundary conditions
      !
      !
      !                      wall b.c.
      !                          |
      !>                 //////////////////
      !>                 |               //
      !>                 |               //
      !>    reflection ->| computational //  <-  wall b.c.
      !>     b.c.        |    domain     //
      !>                 |               //
      !>                 |               //
      !                  //////////////////
      !                          |
      !                      wall b.c.
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of 
      !> boundary layers using wall and reflection boundary conditions
      !
      !> @date
      ! 27_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_x_reflection_y_par_module

        use mpi                 

        use mpi_mg_bc_class, only :
     $       mpi_mg_bc

        use mpi_process_class, only :
     $       mpi_process

        use mpi_requests_module, only :
     $       create_requests_for_one_direction

        use parameters_constant, only :
     $       x_direction, y_direction,
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,npx,npy,bc_size

        use parameters_kind, only :
     $       ikind, rkind

        use sd_operators_class, only :
     $       sd_operators

        use wall_xy_module, only :
     $       compute_wall_flux_x,
     $       compute_wall_flux_y


        implicit none


        private
        public :: only_compute_along_x,
     $            compute_and_exchange_along_x,
     $            fluxes_only_compute_along_x,
     $            fluxes_compute_and_exchange_along_x
        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary layers along the
        !> x-direction (E and W) using wall and reflection b.c.
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> table containing the gridpoint data
        !
        !>@param p_model
        !> physical model
        !--------------------------------------------------------------
        subroutine only_compute_along_x(
     $       nodes, prefactor_reflection, prefactor_wall)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(ne), intent(in)  :: prefactor_reflection
          integer    , dimension(ne), intent(in)  :: prefactor_wall

          
          !< prefactor : equal to -1 or +1 depending on the variable
          !>               type: vector_x/vector_y or not
          integer(ikind)         :: i,j
          integer                :: k


          !< E layer: wall b.c.
          !> W layer: reflection b.c.
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                  prefactor_reflection(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                  prefactor_wall(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do

        end subroutine only_compute_along_x


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
        !>@param p_model
        !> physical model
        !
        !>@param card_pt
        !> cardinal point identifying the direction in which the data
        !> are sent
        !--------------------------------------------------------------
        subroutine compute_and_exchange_along_x(
     $     this, comm_2d, usr_rank, nodes,
     $     prefactor_reflection, prefactor_wall, card_pt)

          implicit none

          class(mpi_mg_bc)                , intent(in)    :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer, dimension(ne), intent(in) :: prefactor_reflection
          integer, dimension(ne), intent(in) :: prefactor_wall
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


          !< create two requests in one direction: one for sending
          !> and the other one for receiving data from the same
          !> processor. This processor is identified by the cardinal
          !> direction 'card_pt'
          !DEC$ FORCEINLINE RECURSIVE
          mpi_requests = create_requests_for_one_direction(
     $         this, comm_2d, usr_rank, nodes, card_pt)

          !< overlap some communications with computations
          select case(card_pt)

            !< if we send along E, we compute W
            case(E)
               do k=1,ne
                  do j=1+bc_size, ny-bc_size
                     !DEC$ IVDEP
                     do i=1,bc_size
                   
                        nodes(i,j,k) = 
     $                       prefactor_reflection(k)*
     $                       nodes(2*bc_size+1-i,j,k)
                        
                     end do
                  end do
               end do

            case(W)
               do k=1,ne
                  do j=1+bc_size, ny-bc_size
                     !DEC$ IVDEP
                     do i=1,bc_size
                   
                        nodes(nx-bc_size+i,j,k) = 
     $                       prefactor_wall(k)*
     $                       nodes(nx-bc_size-i+1,j,k)
                   
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
        subroutine fluxes_only_compute_along_x(
     $     nodes, dx, dy, s_op, flux_x)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s_op
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          integer(ikind) :: id

          !< provide the x-indices modified
          id=nx+1-bc_size

          !< modify the fluxes along the x-direction
          !> at the E and W borders
          !> E border: i= nx-bc_size+1
          !DEC$ FORCEINLINE RECURSIVE
          call compute_wall_flux_x(nodes,dx,dy,s_op,id,flux_x)

        end subroutine fluxes_only_compute_along_x


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
     $     nodes, dx, dy, s_op, card_pt, flux_x)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s_op
          integer                           , intent(in)    :: card_pt
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          integer           :: i

          !select the x-indice of the flux modified: only W
          if(card_pt.eq.W) then
             i=nx-bc_size+1
             call compute_wall_flux_x(nodes,dx,dy,s_op,i,flux_x)
          end if

        end subroutine fluxes_compute_and_exchange_along_x

      end module wall_x_reflection_y_par_module
