      !> @file
      !> module encapsulating subroutines to update the MPI
      !> derived types used to exchange data between tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to update the MPI
      !> derived types used to exchange data between tiles
      !> If the tile is exchanging with the same processor
      !> for the N and S boundary layers (or the E and W)
      !> then the data are combined into one structure
      !
      !> @date
      ! 26_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_mg_construct

        use mpi
        use mpi_process_class  , only : mpi_process
        use parameters_constant, only : N,S,E,W,only_exchange_proc

        implicit none

        private
        public :: update_mpi_derived_types

        
        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine updating the mpi derived types
        !> for the exchange of data between the tiles
        !
        !> @date
        !> 26_08_2012 - initial version - J.L. Desmarais
        !
        !> @param com_send
        !> MPI derived types for sending data
        !
        !> @param com_recv
        !> MPI derived types for receiving data
        !
        !>@param proc_x_choice
        !> type of procedure for computing the boundary layers
        !> in the x-direction
        !
        !>@param proc_y_choice
        !> type of procedure for computing the boundary layers
        !> in the y-direction
        !--------------------------------------------------------------
        subroutine update_mpi_derived_types(
     $       com_recv,
     $       com_send,
     $       com_rank,
     $       proc_x_choice,
     $       proc_y_choice)

          implicit none

          integer, dimension(4), intent(inout) :: com_recv
          integer, dimension(4), intent(inout) :: com_send
          integer, dimension(4), intent(in)    :: com_rank
          integer              , intent(in)    :: proc_x_choice
          integer              , intent(in)    :: proc_y_choice

          
          
          
          !> send_struc_x type for sending along x  
          !> recv_struc_x type for receiving along x
          !> send_struc_y type for sending along y  
          !> recv_struc_y type for receiving along y
          !---------------------------------------------
          type(mpi_process)                       :: mpi_op
          integer                                 :: send_struc_x 
          integer                                 :: recv_struc_x 
          integer                                 :: send_struc_y 
          integer                                 :: recv_struc_y 
          integer                                 :: ierror
          integer(MPI_ADDRESS_KIND), dimension(2) :: array_displacements


          !< if the tile is only exchanging then it may
          !> exchange with the same processor for the
          !> E and W boundary layers
          !> if this is the case, it is more interesting
          !> to combine the two MPI structures into one
          !---------------------------------------------
          if((proc_x_choice.eq.only_exchange_proc).and.
     $         (com_rank(E).eq.com_rank(W))) then
             
             !< we combine the two MPI subarrays
             !> for the E and W boundary layers
             !> into one structure and create
             !> a new MPI derived type for it
             !----------------------------------------
             !< we create a structure with 2 previous
             !> MPI derived types, there is only one
             !> of each (/1,1), and no extra memory
             !> space between them (/0,0/). The two
             !> derived types combined are the one
             !> for the West boundary layer
             !> (com_send(W)) and the one for the East
             !> boundary layer (com_send(E)). We save
             !> the new structure into send_struct_x
             !> and redirect the error to ierror
             !----------------------------------------

             !< combining the derived types for sending
             array_displacements(1)=0
             array_displacements(2)=0
             call MPI_TYPE_CREATE_STRUCT(
     $            2, (/1,1/), array_displacements,
     $            (/com_send(W),com_send(E)/),
     $            send_struc_x,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< commiting the new derived type
             call MPI_TYPE_COMMIT(send_struc_x,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< combining the derived types for receiving
             array_displacements(1)=0
             array_displacements(2)=0
             call MPI_TYPE_CREATE_STRUCT(
     $            2,(/1,1/),array_displacements,
     $            (/com_recv(W),com_recv(E)/),
     $            recv_struc_x, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< commiting the new derived type
             call MPI_TYPE_COMMIT(recv_struc_x,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if


             !< the new structures are saved into
             !> com_send(E) and com_recv(E)
             !> com_send(W) and com_recv(W) are no
             !> longer used
             com_send(E)=send_struc_x
             com_send(W)=-999

             com_recv(E)=recv_struc_x
             com_recv(W)=-999

          end if


          !< if the tile is only exchanging then it may
          !> exchange with the same processor for the
          !> N and S boundary layers
          !> if this is the case, it is more interesting
          !> to combine the two MPI structures into one
          !---------------------------------------------
          if((proc_y_choice.eq.only_exchange_proc).and.
     $         (com_rank(N).eq.com_rank(S))) then

             
             !< combining the derived types for sending
             array_displacements(1)=0
             array_displacements(2)=0
             call MPI_TYPE_CREATE_STRUCT(
     $            2,(/1,1/),array_displacements,
     $            (/com_send(N),com_send(S)/),
     $            send_struc_y, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< commiting the new derived type
             call MPI_TYPE_COMMIT(send_struc_y,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< combining the derived types for receiving
             array_displacements(1)=0
             array_displacements(2)=0
             call MPI_TYPE_CREATE_STRUCT(
     $            2,(/1,1/),array_displacements,
     $            (/com_recv(S),com_recv(N)/),
     $            recv_struc_y, ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if

             !< commiting the new derived type
             call MPI_TYPE_COMMIT(recv_struc_y,ierror)
             if(ierror.ne.MPI_SUCCESS) then
                call mpi_op%finalize_mpi()
                print '(''update_mpi_derived_types'')'
                stop 'MPI_TYPE_CREATE_STRUCT failed'
             end if


             !< the new structures are saved into
             !> com_send(N) and com_recv(N)
             !> com_send(S) and com_recv(S) are no
             !> longer used
             com_send(N)=send_struc_y
             com_send(S)=-999

             com_recv(N)=recv_struc_y
             com_recv(S)=-999

          end if

        end subroutine update_mpi_derived_types

      end module mpi_mg_construct
