      !> @file
      !> test file for the module 'mpi_mg_construct'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the initialization of the mpi construct
      !> by exchanging information between two tiles
      !
      !> @date
      ! 26_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_mpi_mg_construct
      
        use mpi
        use mpi_mg_bc_ext_class, only : mpi_mg_bc_ext
        use mpi_mg_construct   , only : update_mpi_derived_types
        use mpi_process_class  , only : mpi_process
        use parameters_constant, only : periodic_xy_choice,N,E
        use parameters_input   , only : nx,ny,ne,npx,npy,bc_choice
        use parameters_kind    , only : ikind, rkind

        implicit none


        !< operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        type(mpi_process)  :: mpi_op
        type(mpi_mg_bc_ext):: mpi_mg


        !< intermediate variables
        integer                             :: comm_2d
        integer                             :: usr_rank
        integer                             :: ierror,sendtag,recvtag
        integer, dimension(MPI_STATUS_SIZE) :: status
        logical, parameter                  :: test=.true.
        logical                             :: test_validated

        
        !< the test is designed for (npx,npy)=(2,2)
        !> and periodic boundary conditions
        if((npx.ne.2).or.(npy.ne.2).or.
     $       (bc_choice.ne.periodic_xy_choice)) then
           print '(''the test needs (npx,npy)=(2,2)'')'
           stop 'and bc_choice=periodic_xy_choice'
        end if


        !< initialization of the mpi process
        call mpi_op%ini_mpi()


        !< initialization of the cartesian communicator
        call mpi_op%ini_cartesian_communicator(comm_2d, usr_rank)

        
        !< initialization of 'mpi_messenger_bc'
        !> with the update of the MPI derived types
        call mpi_mg%ini(comm_2d)


        !< test the update of the MPI derived types
        call update_mpi_derived_types(
     $       mpi_mg%com_recv, mpi_mg%com_send,
     $       mpi_mg%com_rank,
     $       mpi_mg%proc_x_choice, mpi_mg%proc_y_choice)
        

        !< initialize the data saved in f_tested
        nodes = ini_data(usr_rank)


        !< test the exchange of data in the x-direction
        !> to check if the MPI structure works
        select case(usr_rank)

          case(0)

             sendtag = 123
             call MPI_SEND(
     $            nodes, 1, mpi_mg%com_send(E), 2, sendtag,
     $            comm_2d, ierror)
             
          case(2)
             recvtag = 123
             call MPI_RECV(
     $            nodes, 1, mpi_mg%com_recv(E), 0, recvtag,
     $            comm_2d, status, ierror)
             
          case(1)
             
             recvtag=124
             call MPI_RECV(
     $            nodes, 1, mpi_mg%com_recv(E), 3, recvtag,
     $            comm_2d, status, ierror)
             
          case(3)
             sendtag=124
             call MPI_SEND(
     $            nodes, 1, mpi_mg%com_send(E), 1, sendtag,
     $            comm_2d, ierror)

          case default
             call mpi_op%finalize_mpi()
             stop 'problem with nb of procs'
        end select

        if(ierror.ne.MPI_SUCCESS) then
           call mpi_op%finalize_mpi()
           stop 'MPI_SEND or RECV Failed'
        end if

        if(.not.test) then
           call write_data('test_constructx',usr_rank,nodes)
        else
           test_validated = compare_data(
     $          'test_constructx',usr_rank,nodes)

           print '(''Proc '', I1, '' : exchange_x: '', L1)',
     $          usr_rank, test_validated
        end if
        

        !< test the exchange of data in the y-direction
        !> to check if the MPI structure works
        select case(usr_rank)

          case(0)

             sendtag = 123
             call MPI_SEND(
     $            nodes, 1, mpi_mg%com_send(N), 1, sendtag,
     $            comm_2d, ierror)
             
          case(2)
             recvtag = 124
             call MPI_RECV(
     $            nodes, 1, mpi_mg%com_recv(N), 3, recvtag,
     $            comm_2d, status, ierror)
             
          case(1)
             
             recvtag=123
             call MPI_RECV(
     $            nodes, 1, mpi_mg%com_recv(N), 0, recvtag,
     $            comm_2d, status, ierror)
             
          case(3)
             sendtag=124
             call MPI_SEND(
     $            nodes, 1, mpi_mg%com_send(N), 2, sendtag,
     $            comm_2d, ierror)

          case default
             call mpi_op%finalize_mpi()
             stop 'problem with nb of procs'
        end select

        if(ierror.ne.MPI_SUCCESS) then
           call mpi_op%finalize_mpi()
           stop 'MPI_SEND or RECV Failed'
        end if

        if(.not.test) then
           call write_data('test_constructy',usr_rank,nodes)
        else
           test_validated = compare_data(
     $          'test_constructy',usr_rank,nodes)

           print '(''Proc '', I1, '' : exchange_y: '', L1)',
     $          usr_rank, test_validated
        end if

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the initial data for the nodes
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param proc_rank
        !> rank of the processor computing the tile
        !
        !>@param nodes
        !> nodes computed on the tile
        !--------------------------------------------------------------
        function ini_data(proc_rank) result(nodes)

          implicit none

          integer             , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne) :: nodes

          integer(ikind) :: i,j
          integer        :: k

          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k)=compute_ini_data(proc_rank,i,j,k)
                end do
             end do
          end do

        end function ini_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the initial data for the nodes
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param proc_rank
        !> rank of the processor computing the tile
        !
        !>@param i
        !> integer identifying the first coordinate
        !
        !>@param j
        !> integer identifying the second coordinate
        !
        !>@param k
        !> integer identifying the third coordinate
        !
        !>@param var
        !> nodes computed on the tile
        !--------------------------------------------------------------
        function compute_ini_data(proc_rank,i,j,k) result(var)
          implicit none

          integer(ikind), intent(in) :: i,j
          integer       , intent(in) :: proc_rank,k
          real(rkind)         :: var

          var = 1000*proc_rank+100*k+10*j+i
          !var = proc_rank

        end function compute_ini_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> write the data computed on the tile
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param proc_rank
        !> rank of the processor computing the tile
        !
        !>@param nodes
        !> nodes computed on the tile
        !--------------------------------------------------------------
        subroutine write_data(filename_base,proc_rank,nodes)
          implicit none

          character(len=15)                , intent(in) :: filename_base
          integer                         , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes


          character(len=22) :: filename
          integer(ikind) :: i,j


          write(filename,'(A15,''_'',I1,''.txt'')')
     $         filename_base, proc_rank

          open(unit=11,
     $         file=filename,
     $         status='unknown',
     $         position='rewind')

          do j=1, ny
             do i=1,nx
                write(11,'(4F20.6)')
     $               nodes(i,j,1), nodes(i,j,2),
     $               nodes(i,j,3), nodes(i,j,4)
             end do
          end do

          close(11)

        end subroutine write_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> read the data saved in a text file
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param proc_rank
        !> rank of the processor computing the tile
        !
        !>@param nodes
        !> nodes computed on the tile
        !--------------------------------------------------------------
        subroutine read_data(filename_base,proc_rank,nodes)
          implicit none

          character(len=15)                , intent(in)   :: filename_base
          integer                         , intent(in)   :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(inout):: nodes


          character(len=33) :: filename

          integer(ikind) :: i,j

          write(filename, '(''./data_test/'', A15,''_'',I1,''.txt'')')
     $         filename_base, proc_rank

          open(unit=11,
     $         file=filename,
     $         status='unknown',
     $         position='rewind')

          do j=1, ny
             do i=1,nx
                read(11,'(4F20.6)')
     $               nodes(i,j,1), nodes(i,j,2),
     $               nodes(i,j,3), nodes(i,j,4)
             end do
          end do

          close(11)

        end subroutine read_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compare the data with previously computed data
        !> that have been validated by the user
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param proc_rank
        !> rank of the processor computing the tile
        !
        !>@param nodes
        !> nodes computed on the tile
        !
        !>@param test_validated
        !> logical indicating if the two sets of data correspond
        !--------------------------------------------------------------
        function compare_data(filename_base,proc_rank,nodes)
     $     result(test_validated)
          implicit none

          character(len=15)               , intent(in) :: filename_base
          integer                         , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          logical                                      :: test_validated

          real(rkind), dimension(nx,ny,ne) :: test_nodes

          integer(ikind) :: i,j
          integer        :: k

          call read_data(filename_base,proc_rank,test_nodes)

          k=1
          test_validated=.true.
          do while (test_validated.and.(k.le.ne))
             j=1
             do while (test_validated.and.(j.le.ny))
                i=1
                do while (test_validated.and.(i.le.nx))
                   test_validated=test_nodes(i,j,k).eq.nodes(i,j,k)
                   i=i+1
                end do
                j=j+1
             end do
             k=k+1
          end do

        end function compare_data

      end program test_mpi_mg_construct
