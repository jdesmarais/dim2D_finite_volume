      !> @file
      !> test file for the object 'mpi_mg_bc'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the initialization of the mpi attributes of 'mpi_mg_bc'
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_mpi_mg_bc

        use cg_operators_class , only : cg_operators
        use field_par_class    , only : field_par
        use mpi
        use mpi_mg_bc_class    , only : mpi_mg_bc
        use mpi_process_class  , only : mpi_process
        use parameters_constant, only : periodic_xy_choice,N,S,E,W
        use parameters_input   , only : nx,ny,ne,npx,npy,bc_choice,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none


        !< operators tested
        type(field_par)    :: f_tested
        type(mpi_process)  :: mpi_op
        type(mpi_mg_bc)    :: mpi_mg
        type(cg_operators) :: s_op


        !< intermediate variables
        integer(ikind)                      :: i,j
        integer                             :: k
        integer                             :: ierror,sendtag,recvtag
        integer, dimension(MPI_STATUS_SIZE) :: status
        logical                             :: test_validated


        !< test data
        integer, dimension(4)  :: test_com_rank
        real(rkind), dimension(:,:,:), allocatable :: test_1x_exchange
        real(rkind), dimension(:,:,:), allocatable :: test_1y_exchange

        
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
        call f_tested%ini_cartesian_communicator()


        !< enter the test data for com_rank
        select case(f_tested%usr_rank)
          case(0)
             test_com_rank=[1,1,2,2]
          case(1)
             test_com_rank=[0,0,3,3]
          case(2)
             test_com_rank=[3,3,0,0]
          case(3)
             test_com_rank=[2,2,1,1]
          case default
             call mpi_op%finalize_mpi()
             stop 'proc rank not correct'
        end select

        !< initialize test data for first x exchange
        allocate(test_1x_exchange(bc_size, ny-2*bc_size,ne))

        do k=1,ne
           do j=1, size(test_1x_exchange,2)
              do i=1, size(test_1x_exchange,1)
                 
                 select case(f_tested%usr_rank)
                 
                   case(1)
                      test_1x_exchange(i,j,k) = compute_ini_data(
     $                     3,bc_size+i,bc_size+j,k)
                   case(2)
                      test_1x_exchange(i,j,k) = compute_ini_data(
     $                     0,nx-2*bc_size+i,bc_size+j,k)

                 end select

              end do
           end do
        end do


        !< initialize test data for first y exchange
        allocate(test_1y_exchange(nx,bc_size,ne))

        do k=1,ne
           do j=1, size(test_1y_exchange,2)
              do i=1, size(test_1y_exchange,1)
                 
                 select case(f_tested%usr_rank)
                 
                   case(1)
                      test_1y_exchange(i,j,k) = compute_ini_data(
     $                     0,i,ny-2*bc_size+j,k)
                   case(2)
                      test_1y_exchange(i,j,k) = compute_ini_data(
     $                     3,i,bc_size+j,k)

                 end select

              end do
           end do
        end do


        !< initialize the data in field_par
        !> depending on the processor id
        do k=1, ne
           do j=1, ny
              do i=1,nx
                 f_tested%nodes(i,j,k)=compute_ini_data(f_tested%usr_rank,i,j,k)
              end do
           end do
        end do


        !< test the initialization of 'mpi_messenger_bc'
        call mpi_mg%initialize(f_tested,s_op)


        !< test the com_rank
        test_validated=.true.
        k=1
        do while(test_validated.and.k.le.4)
           test_validated=mpi_mg%com_rank(k).eq.test_com_rank(k)
           k=k+1
        end do
        print '(''proc, '', I1, '' test_com_rank: '', L1)',
     $       f_tested%usr_rank, test_validated


        !< test the exchange of data in the first x-direction
        test_validated=.true.
        select case(f_tested%usr_rank)

          case(0)

             sendtag = 123
             call MPI_SEND(
     $            f_tested%nodes, 1, mpi_mg%com_send(E), 2, sendtag,
     $            f_tested%comm_2d, ierror)
             
          case(2)
             recvtag = 123
             call MPI_RECV(
     $            f_tested%nodes, 1, mpi_mg%com_recv(W), 0, recvtag,
     $            f_tested%comm_2d, status, ierror)
             
          case(1)
             
             recvtag=124
             call MPI_RECV(
     $            f_tested%nodes, 1, mpi_mg%com_recv(E), 3, recvtag,
     $            f_tested%comm_2d, status, ierror)
             
          case(3)
             sendtag=124
             call MPI_SEND(
     $            f_tested%nodes, 1, mpi_mg%com_send(W), 1, sendtag,
     $            f_tested%comm_2d, ierror)             

          case default
             call mpi_op%finalize_mpi()
             stop 'problem with nb of procs'
        end select
        if(ierror.ne.MPI_SUCCESS) then
           call mpi_op%finalize_mpi()
           stop 'MPI_SENDRECV Failed'
        end if

        
        !< check if the exchange along x worked
        test_validated=.true.
        do k=1, ne
           j=1
           do while(test_validated.and.(j.le.size(test_1x_exchange,2)))
              i=1
              do while(test_validated.and.(i.le.size(test_1x_exchange,1)))
                 select case(f_tested%usr_rank)
                    case(1)
                       test_validated=
     $                      f_tested%nodes(i+nx-bc_size,j+bc_size,k).eq.
     $                      test_1x_exchange(i,j,k)
                    case(2)
                       test_validated=
     $                      f_tested%nodes(i,j+bc_size,k).eq.
     $                      test_1x_exchange(i,j,k)
                 end select
                 i=i+1
              end do
              j=j+1
           end do
        end do
        print '(''proc, '', I1, '' test_1x_exchange: '', L1)',
     $       f_tested%usr_rank, test_validated



        !< test the exchange of data in the first y-direction
        test_validated=.true.
        select case(f_tested%usr_rank)

          case(0)

             sendtag = 123
             call MPI_SEND(
     $            f_tested%nodes, 1, mpi_mg%com_send(N), 1, sendtag,
     $            f_tested%comm_2d, ierror)
             
          case(2)
             recvtag = 124
             call MPI_RECV(
     $            f_tested%nodes, 1, mpi_mg%com_recv(N), 3, recvtag,
     $            f_tested%comm_2d, status, ierror)
             
          case(1)
             
             recvtag=123
             call MPI_RECV(
     $            f_tested%nodes, 1, mpi_mg%com_recv(S), 0, recvtag,
     $            f_tested%comm_2d, status, ierror)
             
          case(3)
             sendtag=124
             call MPI_SEND(
     $            f_tested%nodes, 1, mpi_mg%com_send(S), 2, sendtag,
     $            f_tested%comm_2d, ierror)

          case default
             call mpi_op%finalize_mpi()
             stop 'problem with nb of procs'
        end select
        if(ierror.ne.MPI_SUCCESS) then
           call mpi_op%finalize_mpi()
           stop 'MPI_SENDRECV Failed'
        end if


        !< check if the exchange along x worked
        test_validated=.true.
        do k=1, ne
           j=1
           do while(test_validated.and.(j.le.size(test_1y_exchange,2)))
              i=1
              do while(test_validated.and.(i.le.size(test_1y_exchange,1)))
                 select case(f_tested%usr_rank)
                    case(1)
                       test_validated=
     $                      f_tested%nodes(i,j,k).eq.
     $                      test_1y_exchange(i,j,k)
                    case(2)
                       test_validated=
     $                      f_tested%nodes(i,ny-bc_size+j,k).eq.
     $                      test_1y_exchange(i,j,k)
                 end select
                 i=i+1
              end do
              j=j+1
           end do
        end do
        print '(''proc, '', I1, '' test_1y_exchange: '', L1)',
     $       f_tested%usr_rank, test_validated


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()


        contains

        function compute_ini_data(proc_rank,i,j,k) result(var)
          implicit none

          integer, intent(in) :: proc_rank,i,j,k
          real(rkind)         :: var

          var = proc_rank

        end function compute_ini_data

      end program test_mpi_mg_bc
