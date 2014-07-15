      !> @file
      !> test file for the module 'periodic_xy_par_module'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the computation and exchange subroutines for
      !> periodic xy
      !
      !> @date
      ! 27_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       program test_periodic_xy_par_module

        use mpi_mg_bc_ext_class   , only : mpi_mg_bc_ext
        use mpi_process_class     , only : mpi_process
        use parameters_constant   , only : periodic_xy_choice,
     $                                     x_direction, y_direction
        use parameters_input      , only : nx,ny,ne,npx,npy,bc_choice
        use parameters_kind       , only : ikind, rkind
        use periodic_xy_par_module, only : only_compute_along_x,
     $                                     only_compute_along_y,
     $                                     only_exchange

        implicit none


        !< operators tested
        integer                          :: comm_2d
        integer                          :: usr_rank
        real(rkind), dimension(nx,ny,ne) :: nodes
        type(mpi_process)                :: mpi_op
        type(mpi_mg_bc_ext)              :: mpi_mg


        !< intermediate variables
        integer(ikind)     :: i,j
        integer            :: k
        logical, parameter :: test=.true.
        logical            :: test_validated


        !< the test is designed for (npx,npy)=(2,2)
        !> and periodic boundary conditions
        if((npx.ne.2).or.(npy.ne.2).or.(nx.ne.10).or.(ny.ne.10).or.
     $       (ne.ne.4).or.(bc_choice.ne.periodic_xy_choice)) then
           print '(''the test needs (npx,npy,nx,ny)=(2,2,10,10)'')'
           stop 'and bc_choice=periodic_xy_choice'
        end if


        !< initialize the mpi process
        call mpi_op%ini_mpi()


        !< initialization of the cartesian communicator
        call mpi_op%ini_cartesian_communicator(comm_2d, usr_rank)


        !< initialize the data
        nodes = ini_data(usr_rank)


        !< initialize the 'mpi_mg_bc_ext' object
        !> and especially update the mpi derived types 
        !> with mpi structures
        call mpi_mg%ini(comm_2d)


        !< test only_compute_along_x
        call only_compute_along_x(nodes)
        if(.not.test) then
           call write_data('test_pcx',usr_rank,nodes)
        else
           test_validated = compare_data('test_pcx',usr_rank,nodes)
           print '(''Proc '', I1, '': test only compute_x : '',L1)',
     $          usr_rank,test_validated
        end if

        
        !< reinitialize the data
        nodes = ini_data(usr_rank)


        !< test only_compute_along_y
        call only_compute_along_y(nodes)
        if(.not.test) then
           call write_data('test_pcy',usr_rank,nodes)
        else
           test_validated = compare_data('test_pcy',usr_rank,nodes)
           print '(''Proc '', I1, '': test only compute_y : '',L1)',
     $          usr_rank,test_validated
        end if


        !< reinitialize the data
        nodes = ini_data(usr_rank)


        !< test only_exchange(x)
        call only_exchange(
     $       mpi_mg, comm_2d, usr_rank, nodes, x_direction)
        if(.not.test) then
           call write_data('test_pEx',usr_rank,nodes)
        else
           test_validated = compare_data('test_pEx',usr_rank,nodes)
           print '(''Proc '', I1, '': test only_exchange_x : '',L1)',
     $          usr_rank,test_validated
        end if


        !< reinitialize the data
        nodes = ini_data(usr_rank)


        !< test only_exchange(y)
        call only_exchange(
     $       mpi_mg, comm_2d, usr_rank, nodes, y_direction)
        if(.not.test) then
           call write_data('test_pEy',usr_rank,nodes)
        else
           test_validated = compare_data('test_pEy',usr_rank,nodes)
           print '(''Proc '', I1, '': test only_exchange_y : '',L1)',
     $          usr_rank,test_validated
        end if


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()

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

          integer, intent(in) :: proc_rank,i,j,k
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

          character(len=8)                , intent(in) :: filename_base
          integer                         , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes


          character(len=14) ::filename
          

          write(filename,'(A8,''_'',I1,''.txt'')')
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

          character(len=8)                , intent(in)   :: filename_base
          integer                         , intent(in)   :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(inout):: nodes


          character(len=26) :: filename


          write(filename, '(''./data_test/'', A8,''_'',I1,''.txt'')')
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

          character(len=8)                , intent(in) :: filename_base
          integer                         , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          logical                                      :: test_validated

          real(rkind), dimension(nx,ny,ne) :: test_nodes

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

      end program test_periodic_xy_par_module
