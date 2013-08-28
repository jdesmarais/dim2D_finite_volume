      !> @file
      !> program testing the subroutines defined in
      !> 'reflection_xy_module'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> program testing the subroutines defined in
      !> 'reflection_xy_module'
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_reflection_xy_par_module

        use cg_operators_class      , only : cg_operators
        use dim2d_eq_class          , only : dim2d_eq
        use field_par_class         , only : field_par
        use mpi_mg_bc_class         , only : mpi_mg_bc
        use mpi_process_class       , only : mpi_process
        use parameters_constant     , only : reflection_xy_choice,
     $                                       N,S,E,W,
     $                                       x_direction, y_direction
        use parameters_input        , only : nx,ny,ne,npx,npy,bc_choice
        use parameters_kind         , only : ikind, rkind
        use reflection_xy_par_module, only : only_compute_along_x,
     $     				     only_compute_along_y,
     $                                       only_exchange,
     $                                       compute_and_exchange_along_x,
     $                                       compute_and_exchange_along_y

        implicit none


        !< operators tested
        type(field_par)                  :: f_tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        type(cg_operators)               :: s_op
        type(dim2d_eq)                   :: p_model
        type(mpi_process)                :: mpi_op
        type(mpi_mg_bc)                  :: mpi_mg


        !< intermediate variables
        integer(ikind)     :: i,j
        integer            :: k,bc_size
        logical, parameter :: test=.false.
        logical            :: test_validated


        !< the test is designed for (npx,npy)=(2,2)
        !> and periodic boundary conditions
        if((npx.ne.2).or.(npy.ne.2).or.(nx.ne.10).or.(ny.ne.10).or.
     $       (ne.ne.4).or.(bc_choice.ne.reflection_xy_choice)) then
           print '(''the test needs (npx,npy,nx,ny)=(2,2,10,10)'')'
           stop 'and bc_choice=reflection_xy_choice'
        end if


        !< initialization of intermediate variables
        bc_size = s_op%get_bc_size()


        !< initialize the mpi process
        call mpi_op%ini_mpi()


        !< initialization of the cartesian communicator
        call f_tested%ini_cartesian_communicator()


        !< initialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< initialize the 'mpi_mg_bc' object
        call mpi_mg%initialize(f_tested, s_op)


        !< test only_compute_along_x
        call only_compute_along_x(nodes,bc_size,p_model)
        if(.not.test) then
           call write_data('test_cx',f_tested%usr_rank,nodes)
        else
           test_validated = compare_data('test_cx',f_tested%usr_rank,nodes)
           print '(''Proc '', I1, '': test only compute_x : '',L1)',
     $          f_tested%usr_rank,test_validated
        end if

        
        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test only_compute_along_y
        call only_compute_along_y(nodes,bc_size,p_model)
        if(.not.test) then
           call write_data('test_cy',f_tested%usr_rank,nodes)
        else
           test_validated = compare_data('test_cy',f_tested%usr_rank,nodes)
           print '(''Proc '', I1, '': test only compute_y : '',L1)',
     $          f_tested%usr_rank,test_validated
        end if


        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test compute_and_exchange_along_x()
        select case(f_tested%usr_rank)
          case(0,1)
             !< test compute_and_exchange_along_x(E)
             call compute_and_exchange_along_x(
     $            mpi_mg, f_tested, nodes, bc_size, p_model, E)

          case(2,3)
             !< test compute_and_exchange_along_x(W)
             call compute_and_exchange_along_x(
     $            mpi_mg, f_tested, nodes, bc_size, p_model, W)

          case default
             call mpi_op%finalize_mpi()
             stop 'usr_rank not recognized'
        end select
        if(.not.test) then
           call write_data('test_ex',f_tested%usr_rank,nodes)
        else
           test_validated = compare_data('test_ex',f_tested%usr_rank,nodes)
           print '(''Proc '', I1, '': test compute_exchange_x : '',L1)',
     $          f_tested%usr_rank,test_validated
        end if


        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test compute_and_exchange_along_y()
        select case(f_tested%usr_rank)
          case(0,2)
             !< test compute_and_exchange_along_y(N)
             call compute_and_exchange_along_y(
     $            mpi_mg, f_tested, nodes, bc_size, p_model, N)

          case(1,3)
             !< test compute_and_exchange_along_y(S)
             call compute_and_exchange_along_y(
     $            mpi_mg, f_tested, nodes, bc_size, p_model, S)

          case default
             call mpi_op%finalize_mpi()
             stop 'usr_rank not recognized'
        end select
        if(.not.test) then
           call write_data('test_ey',f_tested%usr_rank,nodes)
        else
           test_validated = compare_data('test_ey',f_tested%usr_rank,nodes)
           print '(''Proc '', I1, '': test compute_exchange_y : '',L1)',
     $          f_tested%usr_rank,test_validated
        end if


        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)

        !< the test only_exchange_x and only_exchange_y
        !> cannnot be performed with npx=2 and npy=2

c$$$        !< test only_exchange(x)
c$$$        call only_exchange(
c$$$     $       mpi_mg, f_tested, nodes, x_direction)
c$$$        if(.not.test) then
c$$$           call write_data('test_Ex',f_tested%usr_rank,nodes)
c$$$        else
c$$$           test_validated = compare_data('test_Ex',f_tested%usr_rank,nodes)
c$$$           print '(''Proc '', I1, '': test only_exchange_x : '',L1)',
c$$$     $          f_tested%usr_rank,test_validated
c$$$        end if
c$$$
c$$$
c$$$        !< reinitialize the data
c$$$        nodes = ini_data(f_tested%usr_rank)
c$$$
c$$$
c$$$        !< test only_exchange(y)
c$$$        call only_exchange(
c$$$     $       mpi_mg, f_tested, nodes, y_direction)
c$$$        if(.not.test) then
c$$$           call write_data('test_Ey',f_tested%usr_rank,nodes)
c$$$        else
c$$$           test_validated = compare_data('test_Ey',f_tested%usr_rank,nodes)
c$$$           print '(''Proc '', I1, '': test only_exchange_y : '',L1)',
c$$$     $          f_tested%usr_rank,test_validated
c$$$        end if


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

          character(len=7)                 , intent(in) :: filename_base
          integer                         , intent(in) :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes


          character(len=14) ::filename


          write(filename,'(A7,''_'',I1,''.txt'')')
     $         filename_base, proc_rank

          open(unit=11,
     $         file=filename,
     $         status='unknown',
     $         position='rewind')

          do j=1, ny
             do i=1,nx
                write(11,'(4F10.6)')
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

          character(len=7)                , intent(in)   :: filename_base
          integer                         , intent(in)   :: proc_rank
          real(rkind), dimension(nx,ny,ne), intent(inout):: nodes


          character(len=26) :: filename


          write(filename, '(''./data_test/'', A7,''_'',I1,''.txt'')')
     $         filename_base, proc_rank

          open(unit=11,
     $         file=filename,
     $         status='unknown',
     $         position='rewind')

          do j=1, ny
             do i=1,nx
                read(11,'(4F10.6)')
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

          character(len=7)                , intent(in) :: filename_base
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

      end program test_reflection_xy_par_module
