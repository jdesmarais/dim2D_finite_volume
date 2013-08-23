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
        integer(ikind) :: i,j
        integer        :: k,bc_size


        !< test data
        !real(rkind), dimension(nx,ny,ne) :: test_only_compute_along_x
        !real(rkind), dimension(nx,ny,ne) :: test_only_compute_along_y
        !real(rkind), dimension(nx,ny,ne) :: test_only_exchange_along_x
        !real(rkind), dimension(nx,ny,ne) :: test_only_exchange_along_y
        !real(rkind), dimension(nx,ny,ne) :: test_compute_and_exchange_N
        !real(rkind), dimension(nx,ny,ne) :: test_compute_and_exchange_S
        !real(rkind), dimension(nx,ny,ne) :: test_compute_and_exchange_E
        !real(rkind), dimension(nx,ny,ne) :: test_compute_and_exchange_W


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

        
        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test only_compute_along_y
        call only_compute_along_y(nodes,bc_size,p_model)


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


        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test only_exchange(x)
        call only_exchange(
     $       mpi_mg, f_tested, nodes, x_direction)


        !< reinitialize the data
        nodes = ini_data(f_tested%usr_rank)


        !< test only_exchange(y)
        call only_exchange(
     $       mpi_mg, f_tested, nodes, y_direction)


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()

        contains

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


        function compute_ini_data(proc_rank,i,j,k) result(var)
          implicit none

          integer, intent(in) :: proc_rank,i,j,k
          real(rkind)         :: var

          var = proc_rank

        end function compute_ini_data

      end program test_reflection_xy_par_module
