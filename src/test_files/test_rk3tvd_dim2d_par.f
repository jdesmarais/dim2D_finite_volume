      !> @file
      !> test file for the object 'rk3tvd' with
      !> the Diffuse Interface Model in 2D
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> this file is useful to create 2 output files
      !> (initial conditions and after 1 integration step)
      !> to compare the results of the serial and the 
      !> parallel code
      !
      !> @date
      ! 27_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_rk3tvd_dim2d_par

        use bc_operators_par_class, only : bc_operators_par
        use cg_operators_class    , only : cg_operators
        use dim2d_eq_class        , only : dim2d_eq
        use field_par_class       , only : field_par
        use fv_operators_par_class, only : fv_operators_par
        use mpi
        use mpi_process_class     , only : mpi_process
        use parameters_constant   , only : periodic_xy_choice
        use parameters_input      , only : nx,ny,ne,npx,npy,
     $                                     x_min,x_max,y_min,y_max,
     $                                     bc_choice
        use parameters_kind       , only : ikind,rkind
        use rk3tvd_par_class      , only : rk3tvd_par

        implicit none

        
        !< operators tested
        type(field_par)        :: field_tested
        type(cg_operators)     :: sd
        type(dim2d_eq)         :: p_model
        type(fv_operators_par) :: td
        type(bc_operators_par) :: bc_used
        type(rk3tvd_par)       :: ti
        type(mpi_process)      :: mpi_op
        real(rkind), parameter :: dt=1.0


        !< CPU recorded times
        real    :: time1, time2


        !< intermediate variables
        integer     :: bc_size
        logical     :: test_coordinates
        logical     :: test_validated


        !< if nx.ne.20, ny.ne.20 then the test cannot be done
        if((nx.ne.12).or.(ny.ne.12).or.(ne.ne.4)
     $       .or.(npx.ne.2).or.(npy.ne.2)) then
           stop 'the test needs: (nx,ny,ne,npx,npy)=(12,12,4,2,2)'
        end if
        

        !< get the initial CPU time
        call CPU_TIME(time1)


        !< initialization of the mpi process
        call mpi_op%ini_mpi()


        !< initialize the tables for the field
        test_coordinates = x_min.eq.0
        test_coordinates = test_coordinates.and.(x_max.eq.1)
        test_coordinates = test_coordinates.and.(y_min.eq.0)
        test_coordinates = test_coordinates.and.(y_max.eq.1)

        if(.not.test_coordinates) then
           stop 'the test needs: (x_min,x_max,y_min,y_max)=(0,1,0,1)'
        end if

        if(bc_choice.ne.periodic_xy_choice) then
           stop 'the test needs: bc_choice=periodic_xy_choice'
        end if

        bc_size = sd%get_bc_size()

        call field_tested%ini_cartesian_communicator()
        call field_tested%ini_coordinates(bc_size)
        call p_model%apply_ic(field_tested)

        call bc_used%initialize(field_tested,sd)
        call bc_used%apply_bc_on_nodes(
     $       field_tested, field_tested%nodes, sd, p_model)


        !< compare the initial conditions data
        test_validated = compare_data('test_dim2d_rk0', field_tested)
        print '(''Proc '', I1, '' : test initial data: '', L1)',
     $       field_tested%usr_rank, test_validated


        !< integrate the field for dt
        call ti%integrate(field_tested,sd,p_model,td,bc_used,dt)


        !< write the data after one integration step
        test_validated = compare_data('test_dim2d_rk1', field_tested)
        print '(''Proc '', I1, '' : test first step data: '', L1)',
     $       field_tested%usr_rank, test_validated


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()


        !< get the total time needed to run the test
        call CPU_TIME(time2)
        print '(''time elapsed: '', F10.6)', time2-time1

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the test data to which the parallel computed data
        !> are compared
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param test_nodes
        !> nodes computed by the procedural program
        !
        !>@param test_x_map
        !> x_map corresponding to the procedural program
        !
        !>@param test_y_map
        !> y_map corresponding to the procedural program
        !--------------------------------------------------------------
        subroutine get_test_data(
     $       filename_base,test_nodes,test_x_map,test_y_map)

          implicit none

          character(len=14)               , intent(in)   :: filename_base
          real(rkind), dimension(20,20,ne), intent(inout):: test_nodes
          real(rkind), dimension(20)      , intent(inout):: test_x_map
          real(rkind), dimension(20)      , intent(inout):: test_y_map

          character(len=30) :: filename

          integer(ikind) :: i,j

          write(filename, '(''./data_test/'', A14,''.txt'')')
     $         filename_base

          open(unit=11,
     $         file=filename,
     $         action='read',
     $         form='formatted',
     $         status='unknown',
     $         position='rewind')

          do j=1, 20
             do i=1, 20
                read(11,'(6F20.6)')
     $               test_x_map(i)    , test_y_map(j),
     $               test_nodes(i,j,1), test_nodes(i,j,2),
     $               test_nodes(i,j,3), test_nodes(i,j,4)
             end do
          end do

          close(11)

        end subroutine get_test_data



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the test data to which the parallel computed data
        !> are compared
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param test_nodes
        !> nodes computed by the procedural program
        !
        !>@param test_x_map
        !> x_map corresponding to the procedural program
        !
        !>@param test_y_map
        !> y_map corresponding to the procedural program
        !--------------------------------------------------------------
        function is_test_validated(var,cst) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: test_validated

          test_validated=(
     $         int(var*1000.)-
     $         sign(int(abs(cst*1000.)),int(cst*1000.))).eq.0
          
        end function is_test_validated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the test data to which the parallel computed data
        !> are compared
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename_base
        !> character giving the name of the file
        !
        !>@param test_nodes
        !> nodes computed by the procedural program
        !
        !>@param test_x_map
        !> x_map corresponding to the procedural program
        !
        !>@param test_y_map
        !> y_map corresponding to the procedural program
        !--------------------------------------------------------------
        function compare_data(
     $     filename_base, f_tested)
     $     result(test_validated)

          implicit none

          character(len=14), intent(in) :: filename_base
          class(field_par) , intent(in) :: f_tested
          logical                       :: test_validated


          type(mpi_process)     :: mpi_op
          integer               :: dims_nb
          integer, dimension(2) :: cart_coord
          integer               :: ierror
          integer               :: offset_i, offset_j
          integer(ikind)        :: i,j
          integer               :: k
          


          real(rkind), dimension(20,20,ne) :: test_nodes
          real(rkind), dimension(20)       :: test_x_map
          real(rkind), dimension(20)       :: test_y_map


          !< get the test data
          call get_test_data(filename_base,test_nodes,test_x_map,test_y_map)
          

          !< get the cartesian coordinates of the field
          dims_nb=2
          call MPI_CART_COORDS(
     $         f_tested%comm_2d, f_tested%usr_rank,
     $         dims_nb, cart_coord,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'ini_coordinates: MPI_CART_RANK failed'
          end if


          !< depending on the cartesian coordinates of the tile
          !> the subarray tested will be different
          offset_i = cart_coord(1)*(nx-2*bc_size)+1
          offset_j = cart_coord(2)*(ny-2*bc_size)+1


          test_validated=.true.
          k=1
          do while (test_validated.and.(k.le.2))
             j=1
             do while (test_validated.and.(j.le.ny))
                i=1
                do while (test_validated.and.(i.le.nx))
                   test_validated =is_test_validated(
     $                  f_tested%nodes(i,j,k),
     $                  test_nodes(offset_i+i-1,offset_j+j-1,k))
                   test_validated = test_validated.and.is_test_validated(
     $                  f_tested%x_map(i),
     $                  test_x_map(offset_i+i-1))
                   test_validated = test_validated.and.is_test_validated(
     $                  f_tested%y_map(j),
     $                  test_y_map(offset_j+j-1))
                   i=i+1
                end do
                j=j+1
             end do
             k=k+1
          end do

          if(.not.test_validated) then
             print '(''(i,j,k)='', I2,1X,I2,1X,I2)',
     $            i,j,k
          end if

         end function compare_data

      end program test_rk3tvd_dim2d_par
