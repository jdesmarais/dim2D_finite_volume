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

        use field_abstract_par_class, only : field_abstract_par
        use mpi
        use mpi_process_class       , only : mpi_process
        use parameters_constant     , only : periodic_xy_choice
        use parameters_input        , only : nx,ny,ne,npx,npy,bc_size,
     $                                       x_min,x_max,y_min,y_max,
     $                                       bc_choice
        use parameters_kind         , only : ikind,rkind
        use td_integrator_par_class , only : td_integrator_par

        implicit none

        
        !< operators tested
        type(field_abstract_par) :: field_tested
        type(td_integrator_par)  :: ti
        real(rkind), parameter   :: dt=1.0
        type(mpi_process)        :: mpi_op


        !< CPU recorded times
        real    :: time1, time2


        !< intermediate variables
        logical     :: test_coordinates
        logical     :: test_validated


        !< if nx.ne.12, ny.ne.12 then the test cannot be done
        if((nx.ne.22).or.(ny.ne.22).or.(ne.ne.4)
     $       .or.(npx.ne.2).or.(npy.ne.2)) then
           stop 'the test needs: (nx,ny,ne,npx,npy)=(18,18,4,2,2)'
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

        call field_tested%ini()
        call field_tested%apply_bc_on_nodes()


        !< compare the initial conditions data
        test_validated = compare_data(
     $       'test_dim2d_rk0', 
     $       field_tested%get_comm_2d(),
     $       field_tested%get_usr_rank(),
     $       field_tested%get_nodes(),
     $       field_tested%get_x_map(),
     $       field_tested%get_y_map())

        print '(''Proc '', I1, '' : test initial data: '', L1)',
     $       field_tested%get_usr_rank(), test_validated


        !< integrate the field for dt
        call ti%integrate(field_tested,dt)


        !< write the data after one integration step
        test_validated = compare_data(
     $       'test_dim2d_rk1',
     $       field_tested%get_comm_2d(),
     $       field_tested%get_usr_rank(),
     $       field_tested%get_nodes(),
     $       field_tested%get_x_map(),
     $       field_tested%get_y_map())

        print '(''Proc '', I1, '' : test first step data: '', L1)',
     $       field_tested%get_usr_rank(), test_validated


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()


        !< get the total time needed to run the test
        call CPU_TIME(time2)
        print '(''time elapsed: '', F10.6)', time2-time1


        !print *, field_tested%get_x_map()

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
          real(rkind), dimension(40,40,ne), intent(inout):: test_nodes
          real(rkind), dimension(40)      , intent(inout):: test_x_map
          real(rkind), dimension(40)      , intent(inout):: test_y_map

          character(len=31) :: filename

          integer(ikind) :: i,j

          write(filename, '(''../data_test/'', A14,''.txt'')')
     $         filename_base

          open(unit=11,
     $         file=filename,
     $         action='read',
     $         form='formatted',
     $         status='unknown',
     $         position='rewind')

          do j=1, 40
             do i=1, 40
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

          test_validated=abs(
     $         int(var*1000.)-
     $         sign(int(abs(cst*1000.)),int(cst*1000.))).le.1

          if(.not.test_validated) then
             print *, int(var*1000), sign(int(abs(cst*1000.)),int(cst*1000.))
          end if
          
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
     $     filename_base,
     $     comm_2d, usr_rank,
     $     nodes, x_map, y_map)
     $     result(test_validated)

          implicit none

          character(len=14)               , intent(in) :: filename_base
          integer                         , intent(in) :: comm_2d
          integer                         , intent(in) :: usr_rank
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          logical                                      :: test_validated


          type(mpi_process)     :: mpi_op
          integer               :: dims_nb
          integer, dimension(2) :: cart_coord
          integer               :: ierror
          integer               :: offset_i, offset_j
          integer(ikind)        :: i,j
          integer               :: k
          


          real(rkind), dimension(40,40,ne) :: test_nodes
          real(rkind), dimension(40)       :: test_x_map
          real(rkind), dimension(40)       :: test_y_map


          !< get the test data
          call get_test_data(
     $         filename_base,test_nodes,test_x_map,test_y_map)


          !< get the cartesian coordinates of the field
          dims_nb=2
          call MPI_CART_COORDS(
     $         comm_2d, usr_rank,
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
     $                  nodes(i,j,k),
     $                  test_nodes(offset_i+i-1,offset_j+j-1,k))
                   test_validated = test_validated.and.is_test_validated(
     $                  x_map(i),
     $                  test_x_map(offset_i+i-1))
                   test_validated = test_validated.and.is_test_validated(
     $                  y_map(j),
     $                  test_y_map(offset_j+j-1))
                   i=i+1
                end do
                j=j+1
             end do
             k=k+1
          end do

          if(.not.test_validated) then
             print '(''x_map('',I2,'')='',1X,F8.3,3X,F8.3,1X,58X,
     $               ''y_map('',I2,'')='',1X,F8.3,3X,F8.3,1X,58X,
     $               ''nodes('',I2,'','',I2,'','',I2,'')='',1X,F8.3,3X,F8.3
     $            )',
     $            i, x_map(i), test_x_map(offset_i+i-1),
     $            j, y_map(j), test_y_map(offset_j+j-1),
     $            i,j,k, nodes(i,j,k), test_nodes(offset_i+i-1,offset_j+j-1,k)
          end if

         end function compare_data

      end program test_rk3tvd_dim2d_par
