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
      !> test data are saved in $test_files/data_test
      !> initial conditions: drop_retraction
      !> pmodel_eq: dim2d_eq
      !
      !> @date
      ! 27_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_rk3tvd_dim2d

        use field_abstract_class, only : field_abstract
        use parameters_input    , only : npx,npy,nx,ny,ne
        use parameters_kind     , only : ikind,rkind
        use td_integrator_class , only : td_integrator

        implicit none


        !< operators tested
        type(field_abstract)   :: field_tested
        type(td_integrator)    :: ti
        real(rkind), parameter :: dt=1.0


        !< CPU recorded times
        real    :: time1, time2


        !< intermediate variables
        real(rkind) :: x_min, x_max, y_min, y_max


        !< if nx.ne.20, ny.ne.20 then the test cannot be done
        if((nx.ne.40).or.(ny.ne.40).or.(ne.ne.4)
     $       .or.(npx.ne.1).or.(npy.ne.1)) then
           stop 'the test needs: (nx,ny,ne,npx,npy)=(20,20,4,1,1)'
        end if
        

        !< get the initial CPU time
        call CPU_TIME(time1)


        !< initialize the tables for the field
        x_min   = 0.
        x_max   = 1.
        y_min   = 0.
        y_max   = 1.

        call field_tested%ini()
        call field_tested%apply_bc_on_nodes()


        !< write the data as initial conditions
        call write_data(
     $       'test_dim2d_rk0',
     $       field_tested%get_nodes(),
     $       field_tested%get_x_map(),
     $       field_tested%get_y_map())


        !< integrate the field for dt
        call ti%integrate(field_tested,dt)


        !< write the data after one integration step
        call write_data(
     $       'test_dim2d_rk1',
     $       field_tested%get_nodes(),
     $       field_tested%get_x_map(),
     $       field_tested%get_y_map())


        !< get the total time needed to run the test
        call CPU_TIME(time2)
        print '(''time elapsed: '', F10.6)', time2-time1


        contains

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
        subroutine write_data(filename_base,nodes,x_map,y_map)
          implicit none

          character(len=14)               , intent(in) :: filename_base
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map


          character(len=18) :: filename
          integer(ikind)    :: i,j


          write(filename,'(A14,''.txt'')')
     $         filename_base

          open(unit=11,
     $         file=filename,
     $         status='unknown',
     $         position='rewind')

          do j=1, ny
             do i=1,nx
                write(11,'(6F20.6)')
     $               x_map(i), y_map(j),
     $               nodes(i,j,1), nodes(i,j,2),
     $               nodes(i,j,3), nodes(i,j,4)
             end do
          end do

          close(11)

        end subroutine write_data

      end program test_rk3tvd_dim2d
