      !> @file
      !> test file for the object 'nf90_operators'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the subroutines writing the netcdf files
      !
      !> @date
      ! 14_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_nf90_operators

        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : ikind, rkind
        use pmodel_eq_class   , only : pmodel_eq
        use io_operators_class, only : io_operators

        implicit none

        
        !<operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind), dimension(nx)       :: x_map
        real(rkind), dimension(ny)       :: y_map
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        type(pmodel_eq)                  :: p_model
        real(rkind), parameter           :: time=3.0
        type(io_operators)               :: nf90_writer


        !<CPU recorded times
        real :: time1, time2

        !<test parameters
        integer(ikind) :: i,j

        
        !<if nx<4, ny<4 then the test cannot be done
        if(ne.ne.4) then
           stop 'the test needs: ne=4'
        end if
        

        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        dx=0.6
        dy=0.7

        do j=1, ny
           do i=1, nx
              nodes(i,j,1) = i + (j-1)*nx
           end do
        end do

        do i=1, nx
           x_map(i)=(i-1)*dx
        end do

        do j=1, ny
           y_map(j)=(j-1)*dy
        end do


        !<write the data
        call nf90_writer%ini()
        call nf90_writer%write_data(nodes,x_map,y_map,p_model,time)

        !<get the initial CPU time
        call CPU_TIME(time2)


        print '(''time elapsed: '', F10.6)', time2-time1
        print '(''please check output data file data0.nc'')'

      end program test_nf90_operators
