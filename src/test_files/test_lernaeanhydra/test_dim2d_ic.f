      !> @file
      !> test file for initialization subroutines
      !> for the object 'dim2d_eq'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the initialization subroutines for the
      !> object 'dim2d_eq' by writing the
      !> initial state on a netcdf output file
      !
      !> @date
      ! 19_08_2013 - initial version              - J.L. Desmarais
      ! 15_07_2014 - composition over inheritance - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_dim2d_ic

        use io_operators_class, only : io_operators
        use pmodel_eq_class   , only : pmodel_eq
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : ikind, rkind

        implicit none
        
        
        !<operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind), dimension(nx)       :: x_map
        real(rkind), dimension(ny)       :: y_map
        real(rkind)                      :: time
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        real(rkind)                      :: x_min
        real(rkind)                      :: y_min
        type(pmodel_eq)                  :: p_model
        type(io_operators)               :: nf90_writer
        integer(ikind) :: i,j

        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter :: detailled=.true.
        

        !<warning
        if(ne.ne.4) then
           stop 'ne=4 is required'
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        dx=0.01
        dy=0.01
        x_min = 0
        y_min = 0

        do i=1, nx
           x_map(i)=x_min + (i-1)*dx
        end do

        do j=1, ny
           y_map(j)=y_min + (j-1)*dy
        end do


        !<test the operators
        time=0
        call p_model%apply_ic(nodes,x_map,y_map)
        

        !<write the output data
        call nf90_writer%ini()
        call nf90_writer%write_data(nodes,x_map,y_map,p_model,time)
        print '(''please check output data file data0.nc'')'

        !<get the last CPU time
        call CPU_TIME(time2)
        print *, 'time elapsed: ', time2-time1


      end program test_dim2d_ic
