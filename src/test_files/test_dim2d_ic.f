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
      ! 19_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_dim2d_ic

        use dim2d_eq_class         , only : dim2d_eq
        use field_class            , only : field
        use nf90_operators_wr_class, only : nf90_operators_wr
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : ikind, rkind

        implicit none
        
        
        !<operators tested
        type(field)               :: field_tested
        type(dim2d_eq)            :: p_model
        real(rkind)               :: time
        type(nf90_operators_wr)   :: nf90_writer

        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter        :: detailled=.true.
        integer(ikind)            :: i,j
        real(rkind) :: x_min, y_min
        integer :: bc_size
        

        !<warning
        if(ne.ne.4) then
           stop 'ne=4 is required'
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        field_tested%dx=0.01
        field_tested%dy=0.01
        x_min = -4
        y_min = -4

        do i=1, nx
           field_tested%x_map(i)=x_min + (i-1)*field_tested%dx
        end do

        do j=1, ny
           field_tested%y_map(j)=y_min + (j-1)*field_tested%dy
        end do


        !<test the operators
        time=0
        call p_model%apply_ic(field_tested)
        

        !<write the output data
        call nf90_writer%initialize()
        bc_size=2
        call nf90_writer%write_data(field_tested,p_model,bc_size,time)
        print '(''please check output data file data0.nc'')'

        !<get the last CPU time
        call CPU_TIME(time2)
        print *, 'time elapsed: ', time2-time1


      end program test_dim2d_ic
