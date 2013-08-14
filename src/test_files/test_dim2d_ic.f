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
        use parameters_kind        , only : ikind, rkind

        implicit none
        
        
        !<operators tested
        type(field)               :: field_tested
        integer(ikind), parameter :: nx=100
        integer(ikind), parameter :: ny=100
        integer       , parameter :: ne=4
        type(dim2d_eq)            :: p_model
        real(rkind)               :: time
        type(nf90_operators_wr)   :: nf90_writer

        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter        :: detailled=.true.
        integer(ikind)            :: i,j
        

        !<get the initial CPU time
        call CPU_TIME(time1)


        !<allocate the tables for the field
        call field_tested%allocate_tables(nx,ny,ne)


        !<initialize the tables for the field
        field_tested%dx=0.01
        field_tested%dy=0.01

        do i=1, size(field_tested%x_map)
           field_tested%x_map(i)=(i-1)*field_tested%dx
        end do

        do j=1, size(field_tested%y_map)
           field_tested%y_map(j)=(j-1)*field_tested%dy
        end do


        !<test the operators
        time=0
        call p_model%apply_ic(field_tested)
        

        !<write the output data
        call nf90_writer%initialize()
        call nf90_writer%write_data(field_tested,p_model,time)


        !<get the last CPU time
        call CPU_TIME(time2)
        print *, 'time elapsed: ', time2-time1


      end program test_dim2d_ic
