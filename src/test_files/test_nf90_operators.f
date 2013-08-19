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

        use field_class            , only : field
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : ikind, rkind
        use dim2d_eq_class         , only : dim2d_eq
        use nf90_operators_wr_class, only : nf90_operators_wr

        implicit none

        
        !<operators tested
        type(field) :: field_tested
        type(dim2d_eq)         :: p_model
        real(rkind), parameter :: time=3.0
        type(nf90_operators_wr):: nf90_writer


        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        integer(ikind)             :: i,j
        logical                    :: test_validated

        
        !<if nx<4, ny<4 then the test cannot be done
        if(ne.ne.4) then
           stop 'the test needs: ne=4'
        end if
        

        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        field_tested%dx=0.6
        field_tested%dy=0.7

        do j=1, ny
           do i=1, nx
              field_tested%nodes(i,j,1) = i + (j-1)*nx
           end do
        end do

        do i=1, nx
           field_tested%x_map(i)=(i-1)*field_tested%dx
        end do

        do j=1, ny
           field_tested%y_map(j)=(j-1)*field_tested%dy
        end do


        !<write the data
        call nf90_writer%initialize()
        call nf90_writer%write_data(field_tested,p_model,time)

      end program test_nf90_operators
