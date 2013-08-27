      !> @file
      !> test file for the object 'field'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the allocation of the attributes of 'field'
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_field

        use parameters_input, only : nx, ny, ne
        use parameters_kind , only : rkind
        use field_class     , only : field


        implicit none


        type(field) :: field_tested

        logical :: test_validated


        real(rkind) :: x_min, x_max, y_min, y_max
        integer     :: bc_size

        x_min   = 0
        x_max   = 5
        y_min   = 0
        y_max   = 5
        bc_size = 2


        call field_tested%ini_coordinates(x_min,x_max,y_min,y_max,bc_size)
        

        print '(''allocate_tables test'')'
        print '(''--------------------'')'
        test_validated =
     $       (size(field_tested%nodes,1).eq.nx).and.
     $       (size(field_tested%nodes,2).eq.ny).and.
     $       (size(field_tested%nodes,3).eq.ne).and.
     $       (size(field_tested%x_map,1).eq.nx).and.
     $       (size(field_tested%y_map,1).eq.ny)
        print '(''test_validated: '', 1L)', test_validated
        print '('''')'

      end program test_field
