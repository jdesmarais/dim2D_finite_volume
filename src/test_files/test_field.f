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

        use parameters_kind, only : ikind
        use field_class    , only : field


        implicit none


        type(field) :: field_tested

        integer(ikind), parameter :: nx=100
        integer(ikind), parameter :: ny=100
        integer, parameter :: ne=3
        logical :: test_validated


        print '(''allocate_tables test'')'
        print '(''--------------------'')'
        !DEC$ INLINE RECURSIVE
        call field_tested%allocate_tables(nx,ny,ne)
        test_validated =
     $       allocated(field_tested%nodes).and.
     $       allocated(field_tested%x_map).and.
     $       allocated(field_tested%y_map)
        print '(''test_validated: '', 1L)', test_validated
        print '('''')'

      end program test_field
