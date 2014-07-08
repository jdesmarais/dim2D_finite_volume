      program test_dbf_list_prog

        use dbf_list_class , only : dbf_list
        use parameters_kind, only : ikind

        implicit none

        type(dbf_list)               :: list_tested
        integer(ikind), dimension(2) :: coords
        
        
        call list_tested%ini()


        coords(1) = 1
        coords(2) = 2
        call list_tested%add_to_list(coords)

        coords(1) = 3
        coords(2) = 4
        call list_tested%add_to_list(coords)

        call list_tested%print_on_screen()


      end program test_dbf_list_prog
