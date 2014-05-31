      program test_dbf_element_prog

        use dbf_element_class, only : dbf_element
        use parameters_kind  , only : ikind

        implicit none

        
        type(dbf_element), pointer   :: element_used
        type(dbf_element), pointer   :: element_used2
        integer(ikind), dimension(2) :: coords


        !first element
        print '(''test ini()'')'
        print '(''----------'')'
        coords(1) = 1
        coords(2) = 2
        allocate(element_used)
        call element_used%ini(coords)
        call element_used%print_on_screen()
        print '()'

        !set next
        coords(1) = 3
        coords(2) = 4
        allocate(element_used2)
        call element_used2%ini(coords)
        call element_used%set_next(element_used2)

        !set prev
        coords(1) = 5
        coords(2) = 6
        allocate(element_used2)
        call element_used2%ini(coords)
        call element_used%set_prev(element_used2)

        print '(''test set_prev() and set_next()'')'
        print '(''------------------------------'')'
        call element_used%print_on_screen()
        print '()'

      end program test_dbf_element_prog
