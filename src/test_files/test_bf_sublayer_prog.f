      program test_bf_sublayer_prog

        use bf_sublayer_class, only : bf_sublayer

        implicit none


        type(bf_sublayer), pointer :: bf_sublayer0

        
        !initialize the bf sublayer
        print '()'
        print '(''test ini() : initialization of bf_sublayer'')'
        print '(''------------------------------------------'')'
        allocate(bf_sublayer0)
        call bf_sublayer0%ini(1)
        print '(''prev associated?: '',L1)', associated(bf_sublayer0%get_prev())
        print '(''next associated?: '',L1)', associated(bf_sublayer0%get_next())
        print '(''neighbor1 associated?: '',L1)', associated(bf_sublayer0%get_neighbor(1))
        print '(''neighbor2 associated?: '',L1)', associated(bf_sublayer0%get_neighbor(2))
        print '()'

      end program test_bf_sublayer_prog
