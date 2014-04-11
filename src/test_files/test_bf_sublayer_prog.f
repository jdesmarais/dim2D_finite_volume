      program test_bf_sublayer_prog

        use bf_sublayer_class, only : bf_sublayer

        implicit none


        type(bf_sublayer), pointer :: bf_sublayer0

        
        !initialize the bf sublayer
        print '()'
        print '(''test ini() : initialization of bf_sublayer'')'
        print '(''------------------------------------------'')'
        allocate(bf_sublayer0)
        call bf_sublayer0%ini()
        print '(''prev associated?: '',L1)', associated(bf_sublayer0%prev)
        print '(''next associated?: '',L1)', associated(bf_sublayer0%next)
        print '()'

      end program test_bf_sublayer_prog
