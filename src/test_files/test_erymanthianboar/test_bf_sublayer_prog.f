      program test_bf_sublayer_prog

        use bf_sublayer_class, only : bf_sublayer
        use parameters_kind  , only : rkind

        implicit none


        type(bf_sublayer), pointer :: bf_sublayer0
        real(rkind) :: dx,dy

        
        !initialize the bf sublayer
        print '()'
        print '(''test ini() : initialization of bf_sublayer'')'
        print '(''------------------------------------------'')'
        allocate(bf_sublayer0)
        call bf_sublayer0%ini(1,dx,dy)
        print '(''prev associated?: '',L1)', associated(bf_sublayer0%get_prev())
        print '(''next associated?: '',L1)', associated(bf_sublayer0%get_next())
        print '()'

      end program test_bf_sublayer_prog
