      program test_poinsot_xy_bf

        use test_openbc_local_operators_module, only :
     $     test_openbc_local_operators

        implicit none

        logical :: detailled

        detailled = .false.

        call test_openbc_local_operators(detailled)        

      end program test_poinsot_xy_bf
