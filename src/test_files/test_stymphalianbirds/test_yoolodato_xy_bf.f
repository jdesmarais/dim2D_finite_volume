      ! compare the application of the optimized function
      !  - bc_operators%apply_bc_on_time_dev
      ! 
      ! and the application of the functions for the
      ! boundary layers of the buffer layers
      ! - bc_operators%apply_bc_on_timedev_x_edge
      ! - bc_operators%apply_bc_on_timedev_y_edge
      ! - bc_operators%apply_bc_on_timedev_xy_corner
      ! 
      ! for the Yoo and Lodato boundary conditions
      !---------------------------------------------------------
      program test_yoolodato_xy_bf

        use test_openbc_local_operators_module, only :
     $     test_openbc_local_operators

        implicit none

        logical :: detailled

        detailled = .false.

        call test_openbc_local_operators(detailled)        

      end program test_yoolodato_xy_bf
