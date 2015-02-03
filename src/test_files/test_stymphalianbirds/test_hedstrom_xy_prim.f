      program test_hedstrom_xy_prim

        use hedstrom_xy_module, only :
     $       compute_x_timedev_with_openbc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0

        implicit none


        real(rkind) :: t
        real(rkind) :: x
        real(rkind) :: y
        
        real(rkind), dimension(3,1,ne) :: nodes_L0
        real(rkind), dimension(3,1,ne) :: nodes_R0


        type(pmodel_eq) :: p_model

        real(rkind) :: dx

        real(rkind), dimension(ne) :: timedev_L0
        real(rkind), dimension(ne) :: timedev_R0

        
        ! initialization of the nodes
        nodes_L0 = reshape((/
     $        1.46d0,  1.27d0,  1.47d0,
     $       0.146d0, 0.143d0, 0.145d0,
     $        0.01d0, 0.002d0,  0.05d0,
     $        4.89d0,  4.87d0,  4.86d0/),
     $       (/3,1,ne/))

        nodes_R0 = reshape((/
     $        1.47d0,  1.27d0,  1.46d0,
     $      -0.145d0,-0.143d0,-0.146d0,
     $        0.05d0, 0.002d0,  0.01d0,
     $        4.86d0,  4.87d0,  4.89d0/),
     $       (/3,1,ne/))

        t  = 0.0
        x  = 1.0
        y  = 1.5
        dx = 0.1

        timedev_L0 = compute_x_timedev_with_openbc(
     $       t,x,y,
     $       nodes_L0, 1,1,
     $       p_model, dx,
     $       gradient_x_x_oneside_L0,
     $       incoming_left)

        timedev_R0 = compute_x_timedev_with_openbc(
     $       t,x,y,
     $       nodes_R0, 3,1,
     $       p_model, dx,
     $       gradient_x_x_oneside_R0,
     $       incoming_right)

        print '(''timedev_L0'')'
        print *, timedev_L0
        print '()'

        print '(''timedev_R0'')'
        print *, timedev_R0
        print '()'

      end program test_hedstrom_xy_prim
