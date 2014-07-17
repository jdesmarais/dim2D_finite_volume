      program test_bf_interface_dcr_update_prog2

        use bf_interface_dcr_class, only : bf_interface_dcr
        use parameters_input      , only : npx,npy,ntx,nty,
     $                                     nx,ny,ne,
     $                                     x_min,x_max,
     $                                     y_min,y_max,
     $                                     dt,
     $                                     search_nb_dt,
     $                                     search_dcr
        use parameters_kind       , only : rkind
        use pmodel_eq_class       , only : pmodel_eq
        use test_bf_layer_module  , only : ini_grdpts_id
        use test_case_class       , only : test_case

        implicit none

        !parametrization of the test case studied
        !----------------------------------------------------------------
        ! test case | nb_bubbles | description                          |
        !----------------------------------------------------------------
        ! 1         | 1          | one bubble moving in the E direction |
        !----------------------------------------------------------------
        ! 2         | 1          | one bubble moving in the W direction |
        !----------------------------------------------------------------
        ! 3         | 1          | one bubble moving in the N direction |
        !----------------------------------------------------------------
        ! 4         | 1          | one bubble moving in the S direction |
        !----------------------------------------------------------------
        ! 5         | 1          | one bubble moving in the NE direction|
        !----------------------------------------------------------------
        ! 6         | 1          | one bubble moving in the SE direction|
        !----------------------------------------------------------------
        ! 7         | 1          | one bubble moving in the NW direction|
        !----------------------------------------------------------------
        ! 8         | 1          | one bubble moving in the SW direction|
        !----------------------------------------------------------------
        ! 1         | 2          | two bubbles moving in opposite       |
        !           |            | direction in the x-direction         |
        !----------------------------------------------------------------
        ! 3         | 2          | two bubbles moving in opposite       |
        !           |            | direction in the y-direction         |
        !----------------------------------------------------------------
        ! 5         | 2          | two bubbles moving in opposite       |
        !           |            | direction in the (y=x)-direction     |
        !----------------------------------------------------------------
        ! 6         | 2          | two bubbles moving in opposite       |
        !           |            | direction in the (y=-x)-direction    |
        !----------------------------------------------------------------
        ! 9         | 1          | growing bubble in every direction    |
        !----------------------------------------------------------------
        integer       , parameter           :: test_case_set=9
        integer       , parameter           :: nb_bubbles=1
        type(test_case)                     :: test_case_used


        !variables tested
        type(bf_interface_dcr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        type(pmodel_eq)                     :: p_model


        !local variables
        integer     :: i
        real(rkind) :: dx
        real(rkind) :: dy
        integer     :: file_index
        integer     :: timestep


        !check the inputs
        if(
     $       (npx.ne.1).or.(npy.ne.1).or.
     $       (ntx.ne.40).or.(nty.ne.40).or.
     $       (x_min.ne.-3.0d0).or.(x_max.ne.0.0d0).or.
     $       (y_min.ne.0.0d0).or.(y_max.ne.3.0d0).or.
     $       (dt.ne.0.02d0).or.(ne.ne.4).or.
     $       (search_nb_dt.ne.1).or.(search_dcr.ne.4)) then

           print '(''the test requires: '')'
           print '(''npx = 1'')'
           print '(''npy = 1'')'
           print '(''ntx = 40'')'
           print '(''nty = 40'')'
           print '(''x_min = -3.0d0'')'
           print '(''x_max = 0.0d0'')'
           print '(''y_min = 0.0d0'')'
           print '(''y_max = 3.0d0'')'
           print '(''dt = 0.02d0'')'
           print '(''ne = 4'')'
           print '(''search_nb_dt = 1'')'
           print '(''search_dcr = 4'')'
           stop ''

        end if


        dx = 1.0d0/nx
        dy = 1.0d0/ny
        file_index = 0
        timestep = 0

        !initialization of the test case
        call test_case_used%ini_with_set(test_case_set, nb_bubbles)

        !initialization of the grdpts_id
        call ini_grdpts_id(grdpts_id)

        !initialization of the interface
        call interface_used%ini()

        !initialization of the nodes
        call test_case_used%update_nodes(
     $       dx,dy,
     $       interior_nodes, interface_used)

        !printing the initial state
        call test_case_used%print_state(
     $       interior_nodes, grdpts_id,
     $       interface_used,
     $       file_index)
        file_index = file_index+1
        timestep   = timestep+1


        !update the buffer layers
        do i=1,55
        
           call interface_used%update_bf_layers_with_detector_dcr(
     $          interior_nodes)
        
           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy, p_model)
        
           call test_case_used%update_nodes(
     $          dx,dy,
     $          interior_nodes, interface_used)
        
!           if(mod(timestep,1).eq.0) then
           if(mod(timestep,5).eq.0) then
              call test_case_used%print_state(
     $             interior_nodes, grdpts_id,
     $             interface_used,
     $             file_index)
              file_index = file_index+1
           end if
        
           timestep = timestep+1
        
        end do

      end program test_bf_interface_dcr_update_prog2
