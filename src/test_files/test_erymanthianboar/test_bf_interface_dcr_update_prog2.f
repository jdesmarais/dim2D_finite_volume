      program test_bf_interface_dcr_update_prog2

        use bf_interface_dcr_class, only : bf_interface_dcr
        use parameters_input      , only : nx,ny,ne
        use parameters_kind       , only : rkind
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
        integer       , parameter           :: test_case_set=6
        integer       , parameter           :: nb_bubbles=2
        type(test_case)                     :: test_case_used


        !variables tested
        type(bf_interface_dcr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id


        !local variables
        integer :: i
        real(rkind) :: dx
        real(rkind) :: dy
        integer :: file_index
        integer :: timestep

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
     $          interior_nodes, dx, dy)
        
           call test_case_used%update_nodes(
     $          dx,dy,
     $          interior_nodes, interface_used)
        
!           if(mod(timestep,1).eq.0) then
!           if(mod(timestep,5).eq.0) then
              call test_case_used%print_state(
     $             interior_nodes, grdpts_id,
     $             interface_used,
     $             file_index)
              file_index = file_index+1
!           end if
        
           timestep = timestep+1
        
        end do

      end program test_bf_interface_dcr_update_prog2
