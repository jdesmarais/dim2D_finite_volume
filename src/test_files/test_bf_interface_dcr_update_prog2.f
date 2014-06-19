      program test_bf_interface_dcr_update_prog2

        use bf_interface_dcr_class, only : bf_interface_dcr
        use parameters_input      , only : nx,ny,ne !,bc_size
        use parameters_kind       , only : rkind
        use test_bf_layer_module  , only : ini_grdpts_id
        use test_case_class       , only : test_case

        implicit none


        !parametrization of the test case studied
        integer       , parameter           :: test_case_set=9
        integer       , parameter           :: nb_bubbles=1
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

        dx = 1.0d0/nx !(nx-2*bc_size)
        dy = 1.0d0/ny !(ny-2*bc_size)
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
        do i=1,50
        
           call interface_used%update_bf_layers_with_detector_dcr(
     $          interior_nodes)
        
           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy)
        
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
