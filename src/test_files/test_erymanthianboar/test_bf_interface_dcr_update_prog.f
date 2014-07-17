      program test_bf_interface_dcr_update_prog

        use bf_interface_dcr_class    , only : bf_interface_dcr
        use parameters_input          , only : nx,ny,ne,bc_size
        use parameters_kind           , only : rkind
        use pmodel_eq_class           , only : pmodel_eq
        use test_bf_layer_module      , only : ini_grdpts_id

        use test_cases_interface_update_module, only : update_nodes,
     $                                                 print_state

        implicit none


        type(bf_interface_dcr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer       , parameter           :: test_case=8
        type(pmodel_eq) :: p_model

        integer :: i

        real(rkind) :: dx
        real(rkind) :: dy

        integer :: timestep
        integer :: file_index

        dx = 1.0d0/(nx-2*bc_size)
        dy = 1.0d0/(ny-2*bc_size)


        timestep = 0
        file_index = 0

        !initialization of the grdpts_id
        call ini_grdpts_id(grdpts_id)

        !initialization of the interface
        call interface_used%ini()

        !initialization of the nodes
        call update_nodes(
     $       test_case,
     $       timestep,dx,dy,
     $       interior_nodes, interface_used)

        !printing the initial state
        call print_state(
     $       interior_nodes, grdpts_id,
     $       interface_used,
     $       file_index)
        file_index = file_index+1

        timestep = timestep+1


        !update the buffer layers
        do i=1,30

           call interface_used%update_bf_layers_with_detector_dcr(
     $          interior_nodes)

           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy, p_model)

           call update_nodes(
     $          test_case,
     $          timestep,dx,dy,
     $          interior_nodes, interface_used)

           if(mod(timestep,4).eq.0) then
              call print_state(
     $             interior_nodes, grdpts_id,
     $             interface_used,
     $             file_index)
              file_index = file_index+1
           end if

           timestep = timestep+1

        end do

        print *, interior_nodes(1,1,1)

      end program test_bf_interface_dcr_update_prog
