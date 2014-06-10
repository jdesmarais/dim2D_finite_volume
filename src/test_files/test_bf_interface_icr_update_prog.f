      program test_bf_interface_icr_update_prog

        use bf_detector_i_list_class, only : bf_detector_i_list
        use bf_interface_icr_class  , only : bf_interface_icr
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : rkind
        use test_bf_layer_module    , only : print_interior_data,
     $                                       ini_grdpts_id
        use test_cases_interface_update_module, only : update_nodes

        implicit none


        type(bf_interface_icr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id

        integer :: i

        real(rkind) :: dx
        real(rkind) :: dy

        integer :: timestep

        dx = 1.0d0/(nx-2*bc_size)
        dy = 1.0d0/(ny-2*bc_size)


        timestep = 0

        !initialization of the grdpts_id
        call ini_grdpts_id(grdpts_id)

        !initialization of the interface
        call interface_used%ini()

        !initialization of the nodes
        call update_nodes(
     $       timestep,dx,dy,
     $       interior_nodes, interface_used)

        !printing the initial state
        call print_state(
     $       interior_nodes, grdpts_id,
     $       interface_used,
     $       timestep)

        timestep = timestep+1


        !update the buffer layers
        do i=1,9

           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy)

           call update_nodes(
     $          timestep,dx,dy,
     $          interior_nodes, interface_used)

           call print_state(
     $          interior_nodes, grdpts_id,
     $          interface_used,
     $          timestep)

           timestep = timestep+1

        end do        

        contains

        subroutine print_state(nodes, grdpts_id, interface_used, index)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface_icr)         , intent(in)    :: interface_used
          integer                         , intent(inout) :: index


          character(len=20) :: i_nodes_filename
          character(len=24) :: i_grdpts_id_filename
          character(len=20) :: i_sizes_filename

          character(len=10) :: bf_nodes_filename
          character(len=13) :: bf_grdpts_id_filename
          character(len=10) :: bf_sizes_filename
          character(len=5) :: bf_nb_sbf_filename


          write(i_nodes_filename,
     $         '(''interior_nodes'',I1,''.dat'')') index
          write(i_grdpts_id_filename,
     $         '(''interior_grdpts_id'',I1,''.dat'')') index
          write(i_sizes_filename,
     $         '(''interior_sizes'',I1,''.dat'')') index
                    
          write(bf_nodes_filename,
     $         '(''nodes'',I1,''.dat'')') index
          write(bf_grdpts_id_filename,
     $         '(''grdpt_id'',I1,''.dat'')') index
          write(bf_sizes_filename,
     $         '(''sizes'',I1,''.dat'')') index
          write(bf_nb_sbf_filename,
     $         '(I1,''.dat'')') index


          call print_interior_data(
     $         nodes,
     $         grdpts_id,
     $         i_nodes_filename,
     $         i_grdpts_id_filename,
     $         i_sizes_filename)
          
          call interface_used%print_binary(
     $         bf_nodes_filename,
     $         bf_grdpts_id_filename,
     $         bf_sizes_filename,
     $         bf_nb_sbf_filename)

          call interface_used%print_idetectors_on_binary(
     $         index)

        end subroutine print_state


      end program test_bf_interface_icr_update_prog
