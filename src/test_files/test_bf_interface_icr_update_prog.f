      program test_bf_interface_icr_update_prog

        use bf_mainlayer_class        , only : bf_mainlayer
        use bf_sublayer_class         , only : bf_sublayer
        use bf_detector_icr_list_class, only : bf_detector_icr_list
        use bf_interface_icr_class    , only : bf_interface_icr
        use parameters_constant       , only : N,S,E,W
        use parameters_input          , only : nx,ny,ne,bc_size
        use parameters_kind           , only : rkind
        use test_bf_layer_module      , only : print_interior_data,
     $                                         ini_grdpts_id

        use test_cases_interface_update_module, only : update_nodes

        implicit none


        type(bf_interface_icr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer       , parameter           :: test_case=8

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

           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy)

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


        !remove a buffer layer
        call test_remove_sublayers(
     $       test_case,
     $       interior_nodes,
     $       grdpts_id,
     $       interface_used,
     $       file_index)

                

        contains

        subroutine test_remove_sublayers(
     $       test_case,
     $       nodes,
     $       grdpts_id,
     $       interface_used,
     $       file_index)

          implicit none
          
          integer                         , intent(in)    :: test_case
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface_icr)         , intent(inout) :: interface_used
          integer                         , intent(inout) :: file_index

          integer, dimension(2) :: mainlayer_id
          integer :: i

          select case(test_case)
            case(1)
               mainlayer_id = [E,W]
            case(2)
               mainlayer_id = [E,W]
            case(3)
               mainlayer_id = [N,W]
            case(4)
               mainlayer_id = [S,W]
            case(5)
               mainlayer_id = [N,E]
            case(6)
               mainlayer_id = [S,E]
            case(7)
               mainlayer_id = [N,W]
            case(8)
               mainlayer_id = [S,W]
            case default
               print '(''test_bf_interface_icr_update_prog'')'
               print '(''test_remove_sublayers'')'
               stop 'test case not recognized'
          end select

          do i=1,2

             call test_remove_sublayer(
     $            nodes,
     $            grdpts_id,
     $            interface_used,
     $            mainlayer_id(i),
     $            file_index)

          end do

        end subroutine test_remove_sublayers


        subroutine test_remove_sublayer(
     $       nodes,
     $       grdpts_id,
     $       interface_used,
     $       mainlayer_id,
     $       file_index)
        
          implicit none
          
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface_icr)         , intent(inout) :: interface_used
          integer                         , intent(in)    :: mainlayer_id
          integer                         , intent(inout) :: file_index

          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr

          !remove the head sublayer
          mainlayer_ptr => interface_used%get_mainlayer(mainlayer_id)
          if(associated(mainlayer_ptr)) then
             sublayer_ptr => mainlayer_ptr%get_head_sublayer()
             call interface_used%remove_sublayer(sublayer_ptr)
          end if

          !print the interface
          call print_state(nodes, grdpts_id, interface_used, file_index)
          file_index = file_index+1

        end subroutine test_remove_sublayer


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
