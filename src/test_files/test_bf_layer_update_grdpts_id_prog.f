      program test_bf_layer_update_grdpts_id_prog

        use ifport

        use bf_interface_class          , only : bf_interface
        use bf_mainlayer_class          , only : bf_mainlayer
        use bf_sublayer_class           , only : bf_sublayer

        use parameters_constant         , only : N,S,E,W
        use parameters_input            , only : nx,ny,ne,bc_size
        use parameters_kind             , only : ikind, rkind

        use test_bf_layer_module        , only : print_interior_data,
     $                                           ini_grdpts_id,
     $                                           ini_nodes
        use test_cases_interface_module , only : ini_interface

        implicit none

        real(rkind), parameter :: background_grdpt = 0.025
        real(rkind), parameter :: identified_grdpt = 0.9

        integer, parameter :: test_case_id = 4
        integer, parameter :: relative_distance = 3
        integer, parameter :: relative_size = 5
        integer, parameter :: over_allocated = 2
        integer, parameter :: random_seed = 87956

        real(rkind), dimension(nx,ny,ne) :: nodes
        integer    , dimension(nx,ny)    :: grdpts_id
        type(bf_interface)               :: interface_used

        integer, dimension(:,:), allocatable :: selected_grdpts_N
        integer, dimension(:,:), allocatable :: selected_grdpts_S
        integer, dimension(:,:), allocatable :: selected_grdpts_E
        integer, dimension(:,:), allocatable :: selected_grdpts_W


        !initialize the nodes for the test
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_interior_data(
     $       nodes,
     $       grdpts_id, 
     $       'interior_nodes.dat',
     $       'interior_grdpts_id.dat',
     $       'interior_sizes.dat')

        !initialize the interface
        call ini_interface(
     $       interface_used,
     $       test_case_id,
     $       nodes,
     $       relative_distance=relative_distance,
     $       relative_size=relative_size,
     $       over_allocated=over_allocated)

        !print interface before
         call interface_used%print_binary(
     $        'nodes_0.dat',
     $        'grdpt_id_0.dat',
     $        'sizes_0.dat',
     $        '_0.dat')

        !for each mainlayer determine the gridpoints
        !that should be turned from bc_interior to
        !interior
        call determine_selected_grdpts(
     $       random_seed,
     $       N,
     $       relative_size,
     $       over_allocated,
     $       selected_grdpts_N)
        call determine_selected_grdpts(
     $       random_seed,
     $       S,
     $       relative_size,
     $       over_allocated,
     $       selected_grdpts_S)
        call determine_selected_grdpts(
     $       random_seed,
     $       E,
     $       relative_size,
     $       over_allocated,
     $       selected_grdpts_E)
        call determine_selected_grdpts(
     $       random_seed,
     $       W,
     $       relative_size,
     $       over_allocated,
     $       selected_grdpts_W)


        !identify the selected gridpoints
        !in the mainlayer_gridpoints
        call identify_selected_grdpts(
     $       interface_used,
     $       N,
     $       selected_grdpts_N)
        
         call identify_selected_grdpts(
     $       interface_used,
     $       S,
     $       selected_grdpts_S)
        
         call identify_selected_grdpts(
     $       interface_used,
     $       E,
     $       selected_grdpts_E)

         call identify_selected_grdpts(
     $       interface_used,
     $       W,
     $       selected_grdpts_W)


         !print interface with selected
         !grdpoints identified
         call interface_used%print_binary(
     $        'nodes_1.dat',
     $        'grdpt_id_1.dat',
     $        'sizes_1.dat',
     $        '_1.dat')


         !turn the local coordinates of the gridpoints
         !into general coordinates
         call from_local_to_general_coord(
     $        interface_used,
     $        N,
     $        selected_grdpts_N)

         call from_local_to_general_coord(
     $        interface_used,
     $        S,
     $        selected_grdpts_S)

         call from_local_to_general_coord(
     $        interface_used,
     $        E,
     $        selected_grdpts_E)

         call from_local_to_general_coord(
     $        interface_used,
     $        W,
     $        selected_grdpts_W)         


         !update the gridpoints with the selected
         !gridpoints asked to be turned from
         !bc_interior_pt to interior_pt
         call update_grdpts_id(
     $     interface_used, N, selected_grdpts_N)
         call update_grdpts_id(
     $     interface_used, S, selected_grdpts_S)
         call update_grdpts_id(
     $     interface_used, E, selected_grdpts_E)
         call update_grdpts_id(
     $     interface_used, W, selected_grdpts_W)


         !print interface after update
         call interface_used%print_binary(
     $        'nodes_2.dat',
     $        'grdpt_id_2.dat',
     $        'sizes_2.dat',
     $        '_2.dat')         


        contains

        subroutine determine_selected_grdpts(
     $       random_seed,
     $       mainlayer_id,
     $       relative_size, over_allocated,
     $       selected_grdpts)

          implicit none

          integer                             , intent(in)  :: random_seed
          integer                             , intent(in)  :: mainlayer_id
          integer                             , intent(in)  :: relative_size
          integer                             , intent(in)  :: over_allocated
          integer, dimension(:,:), allocatable, intent(out) :: selected_grdpts


          integer :: nb_grdpts,i,j
          

          allocate(selected_grdpts(2,4+relative_size))
          nb_grdpts = 0

          call srand(random_seed)

          select case(mainlayer_id)
            case(N)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              bc_size+1,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              bc_size+2,
     $              selected_grdpts, nb_grdpts)

               do i=1, relative_size
                  call add_new_selected_grdpt(
     $                 over_allocated+bc_size+i,
     $                 bc_size+2,
     $                 selected_grdpts, nb_grdpts)
               end do

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+relative_size+1,
     $              bc_size+2,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+relative_size+1,
     $              bc_size+1,
     $              selected_grdpts, nb_grdpts)

            case(S)
               do i=1, relative_size
                  call add_new_selected_grdpt(
     $                 over_allocated+bc_size+i,
     $                 over_allocated+bc_size,
     $                 selected_grdpts, nb_grdpts)
               end do

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+relative_size+1,
     $              over_allocated+bc_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+relative_size+1,
     $              over_allocated+bc_size+1,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              over_allocated+bc_size+1,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              over_allocated+bc_size+2,
     $              selected_grdpts, nb_grdpts)

            case(E)
               do j=1, relative_size
                  call add_new_selected_grdpt(
     $                 bc_size+2,
     $                 over_allocated+bc_size+j,
     $                 selected_grdpts, nb_grdpts)
               end do

               call add_new_selected_grdpt(
     $              bc_size+1,
     $              over_allocated+bc_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              bc_size+2,
     $              over_allocated+bc_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              bc_size+1,
     $              over_allocated+bc_size+relative_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              bc_size+2,
     $              over_allocated+bc_size+relative_size,
     $              selected_grdpts, nb_grdpts)

            case(W)
               do j=1, relative_size
                  call add_new_selected_grdpt(
     $                 over_allocated+bc_size,
     $                 over_allocated+bc_size+j,
     $                 selected_grdpts, nb_grdpts)
               end do

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              over_allocated+bc_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+1,
     $              over_allocated+bc_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size,
     $              over_allocated+bc_size+relative_size,
     $              selected_grdpts, nb_grdpts)

               call add_new_selected_grdpt(
     $              over_allocated+bc_size+1,
     $              over_allocated+bc_size+relative_size,
     $              selected_grdpts, nb_grdpts)

            case default
               print '(''test_bf_layer_update_grdpts_id_prog'')'
               print '(''determine_selected_grdpts'')'
               print '(''wrong mainlayer_id'')'
               
          end select

          call minimize_selected_grdpts(
     $         selected_grdpts, 1) !nb_grdpts

        end subroutine determine_selected_grdpts

      
        subroutine add_new_selected_grdpt(
     $     i,j, selected_grdpts, nb_grdpts)

          implicit none

          integer(ikind)         , intent(in)    :: i,j
          integer, dimension(:,:), intent(out)   :: selected_grdpts
          integer                , intent(inout) :: nb_grdpts

          real :: random_value

          random_value = RAND()

          if(random_value.ge.0.5) then

             nb_grdpts = nb_grdpts+1

             selected_grdpts(1,nb_grdpts) = i
             selected_grdpts(2,nb_grdpts) = j

          end if

        end subroutine add_new_selected_grdpt


        subroutine minimize_selected_grdpts(
     $     selected_grdpts, nb_grdpts)

          implicit none

          integer, dimension(:,:), allocatable, intent(inout) :: selected_grdpts
          integer                             , intent(in)    :: nb_grdpts


          integer, dimension(:,:), allocatable :: temp

          integer :: i,j

          allocate(temp(2,nb_grdpts))
          
          do j=1, nb_grdpts
             do i=1,2
                temp(i,j) = selected_grdpts(i,j)
             end do
          end do
          
          call MOVE_ALLOC(temp, selected_grdpts)

        end subroutine minimize_selected_grdpts


        subroutine identify_selected_grdpts(
     $     interface_used, mainlayer_id, selected_grdpts)

          implicit none

          class(bf_interface)    , intent(inout) :: interface_used
          integer                , intent(in)    :: mainlayer_id
          integer, dimension(:,:), intent(in)    :: selected_grdpts


          type(bf_mainlayer), pointer :: selected_mainlayer
          type(bf_sublayer) , pointer :: selected_sublayer
          integer(ikind)              :: i,j
          integer                     :: k
          integer(ikind), dimension(2):: sizes

          selected_mainlayer => interface_used%get_mainlayer(mainlayer_id)
          selected_sublayer  => selected_mainlayer%get_head_sublayer()
          sizes              =  selected_sublayer%get_sizes()

          do j=1, sizes(2)
             do i=1, sizes(1)
                call selected_sublayer%set_nodes_pt(
     $            i,j,1,background_grdpt)
             end do
          end do

          do k=1, size(selected_grdpts,2)
             i = selected_grdpts(1,k)
             j = selected_grdpts(2,k)
             call selected_sublayer%set_nodes_pt(
     $            i,j,1,identified_grdpt)
          end do

        end subroutine identify_selected_grdpts


        subroutine update_grdpts_id(
     $     interface_used, mainlayer_id, selected_grdpts)

          implicit none

          class(bf_interface)    , intent(inout) :: interface_used
          integer                , intent(in)    :: mainlayer_id
          integer, dimension(:,:), intent(in)    :: selected_grdpts


          type(bf_mainlayer), pointer :: selected_mainlayer
          type(bf_sublayer) , pointer :: selected_sublayer

          selected_mainlayer => interface_used%get_mainlayer(mainlayer_id)
          selected_sublayer  => selected_mainlayer%get_head_sublayer()

          call selected_sublayer%update_grdpts_id(selected_grdpts)

        end subroutine update_grdpts_id


        subroutine from_local_to_general_coord(
     $     interface_used, mainlayer_id, selected_grdpts)

          implicit none

          class(bf_interface)    , intent(inout) :: interface_used
          integer                , intent(in)    :: mainlayer_id
          integer, dimension(:,:), intent(inout) :: selected_grdpts


          type(bf_mainlayer) , pointer :: selected_mainlayer
          type(bf_sublayer)  , pointer :: selected_sublayer
          integer(ikind), dimension(2) :: match_table
          integer :: i,j

          selected_mainlayer => interface_used%get_mainlayer(mainlayer_id)
          selected_sublayer  => selected_mainlayer%get_head_sublayer()
          match_table        =  selected_sublayer%get_general_to_local_coord_tab()

          do j=1, size(selected_grdpts,2)
             do i=1,2
                selected_grdpts(i,j) = selected_grdpts(i,j) + match_table(i)
             end do
          end do

        end subroutine from_local_to_general_coord

      end program test_bf_layer_update_grdpts_id_prog
