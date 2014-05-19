      program test_bf_layer_adapt_corner_prog

        use bf_layer_adapt_corner_module, only : adapt_bf_layer_N_to_corner
        use bf_interface_class          , only : bf_interface
        use bf_sublayer_class           , only : bf_sublayer

        use parameters_constant         , only : N_E,N_W,S_E,S_W
        use parameters_input            , only : nx,ny,ne
        use parameters_kind             , only : ikind, rkind

        use test_bf_layer_module        , only : print_interior_data,
     $                                           ini_grdpts_id,
     $                                           ini_nodes
        use test_cases_interface_module , only : ini_interface

        implicit none


        integer, parameter               :: corner_tested = 1
        integer, parameter               :: corner_order = 1
        integer, parameter               :: test_case_id = 21
        integer, parameter               :: bf_corner_distance = 6
        integer, parameter               :: over_allocated = 0

        integer, dimension(4)            :: corner_table
        integer                          :: corner_id

        real(rkind), dimension(nx,ny,ne) :: nodes
        integer    , dimension(nx,ny)    :: grdpts_id
        type(bf_interface)               :: interface_used

        type(bf_sublayer), pointer       :: sublayer1
        type(bf_sublayer), pointer       :: sublayer2
        integer                          :: i

        !choose the corner tested
        corner_table       = [N_E,N_W,S_E,S_W]
        corner_id          = corner_table(corner_tested)

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
     $       corner_id,
     $       corner_order,
     $       bf_corner_distance,
     $       over_allocated)

        !print interface before
        call interface_used%print_binary(
     $       'nodes_0.dat',
     $       'grdpt_id_0.dat',
     $       'sizes_0.dat',
     $       '_0.dat')

        !adapt the buffer layers neighboring the
        !corner to the new corner buffer layer
        !1) get the neighboring buffer layer of the corner
        !2) for both buffer layers, check if they are compatible
        !3) if they are compatible, adapt them
        do i=1, size(corner_table,1)

           print '(''corner_id: ''  , I2)', corner_table(i)

           !1) get the neighboring sublayers to the corner
           call interface_used%get_neighboring_sublayers(
     $          corner_table(i),
     $          sublayer1, sublayer2)

           !2+3) for each buffer layer, check whether they are
           !compatible with the new corner buffer layer
           call test_sublayer_for_corner(corner_table(i), nodes, sublayer1)
           call test_sublayer_for_corner(corner_table(i), nodes, sublayer2)

           print '(''*********************'')'
           print '('''')'

        end do

        !print interface after
        call interface_used%print_binary(
     $       'nodes_1.dat',
     $       'grdpt_id_1.dat',
     $       'sizes_1.dat',
     $       '_1.dat')

        contains

        subroutine test_sublayer_for_corner(
     $       corner_id, nodes, sublayer_ptr)

          implicit none

          integer                         , intent(in)    :: corner_id
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          type(bf_sublayer), pointer      , intent(inout) :: sublayer_ptr

          logical                          :: compatible
          integer(ikind), dimension(2,2)   :: new_alignment
          logical       , dimension(4)     :: new_neighbors

          if(associated(sublayer_ptr)) then
             
             compatible = sublayer_ptr%is_compatible_with_corner(
     $            corner_id,
     $            new_alignment,
     $            new_neighbors)
             
             print '(''sublayer_id: '', I2)',
     $            sublayer_ptr%get_localization()
             
             if(compatible) then
                print '(''sublayer compatible'')'
                call sublayer_ptr%adapt_to_corner(
     $               nodes, new_alignment, new_neighbors)
             else
                print '(''sublayer incompatible'')'
             end if
             
           else
              print '(''no neighboring layer'')'
           end if
           print '('''')'
           
        end subroutine test_sublayer_for_corner

      end program test_bf_layer_adapt_corner_prog
