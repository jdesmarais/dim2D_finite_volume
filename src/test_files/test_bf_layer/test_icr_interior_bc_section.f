      program test_icr_interior_bc_section

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use check_data_module, only :
     $       is_int_matrix_validated,
     $       is_int_matrix3D_validated

        use icr_interface_class, only :
     $       icr_interface

        use icr_interior_bc_section_module, only :
     $       analyze_interior_square_xy,
     $       analyze_interior_square_bounds,
     $       analyze_interior_edge_x,
     $       analyze_interior_edge_y

        use icr_path_chain_class, only :
     $       icr_path_chain

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       
     $       dct_icr_distance,
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       cpt2not_and_cpt3normal,
     $       cpt2not_and_cpt3not

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_analyze_interior_square_xy(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_interior_square_xy: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_interior_square_bounds(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_interior_square_bounds: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_interior_edge_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_interior_edge_x: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_interior_edge_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_interior_edge_y: '',L1)', test_loc
        print '()'


        contains

        function test_analyze_interior_square_xy(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated 


          type(bf_interface_coords)        :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          integer                        :: k
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used


          test_validated = .true.


          do k=1,11

             !input
             call get_test_param_analyze_square_xy(
     $            k,
     $            bf_interface_used,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            bc_section,
     $            test_nb_pts,
     $            test_pts)

             !output
             call analyze_interior_square_xy(
     $            bf_interface_used,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model,
     $            bc_section)

             !validation
             if(test_nb_pts.gt.0) then
                test_loc = icr_interface_used%paths%get_nb_paths().eq.1
             else
                test_loc = icr_interface_used%paths%get_nb_paths().eq.0
             end if
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_nb_pts('',I2,'') failed '')',k
                print '(''nb_paths: '',I2)', icr_interface_used%paths%get_nb_paths()
             end if

             if(test_loc.and.(test_nb_pts.gt.0)) then
                icr_path_used => icr_interface_used%paths%get_head_path()
                test_loc = is_int_matrix_validated(
     $               icr_path_used%pts(:,1:test_nb_pts),
     $               test_pts(:,1:test_nb_pts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_pts('',I2,'') failed'')',k
                end if
             end if

             icr_path_used => icr_interface_used%paths%get_head_path()
             if(associated(icr_path_used)) then
                call icr_interface_used%paths%remove_path(icr_path_used)
             end if

           end do

        end function test_analyze_interior_square_xy


        subroutine get_test_param_analyze_square_xy(
     $     test_id,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          type(bf_interface_coords)       , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          integer(ikind)                 :: i,j,k
          integer(ikind), dimension(2,3) :: gen_coords
          real(rkind)   , dimension(3)   :: nodes_set
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((D,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          bc_section = [NW_corner_type,1,ny-1,no_overlap,no_overlap]

          gen_coords(:,1) = [ 2+dct_icr_distance, ny-1-dct_icr_distance-1]
          gen_coords(:,2) = [ 2+dct_icr_distance, ny-1-dct_icr_distance  ]
          gen_coords(:,3) = [ 3+dct_icr_distance, ny-1-dct_icr_distance  ]
          

          select case(test_id)

            !all pts activated
            case(1)
               nodes_set = [A,A,A]

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all pts desactivated
            case(2)
               nodes_set = [D,D,D]
               
               test_nb_pts = 0

            !only one point activated
            case(3)
               nodes_set = [A,D,D]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !only one point activated
            case(4)
               nodes_set = [D,A,D]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !only one point activated
            case(5)
               nodes_set = [D,D,A]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all activated in spite of the overlap
            case(6)
               nodes_set = [A,A,A]
               bc_section(5) = N_overlap

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all activated in spite of the overlap
            case(7)
               nodes_set = [A,A,A]
               bc_section(5) = W_overlap

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !no point activated because of the overlap
            case(8)
               nodes_set = [A,A,A]
               bc_section(5) = E_overlap

               test_nb_pts = 0

            !no point activated because of the overlap
            case(9)
               nodes_set = [A,A,A]
               bc_section(5) = S_overlap

               test_nb_pts = 0

            !no point activated because of the overlap
            case(10)
               nodes_set = [A,A,A]
               bc_section(4) = cpt2not_and_cpt3normal

               test_nb_pts = 0

            !no point activated because of the overlap
            case(11)
               nodes_set = [A,A,A]
               bc_section(4) = cpt2not_and_cpt3not

               test_nb_pts = 0

          end select

          do k=1,3
             interior_nodes(gen_coords(1,k),gen_coords(2,k),1) = nodes_set(k)
          end do


        end subroutine get_test_param_analyze_square_xy


        function test_analyze_interior_square_bounds(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer(ikind), dimension(5)     :: bc_section

          integer(ikind), dimension(2,2,2) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: activated_grdpts
          integer                          :: nb_activated_grdpts
          logical                          :: no_activation

          integer(ikind), dimension(2,2,2) :: test_analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: test_activated_grdpts
          integer                          :: test_nb_activated_grdpts
          logical                          :: test_no_activation

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,8

             !input
             call get_test_param_square_bounds(
     $            k,
     $            bc_section,
     $            test_analyzed_grdpts_bounds,
     $            test_activated_grdpts,
     $            test_nb_activated_grdpts,
     $            test_no_activation)

             !output
             call analyze_interior_square_bounds(
     $            bc_section,
     $            analyzed_grdpts_bounds,
     $            activated_grdpts,
     $            nb_activated_grdpts,
     $            no_activation)

             !validation
             !no_activation
             test_loc = no_activation.eqv.test_no_activation
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test no_activation failed'')'
             end if

             if(test_loc.and.(.not.no_activation)) then

                !analyze_grdpts_bounds
                test_loc = is_int_matrix3D_validated(
     $               analyzed_grdpts_bounds,
     $               test_analyzed_grdpts_bounds,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test analyze_grdpts_bounds failed'')'
                end if

                !nb_activated_grdpts
                test_loc = nb_activated_grdpts.eq.test_nb_activated_grdpts
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test nb_activated_grdpts failed'')'
                end if

                !activated_grdpts
                test_loc = is_int_matrix_validated(
     $               activated_grdpts(:,1:test_nb_activated_grdpts),
     $               test_activated_grdpts(:,1:test_nb_activated_grdpts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test activated_grdpts failed'')'
                end if

             end if

          end do

        end function test_analyze_interior_square_bounds


        subroutine get_test_param_square_bounds(
     $     test_id,
     $     bc_section,
     $     test_analyzed_grdpts_bounds,
     $     test_activated_grdpts,
     $     test_nb_activated_grdpts,
     $     test_no_activation)

          implicit none

          integer                         , intent(in)  :: test_id
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer(ikind), dimension(2,2,2), intent(out) :: test_analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(out) :: test_activated_grdpts
          integer                         , intent(out) :: test_nb_activated_grdpts
          logical                         , intent(out) :: test_no_activation

          
          bc_section(4) = no_overlap
          bc_section(5) = no_overlap

          test_no_activation = .false.


          select case(test_id)
            case(1)
               bc_section(1) = NW_corner_type
               bc_section(2) = 1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              ny-1-dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [2,ny-1]


            case(2)
               bc_section(1) = NE_corner_type
               bc_section(2) = nx-1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              nx-1-dct_icr_distance-1,
     $              ny-1-dct_icr_distance,
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [nx-1,ny-1]

            case(3)
               bc_section(1) = SW_corner_type
               bc_section(2) = 1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [2,2]

            case(4)
               bc_section(1) = SE_corner_type
               bc_section(2) = nx-1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [nx-1,2]

            case(5)
               bc_section(1) = NW_edge_type
               bc_section(2) = 1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              0,
     $              ny-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              ny+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [1,ny-1]
               test_activated_grdpts(:,2) = [2,ny-1]
               test_activated_grdpts(:,3) = [2,ny  ]
               

            case(6)
               bc_section(1) = NE_edge_type
               bc_section(2) = nx-1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              nx+1,
     $              ny-1-dct_icr_distance,
     $              
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance+1,
     $              nx-1-dct_icr_distance,
     $              ny+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [nx-1,ny-1]
               test_activated_grdpts(:,2) = [nx  ,ny-1]
               test_activated_grdpts(:,3) = [nx-1,ny  ]

            case(7)
               bc_section(1) = SW_edge_type
               bc_section(2) = 1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              0,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance-1,
     $              
     $              0,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [2,1]
               test_activated_grdpts(:,2) = [1,2]
               test_activated_grdpts(:,3) = [2,2]

            case(8)
               bc_section(1) = SE_edge_type
               bc_section(2) = nx-1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              0,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance-1,
     $              
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              nx+1,
     $              2+dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [nx-1,1]
               test_activated_grdpts(:,2) = [nx-1,2]
               test_activated_grdpts(:,3) = [nx  ,2]

          end select

        end subroutine get_test_param_square_bounds


        function test_analyze_interior_edge_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_coords)        :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          integer                        :: k,l
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used


          test_validated = .true.


          do l=1,2
             do k=1,8

                !input
                call get_test_param_analyze_interior_edge_x(
     $               k,l,
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               bc_section,
     $               test_nb_pts,
     $               test_pts)

                !output
                call analyze_interior_edge_x(
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model,
     $               bc_section)

                !validation
                if(test_nb_pts.gt.0) then
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.1
                else
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.0
                end if
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_nb_pts('',2I2,'') failed '')',k,l
                end if

                if(test_loc.and.(test_nb_pts.gt.0)) then
                   icr_path_used => icr_interface_used%paths%get_head_path()
                   test_loc = is_int_matrix_validated(
     $                  icr_path_used%pts(:,1:test_nb_pts),
     $                  test_pts(:,1:test_nb_pts),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''test_pts('',2I2,'') failed'')',k,l
                   end if
                end if

                icr_path_used => icr_interface_used%paths%get_head_path()
                if(associated(icr_path_used)) then
                   call icr_interface_used%paths%remove_path(icr_path_used)
                end if

             end do
          end do

        end function test_analyze_interior_edge_x


        subroutine get_test_param_analyze_interior_edge_x(
     $     test_id,
     $     edge_type,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          integer                         , intent(in)  :: edge_type
          type(bf_interface_coords)       , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          real(rkind), dimension(ne) :: unactivated_node
          real(rkind), dimension(ne) :: activated_node

          integer(ikind) :: i,j,k
          integer(ikind) :: j_s
          integer :: i_pt
          integer :: i_activated
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((0.0d0,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          j_s = ny/2

          select case(edge_type)
            case(1)
               i_pt = nx-1
               bc_section = [E_edge_type,i_pt  ,j_s,j_s+8,no_overlap]
               i_activated = i_pt-dct_icr_distance
            case(2)
               i_pt = 2
               bc_section = [W_edge_type,i_pt-1,j_s,j_s+8,no_overlap]
               i_activated = i_pt+dct_icr_distance
          end select

          unactivated_node = [ D,0.0d0,0.0d0]
          activated_node   = [ A,0.0d0,0.0d0]
          

          select case(test_id)

            !all pts activated
            case(1)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !all pts desactivated
            case(2)
               interior_nodes(i_activated,j_s:j_s+8,1) = [D,D,D,D,D,D,D,D,D]
               
               test_nb_pts = 0

            !one out of two points activated
            case(3)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,A,D,A,D,A,D,A]
               
               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !one out of three points activated
            case(4)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,D,A,D,D,A,D,D]
               
               test_nb_pts = 8
               test_pts(1,:)   = i_pt
               test_pts(2,1:8) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                            j_s+5,j_s+6,j_s+7]
               

            !one out of four points activated
            case(5)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,D,D,A,D,D,D,A]
               
               test_nb_pts = 7

               test_pts(1,:)   = i_pt
               test_pts(2,1:7) = [j_s,j_s+1,j_s+3,j_s+4,
     $                            j_s+5,j_s+7,j_s+8]

            !all activated in spite of the overlap
            case(6)
               select case(edge_type)
                 case(1)
                    bc_section(5) = E_overlap
                 case(2)
                    bc_section(5) = W_overlap
               end select

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !no point activated because of the overlap
            case(7)
               select case(edge_type)
                 case(1)
                    bc_section(5) = W_overlap
                 case(2)
                    bc_section(5) = E_overlap
               end select

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

            !no point activated because of the overlap
            case(8)
               bc_section(5) = EW_overlap

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

          end select               

        end subroutine get_test_param_analyze_interior_edge_x


        function test_analyze_interior_edge_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_coords)        :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          integer                        :: k,l
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used


          test_validated = .true.


          do l=1,2
             do k=1,8

                !input
                call get_test_param_analyze_interior_edge_y(
     $               k,l,
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               bc_section,
     $               test_nb_pts,
     $               test_pts)

                !output
                call analyze_interior_edge_y(
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model,
     $               bc_section)

                !validation
                if(test_nb_pts.gt.0) then
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.1
                else
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.0
                end if
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_nb_pts('',2I2,'') failed '')',k,l
                end if

                if(test_loc.and.(test_nb_pts.gt.0)) then
                   icr_path_used => icr_interface_used%paths%get_head_path()
                   test_loc = is_int_matrix_validated(
     $                  icr_path_used%pts(:,1:test_nb_pts),
     $                  test_pts(:,1:test_nb_pts),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''test_pts('',2I2,'') failed'')',k,l
                   end if
                end if

                icr_path_used => icr_interface_used%paths%get_head_path()
                if(associated(icr_path_used)) then
                   call icr_interface_used%paths%remove_path(icr_path_used)
                end if

             end do
          end do

        end function test_analyze_interior_edge_y


        subroutine get_test_param_analyze_interior_edge_y(
     $     test_id,
     $     edge_type,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          integer                         , intent(in)  :: edge_type
          type(bf_interface_coords)       , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          real(rkind), dimension(ne) :: unactivated_node
          real(rkind), dimension(ne) :: activated_node

          integer(ikind) :: i,j,k
          integer(ikind) :: i_s
          integer :: j_pt
          integer :: j_activated
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((0.0d0,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          i_s = nx/2

          select case(edge_type)
            case(1)
               j_pt = ny-1
               bc_section = [N_edge_type,i_s,j_pt,i_s+8,no_overlap]
               j_activated = j_pt-dct_icr_distance
            case(2)
               j_pt = 2
               bc_section = [S_edge_type,i_s,j_pt-1,i_s+8,no_overlap]
               j_activated = j_pt+dct_icr_distance
          end select

          unactivated_node = [ D,0.0d0,0.0d0]
          activated_node   = [ A,0.0d0,0.0d0]
          

          select case(test_id)

            !all pts activated
            case(1)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9
               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !all pts desactivated
            case(2)
               interior_nodes(i_s:i_s+8,j_activated,1) = [D,D,D,D,D,D,D,D,D]
               
               test_nb_pts = 0

            !one out of two points activated
            case(3)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,A,D,A,D,A,D,A]
               
               test_nb_pts = 9
               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !one out of three points activated
            case(4)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,D,A,D,D,A,D,D]
               
               test_nb_pts = 8
               test_pts(1,1:8) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7]
               test_pts(2,:) = j_pt

            !one out of four points activated
            case(5)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,D,D,A,D,D,D,A]
               
               test_nb_pts = 7
               test_pts(1,1:7) = [i_s,i_s+1,i_s+3,i_s+4,
     $                            i_s+5,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !all activated in spite of the overlap
            case(6)
               select case(edge_type)
                 case(1)
                    bc_section(5) = N_overlap
                 case(2)
                    bc_section(5) = S_overlap
               end select

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !no point activated because of the overlap
            case(7)
               select case(edge_type)
                 case(1)
                    bc_section(5) = S_overlap
                 case(2)
                    bc_section(5) = N_overlap
               end select

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

           !no point activated because of the overlap
            case(8)
               bc_section(5) = NS_overlap

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

          end select               

        end subroutine get_test_param_analyze_interior_edge_y


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)                :: p_model
          real(rkind), dimension(3)      :: x_map
          real(rkind), dimension(3)      :: y_map
          real(rkind), dimension(3,3,ne) :: nodes
          
          logical :: test_loc


          nodes(2,2,1) = -1.0d0
          test_loc = p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''openbc_undermined if nodes(2,2,1)<0'')'
             stop ''
          end if

          
          nodes(2,2,1) = 1.0d0
          test_loc = .not.p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''.not.openbc_undermined if nodes(2,2,1)>0'')'
             stop ''
          end if

        end subroutine check_inputs

      end program test_icr_interior_bc_section
