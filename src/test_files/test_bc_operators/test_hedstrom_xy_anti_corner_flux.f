      !test the hedstrom_xy_anti_corner_flux_module by : 
      ! 1) checking thta the anti-corner is indeed the combination of
      !    a corner and edge type grid points
      ! 2) checking the symmetry
      program test_hedstrom_xy_anti_corner_flux

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated

        use hedstrom_xy_anti_corner_flux_module, only :
     $       compute_timedev_anti_corner_with_fluxes

        use hedstrom_xy_module, only :
     $       compute_timedev_x_edge_local,
     $       compute_timedev_y_edge_local,
     $       compute_timedev_corner_local

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_bf_layer, only :
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       
     $       no_overlap

        use parameters_constant, only :
     $       vector_x,
     $       vector_y,
     $       left, right

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.
        
        call check_inputs()

        test_loc = test_compute_timedev_anti_corner_with_fluxes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_timedev_anti_corner_with_fluxes: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated
        

        contains


        function test_compute_timedev_anti_corner_with_fluxes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(pmodel_eq)                  :: p_model
          real(rkind)                      :: t
          integer    , dimension(6,6)      :: grdpts_id
          real(rkind), dimension(6)        :: x_map
          real(rkind), dimension(6)        :: y_map
          real(rkind), dimension(6,6,ne)   :: nodes
          real(rkind), dimension(7,6,ne)   :: flux_x
          real(rkind), dimension(6,7,ne)   :: flux_y
          real(rkind), dimension(6,6,ne)   :: timedev
          type(sd_operators_x_oneside_L1)  :: s_x_L1
          type(sd_operators_x_oneside_R1)  :: s_x_R1
          type(sd_operators_y_oneside_L1)  :: s_y_L1
          type(sd_operators_y_oneside_R1)  :: s_y_R1
          real(rkind)                      :: dx
          real(rkind)                      :: dy
          integer    , dimension(5)        :: bc_section
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)   :: bf_alignment

          real(rkind), dimension(2,2,ne) :: timedev_ref
          integer(ikind) :: i,j
          integer        :: k

          test_validated = .true.


          call p_model%initial_conditions%ini_far_field()

          t=0.0d0

          grdpts_id = reshape((/
     $         1,1,1,1,1,1,
     $         1,1,1,1,1,1,
     $         1,1,2,2,2,2,
     $         1,1,2,3,3,3,
     $         1,1,2,3,0,0,
     $         1,1,2,3,0,0/),
     $         (/6,6/))
          
          x_map = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0, 5.5d0 ]
          y_map = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]
          nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.057d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.067d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $           0.02d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))

          nodes(5:6,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,2),k=1,ne)/),
     $         (/2,2,ne/))

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,7), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [NE_edge_type,3,3,no_overlap,no_overlap]


          !compute the NE_edge
          call compute_timedev_anti_corner_with_fluxes(
     $         t,
     $         bf_alignment,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         interior_nodes,
     $         s_x_L1, s_x_R1,
     $         s_y_L1, s_y_R1,
     $         p_model,
     $         bc_section,
     $         flux_x, flux_y,
     $         timedev)
          
          !check that only the time derivatives of the edge are
          !modified
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(NE): failed'')'
          end if

          !compute the reference time derivatives
          !compute the NE_corner at (3,3)
          i=3
          j=3
          timedev_ref(1,1,:) = compute_timedev_corner_local(
     $         t, x_map, y_map, nodes,
     $         p_model,
     $         gradient_x_x_oneside_R0, gradient_y_y_oneside_R0,
     $         dx,dy,
     $         i,j,
     $         right, right)


          !compute the N edge at (4,3)
          i=4
          j=3
          timedev_ref(2,1,:) = compute_timedev_y_edge_local(
     $         t, x_map, y_map, nodes,
     $         p_model,
     $         gradient_y_y_oneside_R0, dy,
     $         i,j,
     $         flux_x,dx,
     $         right)


          !compute the E edge at (3,4)
          i=3
          j=4
          timedev_ref(1,2,:) = compute_timedev_x_edge_local(
     $         t, x_map, y_map, nodes,
     $         p_model,
     $         gradient_x_x_oneside_R0, dx,
     $         i,j,
     $         flux_y, dy,
     $         right)


          !compute the NE_corner at (4,4)
          i=4
          j=4
          timedev_ref(2,2,:) = compute_timedev_corner_local(
     $         t, x_map, y_map, nodes,
     $         p_model,
     $         gradient_x_x_oneside_R0, gradient_y_y_oneside_R0,
     $         dx,dy,
     $         i,j,
     $         right, right)


          !compare both
          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(NE) failed'')'
          end if

          !compute the NW_corner
          call reflect_x(
     $         p_model,
     $         x_map,
     $         grdpts_id,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,7), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [NW_edge_type,3,3,no_overlap,no_overlap]

          call compute_timedev_anti_corner_with_fluxes(
     $         t,
     $         bf_alignment,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         interior_nodes,
     $         s_x_L1, s_x_R1,
     $         s_y_L1, s_y_R1,
     $         p_model,
     $         bc_section,
     $         flux_x, flux_y,
     $         timedev)
          
          !check that only the time derivatives of the edge are
          !modified
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(NW): failed'')'
          end if

          !compare with the reference
          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(NW) failed'')'
          end if


          !compute the SW_corner
          call reflect_y(
     $         p_model,
     $         y_map,
     $         grdpts_id,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,7), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [SW_edge_type,3,3,no_overlap,no_overlap]

          call compute_timedev_anti_corner_with_fluxes(
     $         t,
     $         bf_alignment,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         interior_nodes,
     $         s_x_L1, s_x_R1,
     $         s_y_L1, s_y_R1,
     $         p_model,
     $         bc_section,
     $         flux_x, flux_y,
     $         timedev)
          
          !check that only the time derivatives of the edge are
          !modified
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(SW): failed'')'
          end if

          !compare with the reference
          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(SW) failed'')'
          end if


          !compute the SE_corner
          call reflect_x(
     $         p_model,
     $         x_map,
     $         grdpts_id,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,7), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [SE_edge_type,3,3,no_overlap,no_overlap]

          call compute_timedev_anti_corner_with_fluxes(
     $         t,
     $         bf_alignment,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         interior_nodes,
     $         s_x_L1, s_x_R1,
     $         s_y_L1, s_y_R1,
     $         p_model,
     $         bc_section,
     $         flux_x, flux_y,
     $         timedev)
          
          !check that only the time derivatives of the edge are
          !modified
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(SE): failed'')'
          end if

          !compare with the reference
          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(SE) failed'')'
          end if

        end function test_compute_timedev_anti_corner_with_fluxes


        function check_around_timedev(
     $     timedev,
     $     gen_coords,
     $     cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: timedev
          integer(ikind), dimension(2,2)  , intent(in) :: gen_coords
          real(rkind)                     , intent(in) :: cst
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          integer(ikind) :: i,j
          integer        :: k

          
          test_validated = .true.


          do k=1, size(timedev,3)
             do j=1, size(timedev,2)
                do i=1, size(timedev,1)

                   test_loc = is_real_validated(timedev(i,j,k),cst,.false.)

                   if(.not.test_loc) then
                      test_loc =
     $                     (i.ge.gen_coords(1,1)).and.
     $                     (i.le.gen_coords(1,2)).and.
     $                     (j.ge.gen_coords(2,1)).and.
     $                     (j.le.gen_coords(2,2))
                      
                   end if

                   if(detailled.and.(.not.test_loc)) then
                      print *, timedev(i,j,k)
                      print *, cst
                      print '(''['',3I4,'']'')', i,j,k
                   end if

                   test_validated = test_validated.and.test_loc

                end do
             end do
          end do

        end function check_around_timedev


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: far_field

          call p_model%initial_conditions%ini_far_field()
          far_field = p_model%initial_conditions%get_far_field(0.0d0,1.0d0,1.0d0)

          if(.not.is_real_vector_validated(
     $         far_field,
     $         [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $         .true.)) then
             print '(''the test requires p_model%get_far_field(t,x,y)='')'
             print '(''[1.465102139d0,0.14651021d0,0.0d0,2.84673289d0]'')'
             print '()'
             print '(''T0 should be 0.95'')'
             print '(''flow_direction should be x-direction'')'
             print '(''ic_choice should be newgrdpt_test'')'
             print '()'
             stop ''
             
          end if

        end subroutine check_inputs


        subroutine reflect_x(
     $     p_model,
     $     bf_x_map,
     $     grdpts_id,
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_x_map
          integer    , dimension(6,6)   , intent(inout) :: grdpts_id
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: i_r
          integer :: j
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_x_map_r
          integer    , dimension(6,6)    :: grdpts_id_r
          real(rkind), dimension(6,6,ne) :: bf_nodes_r
          real(rkind), dimension(2,2,ne) :: timedev_ref_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_x) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the x_maps
          do i=1, size(bf_x_map,1)

             i_r = size(bf_x_map)-(i-1)
             bf_x_map_r(i) = -bf_x_map(i_r)

          end do
          bf_x_map = bf_x_map_r


          !reflect the grdpts_id
          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)
                i_r = size(grdpts_id,1)-(i-1)
                grdpts_id_r(i,j) = grdpts_id(i_r,j)
             end do
          end do
          grdpts_id = grdpts_id_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_x) then

                do j=1, size(bf_nodes,2)
                   do i=1, size(bf_nodes,1)
                      i_r = size(bf_nodes,1)-(i-1)
                      bf_nodes_r(i,j,k) = - bf_nodes(i_r,j,k)
                   end do
                end do

             else

                do j=1, size(bf_nodes,2)
                   do i=1, size(bf_nodes,1)
                      i_r = size(bf_nodes,1)-(i-1)
                      bf_nodes_r(i,j,k) = bf_nodes(i_r,j,k)
                   end do
                end do

             end if             

          end do
          bf_nodes = bf_nodes_r


          !reflect the timedev_ref
          do k=1, ne
             if(var_type(k).eq.vector_x) then
                do j=1, size(timedev_ref,2)
                   do i=1, size(timedev_ref,1)
                      i_r = size(timedev_ref,1)-(i-1)
                      timedev_ref_r(i,j,k) = - timedev_ref(i_r,j,k)
                   end do
                end do
             else
                do j=1, size(timedev_ref,2)
                   do i=1, size(timedev_ref,1)
                      i_r = size(timedev_ref,1)-(i-1)
                      timedev_ref_r(i,j,k) = timedev_ref(i_r,j,k)
                   end do
                end do
             end if
          end do
          timedev_ref = timedev_ref_r

        end subroutine reflect_x


        subroutine reflect_y(
     $     p_model,
     $     bf_y_map,
     $     grdpts_id,
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_y_map
          integer    , dimension(6,6)   , intent(inout) :: grdpts_id
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: j
          integer :: j_r
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_y_map_r
          real(rkind), dimension(6,6)    :: grdpts_id_r
          real(rkind), dimension(6,6,ne) :: bf_nodes_r
          real(rkind), dimension(2,2,ne) :: timedev_ref_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_y) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the y_maps
          do j=1, size(bf_y_map,1)

             j_r = size(bf_y_map,1)-(j-1)
             bf_y_map_r(j) = -bf_y_map(j_r)

          end do
          bf_y_map = bf_y_map_r


          !reflect the grdpts_id
          do j=1, size(grdpts_id,2)
             j_r = size(grdpts_id,2)-(j-1)
             do i=1, size(grdpts_id,1)                
                grdpts_id_r(i,j) = grdpts_id(i,j_r)
             end do
          end do
          grdpts_id = grdpts_id_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_y) then

                do j=1, size(bf_nodes,2)

                   j_r = size(bf_nodes,2)-(j-1)

                   do i=1, size(bf_nodes,1)
                      bf_nodes_r(i,j,k) = - bf_nodes(i,j_r,k)
                   end do

                end do

             else

                do j=1, size(bf_nodes,2)

                   j_r = size(bf_nodes,2)-(j-1)

                   do i=1, size(bf_nodes,1)
                      bf_nodes_r(i,j,k) = bf_nodes(i,j_r,k)
                   end do

                end do

             end if             

          end do
          bf_nodes = bf_nodes_r


          !reflect the timedev_ref
          do k=1, ne

             if(var_type(k).eq.vector_y) then 

                do j=1, size(timedev_ref,2)

                   j_r = size(timedev_ref,2)-(j-1)

                   do i=1, size(timedev_ref,1)
                      timedev_ref_r(i,j,k) = - timedev_ref(i,j_r,k)
                   end do
                end do

             else
                
                do j=1, size(timedev_ref,2)

                   j_r = size(timedev_ref,2)-(j-1)

                   do i=1, size(timedev_ref,1)
                      timedev_ref_r(i,j,k) = timedev_ref(i,j_r,k)
                   end do
                end do

             end if
          end do
          timedev_ref = timedev_ref_r

        end subroutine reflect_y

      end program test_hedstrom_xy_anti_corner_flux
