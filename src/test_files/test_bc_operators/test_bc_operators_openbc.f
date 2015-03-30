      program test_bc_operators_openbc

        use bc_operators_class, only :
     $       bc_operators

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated
        
        use hedstrom_xy_module, only :
     $       compute_timedev_corner_local

        use parameters_bf_layer, only :
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       no_overlap

        use parameters_constant, only :
     $       vector_x,
     $       vector_y,
     $       left,
     $       right

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_R0

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_compute_timedev_corner(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_timedev_corner: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated


        contains


        function test_compute_timedev_corner(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          type(bc_operators)             :: bc_operators_openbc_used
          type(pmodel_eq)                :: p_model
          real(rkind)                    :: t
          real(rkind), dimension(6,6,ne) :: nodes
          real(rkind), dimension(6)      :: x_map
          real(rkind), dimension(6)      :: y_map
          real(rkind)                    :: dx
          real(rkind)                    :: dy
          integer    , dimension(5)      :: bc_section
          real(rkind), dimension(6,6,ne) :: timedev

          real(rkind), dimension(2,2,ne) :: timedev_ref

          integer(ikind) :: i,j
          integer        :: k
          logical        :: test_loc


          test_validated = .true.

          !NE_corner
          !============================================================
          !input
          
          call p_model%initial_conditions%ini_far_field()

          t=0.0d0
          
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

          nodes(5:6,1:4,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,4),k=1,ne)/),
     $         (/2,4,ne/))

          nodes(1:6,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,6),j=1,2),k=1,ne)/),
     $         (/6,2,ne/))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [NE_corner_type,3,3,no_overlap,no_overlap]

          !output
          call bc_operators_openbc_used%compute_timedev_corner(
     $         t,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bc_section,
     $         timedev)

          do j=3,4
             do i=3,4
                timedev_ref(i-2,j-2,:) = compute_timedev_corner_local(
     $               t,
     $               x_map,
     $               y_map,
     $               nodes,
     $               p_model,
     $               gradient_x_x_oneside_R0, gradient_y_y_oneside_R0,
     $               dx,dy,
     $               i,j,
     $               right, right)

             end do
          end do

          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(NE): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(NE) failed'')'
          end if


          !NW_corner
          !============================================================
          !input
          call reflect_x(
     $         p_model,
     $         x_map,
     $         nodes,
     $         timedev_ref)

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [NW_corner_type,3,3,no_overlap,no_overlap]

          !output
          call bc_operators_openbc_used%compute_timedev_corner(
     $         t,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bc_section,
     $         timedev)

          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(NW): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(NW) failed'')'
          end if


          !SW_corner
          !============================================================
          !input
          call reflect_y(
     $         p_model,
     $         y_map,
     $         nodes,
     $         timedev_ref)

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [SW_corner_type,3,3,no_overlap,no_overlap]

          !output
          call bc_operators_openbc_used%compute_timedev_corner(
     $         t,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bc_section,
     $         timedev)

          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(SW): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(SW) failed'')'
          end if


          !SE_corner
          !============================================================
          !input
          call reflect_x(
     $         p_model,
     $         x_map,
     $         nodes,
     $         timedev_ref)

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          bc_section = [SE_corner_type,3,3,no_overlap,no_overlap]

          !output
          call bc_operators_openbc_used%compute_timedev_corner(
     $         t,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bc_section,
     $         timedev)

          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,3,4,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(SE): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(SE) failed'')'
          end if

        end function test_compute_timedev_corner


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
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_x_map
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: i_r
          integer :: j
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_x_map_r
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
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_y_map
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: j
          integer :: j_r
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_y_map_r
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

      end program test_bc_operators_openbc
