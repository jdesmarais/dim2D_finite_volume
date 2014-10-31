      !the test consists of integrating two identical
      !initial computational domains and to check that
      !the time integration leads to the same results
      !
      !the only difference between the computational
      !fields is that one field is computed using one
      !unique domain while the second computational
      !field is made of four parts
      !
      !   ________         ________
      !  |        |       |________|
      !  |        |   ?   | |    | |
      !  |        |   =   |_|____|_|
      !  |________|       |________|
      !
      !----------------------------------------------
      program test_field_extended_time_dev

        use bf_interface_class, only :
     $       bf_interface

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use field_class, only :
     $       field

        use field_extended_class, only :
     $       field_extended

        use ic_class, only :
     $       ic

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W,
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       x_min, x_max,
     $       y_min, y_max

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use rk3tvd_steps_module, only :
     $       compute_1st_step, compute_1st_step_nopt

        implicit none

        
        type(field)            :: field_tested
        type(field_extended)   :: field_ext_tested

        real(rkind), parameter :: dt = 1.0

        logical    , parameter :: one_piece   = .true.
        logical    , parameter :: four_pieces = .false.
        logical                :: domain_decomposition

        real(rkind), dimension(nx)       :: x_map
        real(rkind), dimension(ny)       :: y_map
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind), dimension(nx,ny,ne) :: nodes_tmp
        real(rkind), dimension(nx,ny,ne) :: time_dev

        domain_decomposition = four_pieces


        !one domain
        if(domain_decomposition.eqv.one_piece) then

           if(.not.((nx.eq.40).and.(ny.eq.40))) then
              stop 'the one-piece test requires (nx,ny)=(40,40)'
           end if


           !0.0) initialization
           call field_tested%ini()
           call initialize_maps_one_piece(x_map,y_map)
           call initialize_nodes(x_map,y_map,nodes)
           
           call field_tested%set_x_map(x_map)
           call field_tested%set_y_map(y_map)
           call field_tested%set_nodes(nodes)

           call write_data(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_0.dat',
     $          'y_map_0.dat',
     $          'nodes_0.dat',
     $          'sizes_0.dat')


           !1) time derivatives computation
           time_dev = field_tested%compute_time_dev()

           call write_data(
     $          x_map,
     $          y_map,
     $          time_dev,
     $          'x_map_1.dat',
     $          'y_map_1.dat',
     $          'nodes_1.dat',
     $          'sizes_1.dat')


           !2) first step integration computation
           call field_tested%compute_integration_step(
     $          dt, nodes_tmp, time_dev, compute_1st_step)

           x_map = field_tested%get_x_map()
           y_map = field_tested%get_y_map()
           nodes = field_tested%get_nodes()

           call write_data(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_2.dat',
     $          'y_map_2.dat',
     $          'nodes_2.dat',
     $          'sizes_2.dat')


           !3) nodes synchronization
           call field_tested%apply_bc_on_nodes()

           x_map = field_tested%get_x_map()
           y_map = field_tested%get_y_map()
           nodes = field_tested%get_nodes()

           call write_data(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_3.dat',
     $          'y_map_3.dat',
     $          'nodes_3.dat',
     $          'sizes_3.dat')


           !4) integration cycle
           call field_tested%integrate(dt)

           x_map = field_tested%get_x_map()
           y_map = field_tested%get_y_map()
           nodes = field_tested%get_nodes()

           call write_data(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_4.dat',
     $          'y_map_4.dat',
     $          'nodes_4.dat',
     $          'sizes_4.dat')

        end if


        if(domain_decomposition.eqv.four_pieces) then

           if(.not.((nx.eq.24).and.(ny.eq.24))) then
              stop 'the one-piece test requires (nx,ny)=(20,20)'
           end if

           !0) initialization
           call field_ext_tested%ini()

           call initialize_four_pieces(
     $          x_map,
     $          y_map,
     $          nodes,
     $          field_ext_tested)

           call write_data_ext(
     $          field_ext_tested,
     $          'x_map_ext_0.dat',
     $          'y_map_ext_0.dat',
     $          'nodes_ext_0.dat',
     $          'sizes_ext_0.dat')

           x_map = field_ext_tested%get_x_map()
           y_map = field_ext_tested%get_y_map()
           nodes = field_ext_tested%get_nodes()

           call write_grdpts_id(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_ext_int_0.dat',
     $          'y_map_ext_int_0.dat',
     $          'grdptsid_ext_int_0.dat',
     $          'nodes_ext_int_0.dat',
     $          'sizes_ext_int_0.dat')
           
           call print_interface(
     $          field_ext_tested%domain_extension,
     $          0)


           !1) time derivatives computation
           call field_ext_tested%domain_extension%allocate_before_timeInt()

           time_dev = field_ext_tested%compute_time_dev_ext()

           call write_timedev_ext(
     $          field_ext_tested,
     $          time_dev,
     $          'x_map_ext_1.dat',
     $          'y_map_ext_1.dat',
     $          'nodes_ext_1.dat',
     $          'sizes_ext_1.dat')
           
           x_map = field_ext_tested%get_x_map()
           y_map = field_ext_tested%get_y_map()
           nodes = field_ext_tested%get_nodes()

           call write_grdpts_id(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_ext_int_1.dat',
     $          'y_map_ext_int_1.dat',
     $          'grdptsid_ext_int_1.dat',
     $          'nodes_ext_int_1.dat',
     $          'sizes_ext_int_1.dat')
           
           call print_interface(
     $          field_ext_tested%domain_extension,
     $          1,
     $          timedev=.true.)


           !2) first integration step computation
           call field_ext_tested%compute_integration_step_ext(
     $          dt, nodes_tmp, time_dev,
     $          compute_1st_step, compute_1st_step_nopt)

           call write_data_ext(
     $          field_ext_tested,
     $          'x_map_ext_2.dat',
     $          'y_map_ext_2.dat',
     $          'nodes_ext_2.dat',
     $          'sizes_ext_2.dat')

           call field_ext_tested%domain_extension%deallocate_after_timeInt()


           !3) nodes synchronization
           call field_ext_tested%apply_bc_on_nodes()

           call write_data_ext(
     $          field_ext_tested,
     $          'x_map_ext_3.dat',
     $          'y_map_ext_3.dat',
     $          'nodes_ext_3.dat',
     $          'sizes_ext_3.dat')

           x_map = field_ext_tested%get_x_map()
           y_map = field_ext_tested%get_y_map()
           nodes = field_ext_tested%get_nodes()

           call write_grdpts_id(
     $          x_map,
     $          y_map,
     $          nodes,
     $          'x_map_ext_int_3.dat',
     $          'y_map_ext_int_3.dat',
     $          'grdptsid_ext_int_3.dat',
     $          'nodes_ext_int_3.dat',
     $          'sizes_ext_int_3.dat')
           
           call print_interface(
     $          field_ext_tested%domain_extension,
     $          3)


           !4) integration cycle
           call field_ext_tested%integrate(dt)

           call write_data_ext(
     $          field_ext_tested,
     $          'x_map_ext_4.dat',
     $          'y_map_ext_4.dat',
     $          'nodes_ext_4.dat',
     $          'sizes_ext_4.dat')
           

        end if

        contains

        subroutine initialize_maps_one_piece(x_map, y_map)

          implicit none

          real(rkind), dimension(:), intent(inout) :: x_map
          real(rkind), dimension(:), intent(inout) :: y_map

          real(rkind) :: dx
          real(rkind) :: dy

          integer(ikind) :: i
          integer(ikind) :: j

          dx = (x_max-x_min)/nx
          dy = (y_max-y_min)/ny

          do i=1, size(x_map,1)
             x_map(i) = x_min + (i-1)*dx
          end do

          do j=1, size(y_map,1)
             y_map(j) = y_min + (j-1)*dy
          end do

        end subroutine initialize_maps_one_piece


        subroutine initialize_four_pieces(
     $     x_map,
     $     y_map,
     $     nodes,
     $     field_ext_tested)

          implicit none

          real(rkind), dimension(:)    , intent(inout) :: x_map
          real(rkind), dimension(:)    , intent(inout) :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          type(field_extended)         , intent(inout) :: field_ext_tested

          real(rkind)                    :: dx
          real(rkind)                    :: dy
          integer(ikind)                 :: i
          integer(ikind)                 :: j
          integer(ikind), dimension(2,2) :: alignment

          dx = (x_max-x_min)/40
          dy = (y_max-y_min)/40

          do i=1, size(x_map,1)
             x_map(i) = x_min + (i+8-1)*dx
          end do

          do j=1, size(y_map,1)
             y_map(j) = y_min + (j+8-1)*dy
          end do

          call initialize_nodes(
     $         x_map,
     $         y_map,
     $         nodes)
          
          call field_ext_tested%set_x_map(x_map)
          call field_ext_tested%set_y_map(y_map)
          call field_ext_tested%set_nodes(nodes)
          
          !extend the interior domain to have
          !four buffer layers
          !----------------------------------
          !add the north buffer layer
          alignment(1,1) = align_W-7
          alignment(1,2) = align_E+7
          alignment(2,1) = align_N
          alignment(2,2) = align_N+7
          call add_sublayer(
     $         field_ext_tested,
     $         N,
     $         alignment,
     $         x_map, y_map, nodes)

          !add the south buffer layer
          alignment(1,1) = align_W-7
          alignment(1,2) = align_E+7
          alignment(2,1) = align_S-7
          alignment(2,2) = align_S
          call add_sublayer(
     $         field_ext_tested,
     $         S,
     $         alignment,
     $         x_map, y_map, nodes)

          !add the east buffer layer
          alignment(1,1) = align_E
          alignment(1,2) = align_E+7
          alignment(2,1) = align_S+1
          alignment(2,2) = align_N-1
          call add_sublayer(
     $         field_ext_tested,
     $         E,
     $         alignment,
     $         x_map, y_map, nodes)

          !add the west buffer layer
          alignment(1,1) = align_W-7
          alignment(1,2) = align_W
          alignment(2,1) = align_S+1
          alignment(2,2) = align_N-1
          call add_sublayer(
     $         field_ext_tested,
     $         W,
     $         alignment,
     $         x_map, y_map, nodes)

        end subroutine initialize_four_pieces


        subroutine initialize_nodes(
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes)

          implicit none

          real(rkind), dimension(:)    , intent(in)  :: interior_x_map
          real(rkind), dimension(:)    , intent(in)  :: interior_y_map
          real(rkind), dimension(:,:,:), intent(out) :: interior_nodes

          type(ic) :: initial_conditions

          call initial_conditions%apply_ic(
     $         interior_nodes,
     $         interior_x_map,
     $         interior_y_map)

c$$$          integer(ikind) :: i,j
c$$$          
c$$$
c$$$          do j=1, size(interior_y_map,1)
c$$$             do i=1, size(interior_x_map,1)
c$$$                interior_nodes(i,j,:) = ini_cond(
c$$$     $               interior_x_map(i),
c$$$     $               interior_y_map(j))
c$$$             end do
c$$$          end do

        end subroutine initialize_nodes

      
        function ini_cond(x,y) result(var)
        
          implicit none
          
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var

          real(rkind) :: amp
          real(rkind) :: r
          real(rkind) :: T
          real(rkind) :: u0
          real(rkind) :: v0

          amp = 1.0
          r   = SQRT(x**2+y**2)
          T   = 0.5d0*(x_max-x_min)
          u0  = 1.0d0
          v0  = 2.0d0

          var(1) = amp*Cos(r*2*ACOS(-1.0d0)/T)
          var(2) = var(1)*u0
          var(3) = var(1)*v0
          var(4) = 0.5d0*var(1)*(u0**2+v0**2)

        end function ini_cond


        subroutine add_sublayer(
     $     field_ext_tested,
     $     mainlayer_id,
     $     alignment,
     $     x_map, y_map, nodes)

          implicit none

          type(field_extended)            , intent(inout) :: field_ext_tested
          integer                         , intent(in)    :: mainlayer_id
          integer(ikind), dimension(4)    , intent(inout) :: alignment
          real(rkind)   , dimension(:)    , intent(in)    :: x_map
          real(rkind)   , dimension(:)    , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes

          type(bf_sublayer), pointer :: added_sublayer

          real(rkind)   , dimension(:)    , allocatable :: bf_x_map
          real(rkind)   , dimension(:)    , allocatable :: bf_y_map
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          integer(ikind), dimension(:,:)  , allocatable :: bf_grdpts_id

          integer(ikind) :: i,j

          added_sublayer => field_ext_tested%domain_extension%allocate_sublayer(
     $         mainlayer_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         alignment)

          call added_sublayer%get_x_map(bf_x_map)
          call added_sublayer%get_y_map(bf_y_map)

          allocate(bf_nodes(size(bf_x_map,1),size(bf_y_map,1),ne))
          call initialize_nodes(
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes)

          call added_sublayer%set_nodes(bf_nodes)

          allocate(bf_grdpts_id(
     $         size(bf_x_map,1),
     $         size(bf_y_map,1)))

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          
          select case(mainlayer_id)
            case(N)
               
               !----------------------------------
               do j=1, size(bf_grdpts_id,2)-bc_size

                  i=1
                  bf_grdpts_id(i,j) = bc_pt

                  i=bc_size
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=bc_size+1, size(bf_grdpts_id,1)-bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do

                  i=size(bf_grdpts_id,1)-bc_size+1
                  bf_grdpts_id(i,j) = bc_interior_pt

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt

               end do

               j=size(bf_grdpts_id,2)-1
               !----------------------------------
               i=1
               bf_grdpts_id(i,j)=bc_pt

               do i=bc_size,size(bf_grdpts_id,1)-1
                  bf_grdpts_id(i,j) = bc_interior_pt
               end do

               i=size(bf_grdpts_id,1)
               bf_grdpts_id(i,j) = bc_pt

               j=size(bf_grdpts_id,2)
               !----------------------------------
               do i=1, size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt
               end do
               
            case(S)

               j=1
               !----------------------------------
               do i=1, size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt
               end do

               j=bc_size
               !----------------------------------
               i=1
               bf_grdpts_id(i,j)=bc_pt

               do i=bc_size,size(bf_grdpts_id,1)-1
                  bf_grdpts_id(i,j) = bc_interior_pt
               end do

               i=size(bf_grdpts_id,1)
               bf_grdpts_id(i,j) = bc_pt

               !----------------------------------
               do j=bc_size+1, size(bf_grdpts_id,2)

                  i=1
                  bf_grdpts_id(i,j) = bc_pt

                  i=bc_size
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=bc_size+1, size(bf_grdpts_id,1)-bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do

                  i=size(bf_grdpts_id,1)-bc_size+1
                  bf_grdpts_id(i,j) = bc_interior_pt

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt

               end do

            case(E)
               
               do j=1, size(bf_grdpts_id,2)
                  do i=1, size(bf_grdpts_id,1) -bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do
                  
                  i=size(bf_grdpts_id,1)-1
                  bf_grdpts_id(i,j) = bc_interior_pt

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt

               end do

            case(W)

               do j=1, size(bf_grdpts_id,2)

                  i=1
                  bf_grdpts_id(i,j) = bc_pt
                  
                  i=bc_size
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=bc_size+1, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = interior_pt
                  end do

               end do

          end select

          call added_sublayer%set_grdpts_id(bf_grdpts_id)
          
          call field_ext_tested%domain_extension%update_integration_borders(added_sublayer)

        end subroutine add_sublayer


        subroutine write_data(
     $     x_map,
     $     y_map,
     $     nodes,
     $     filename_x_map,
     $     filename_y_map,
     $     filename_nodes,
     $     filename_sizes)

          implicit none

          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          character(*)                 , intent(in) :: filename_x_map
          character(*)                 , intent(in) :: filename_y_map
          character(*)                 , intent(in) :: filename_nodes
          character(*)                 , intent(in) :: filename_sizes

          integer :: ios
          
          !x_map
          open(unit=1,
     $          file=filename_x_map,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) x_map
              close(unit=1)
           else
              stop 'file opening pb'
           end if

          !y_map
          open(unit=1,
     $          file=filename_y_map,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=1, iostat=ios) y_map
             close(unit=1)
          else
             stop 'file opening pb'
          end if
          
          !nodes
          open(unit=1,
     $          file=filename_nodes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nodes
              close(unit=1)
           else
              stop 'file opening pb'
           end if

           !sizes
           open(unit=1,
     $          file=filename_sizes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) size(nodes,1),
     $             size(nodes,2),
     $             size(nodes,3)
              close(unit=1)
           else
              stop 'file opening pb'
           end if
          
        end subroutine write_data


        subroutine write_data_ext(
     $     field_ext_tested,
     $     filename_x_map,
     $     filename_y_map,
     $     filename_nodes,
     $     filename_sizes)
        
          implicit none

          type(field_extended), intent(in) :: field_ext_tested
          character(*)        , intent(in) :: filename_x_map
          character(*)        , intent(in) :: filename_y_map
          character(*)        , intent(in) :: filename_nodes
          character(*)        , intent(in) :: filename_sizes

          real(rkind), dimension(40)       :: full_x_map
          real(rkind), dimension(40)       :: full_y_map
          real(rkind), dimension(40,40,ne) :: full_nodes

          real(rkind), dimension(:)    , allocatable :: bf_x_map
          real(rkind), dimension(:)    , allocatable :: bf_y_map
          real(rkind), dimension(:,:,:), allocatable :: bf_nodes

          real(rkind), dimension(nx)       :: int_x_map
          real(rkind), dimension(ny)       :: int_y_map
          real(rkind), dimension(nx,ny,ne) :: int_nodes
          
          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr


          !fill the x_map, y_map and nodes by combining
          !data from the interior and the buffer layers
          !--------------------------------------------
          !combine data from S
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(S)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_nodes_array(bf_nodes)

          full_y_map(1:10)        = bf_y_map(1:10)
          full_nodes(1:40,1:10,:) = bf_nodes(1:40,1:10,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_nodes)          


          !combine data from W
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(W)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_nodes_array(bf_nodes)

          full_x_map(1:10)         = bf_x_map(1:10)
          full_nodes(1:10,11:30,:) = bf_nodes(1:10,3:22,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_nodes)          


          !combine data from interior
          int_x_map = field_ext_tested%get_x_map()
          int_y_map = field_ext_tested%get_y_map()
          int_nodes = field_ext_tested%get_nodes()

          full_x_map(11:30)         = int_x_map(3:22)
          full_y_map(11:30)         = int_y_map(3:22)
          full_nodes(11:30,11:30,:) = int_nodes(3:22,3:22,:)

          !combine data from E
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(E)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_nodes_array(bf_nodes)

          full_x_map(31:40)         = bf_x_map(3:12)
          full_nodes(31:40,11:30,:) = bf_nodes(3:12,3:22,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_nodes)          


          !combine data from N
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(N)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_nodes_array(bf_nodes)

          full_y_map(31:40)        = bf_y_map(3:12)
          full_nodes(1:40,31:40,:) = bf_nodes(1:40,3:12,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_nodes)


          !write the combined data on external binary files
          call write_data(
     $         full_x_map,
     $         full_y_map,
     $         full_nodes,
     $         filename_x_map,
     $         filename_y_map,
     $         filename_nodes,
     $         filename_sizes)          

        end subroutine write_data_ext


        subroutine write_timedev_ext(
     $     field_ext_tested,
     $     int_timedev,
     $     filename_x_map,
     $     filename_y_map,
     $     filename_timedev,
     $     filename_sizes)
        
          implicit none

          type(field_extended)            , intent(in) :: field_ext_tested
          real(rkind), dimension(nx,ny,ne), intent(in) :: int_timedev
          character(*)                    , intent(in) :: filename_x_map
          character(*)                    , intent(in) :: filename_y_map
          character(*)                    , intent(in) :: filename_timedev
          character(*)                    , intent(in) :: filename_sizes

          real(rkind), dimension(40)       :: full_x_map
          real(rkind), dimension(40)       :: full_y_map
          real(rkind), dimension(40,40,ne) :: full_timedev

          real(rkind), dimension(:)    , allocatable :: bf_x_map
          real(rkind), dimension(:)    , allocatable :: bf_y_map
          real(rkind), dimension(:,:,:), allocatable :: bf_timedev

          real(rkind), dimension(nx)       :: int_x_map
          real(rkind), dimension(ny)       :: int_y_map
          
          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr


          !fill the x_map, y_map and timedev by combining
          !data from the interior and the buffer layers
          !--------------------------------------------
          !combine data from S
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(S)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_time_dev(bf_timedev)

          full_y_map(1:10)          = bf_y_map(1:10)
          full_timedev(1:40,1:10,:) = bf_timedev(1:40,1:10,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_timedev)          


          !combine data from W
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(W)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_time_dev(bf_timedev)

          full_x_map(1:10)           = bf_x_map(1:10)
          full_timedev(1:10,11:30,:) = bf_timedev(1:10,3:22,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_timedev)          


          !combine data from interior
          int_x_map = field_ext_tested%get_x_map()
          int_y_map = field_ext_tested%get_y_map()

          full_x_map(11:30)           = int_x_map(3:22)
          full_y_map(11:30)           = int_y_map(3:22)
          full_timedev(11:30,11:30,:) = int_timedev(3:22,3:22,:)

          !combine data from E
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(E)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_time_dev(bf_timedev)

          full_x_map(31:40)           = bf_x_map(3:12)
          full_timedev(31:40,11:30,:) = bf_timedev(3:12,3:22,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_timedev)          


          !combine data from N
          mainlayer_ptr => field_ext_tested%domain_extension%get_mainlayer(N)
          sublayer_ptr  => mainlayer_ptr%get_head_sublayer()
          call sublayer_ptr%get_x_map(bf_x_map)
          call sublayer_ptr%get_y_map(bf_y_map)
          call sublayer_ptr%get_time_dev(bf_timedev)

          full_y_map(31:40)          = bf_y_map(3:12)
          full_timedev(1:40,31:40,:) = bf_timedev(1:40,3:12,:)

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          deallocate(bf_timedev)


          !write the combined data on external binary files
          call write_data(
     $         full_x_map,
     $         full_y_map,
     $         full_timedev,
     $         filename_x_map,
     $         filename_y_map,
     $         filename_timedev,
     $         filename_sizes)

        end subroutine write_timedev_ext


        subroutine write_grdpts_id(
     $     x_map,
     $     y_map,
     $     nodes,
     $     filename_x_map,
     $     filename_y_map,
     $     filename_grdpts_id,
     $     filename_nodes,
     $     filename_sizes)

          implicit none

          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          character(*)                    , intent(in) :: filename_x_map
          character(*)                    , intent(in) :: filename_y_map
          character(*)                    , intent(in) :: filename_grdpts_id
          character(*)                    , intent(in) :: filename_nodes
          character(*)                    , intent(in) :: filename_sizes


          integer, dimension(nx,ny) :: grdpts_id
          integer :: i,j
          integer :: ios


          !grdpts_id
          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)
                grdpts_id(i,j) = interior_pt
             end do
          end do

          !gridpts_id
          open(unit=1,
     $         file=filename_grdpts_id,
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=1, iostat=ios) grdpts_id
             close(unit=1)
          else
             stop 'file opening pb'
          end if

          call write_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         filename_x_map,
     $         filename_y_map,
     $         filename_nodes,
     $         filename_sizes)           

        end subroutine write_grdpts_id


        subroutine print_interface(
     $     interface_used,
     $     index,
     $     timedev)

          implicit none

          class(bf_interface), intent(in) :: interface_used
          integer            , intent(in) :: index
          logical, optional  , intent(in) :: timedev

          logical :: timedev_op

          integer :: format_index

          character(len=15) :: format_nodes
          character(len=18) :: format_grdpt
          character(len=10) :: format_nbsbl
          
          character(len=11) :: filename_x_map
          character(len=11) :: filename_y_map
          character(len=11) :: filename_nodes
          character(len=14) :: filename_grdpt
          character(len=11) :: filename_sizes
          character(len=6 ) :: filename_nbsbl

          if(present(timedev)) then
             timedev_op = timedev
          else
             timedev_op = .false.
          end if

          !determine the number of integer needed to write the
          !file index
          if(index.le.9) then
             format_index = 1
          else
             if((index.ge.10).and.(index.le.99)) then
                format_index = 2
             else
                print '(''test_bf_interface_prog'')'
                print '(''print_output'')'
                stop 'file_index not supported'
             end if
          end if
          

          !determine the format for the name of the output files
          write(format_nodes, '(''(A5,I'',I1,'',A4)'')') format_index
          write(format_grdpt, '(''(A8,I'',I1,'',A4)'')') format_index
          write(format_nbsbl, '(''(I'',I1,'',A4)'')'  ) format_index
          

          !determine the name of the output files
          write(filename_x_map, format_nodes) 'x_map', index, '.dat'
          write(filename_y_map, format_nodes) 'y_map', index, '.dat'
          write(filename_nodes, format_nodes) 'nodes', index, '.dat'
          write(filename_grdpt, format_grdpt) 'grdpt_id', index, '.dat'
          write(filename_sizes, format_nodes) 'sizes', index, '.dat'
          write(filename_nbsbl, format_nbsbl) index, '.dat'


          !write the outputs
          call interface_used%print_binary(
     $         filename_x_map,
     $         filename_y_map,
     $         filename_nodes,
     $         filename_grdpt,
     $         filename_sizes,
     $         filename_nbsbl,
     $         timedev=timedev_op)

        end subroutine print_interface

      end program test_field_extended_time_dev
