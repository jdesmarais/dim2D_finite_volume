      !test the function td_operators%compute_time_dev_nopt()
      program test_bf_compute

        use bc_operators_class, only :
     $       bc_operators

        use bf_compute_class, only :
     $       bf_compute

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_bf_layer, only :
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt
        
        use parameters_input, only :
     $       x_min, x_max,
     $       y_min, y_max,
     $       nx,ny,ne,
     $       bc_size
        
        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt
        
        use sd_operators_class, only :
     $       sd_operators
        
        use td_operators_class, only :
     $       td_operators
        
        implicit none
        
        type(sd_operators) :: s_op
        type(pmodel_eq)    :: p_model
        type(bc_operators) :: bc_op
        type(td_operators) :: td_op

        real(rkind)                                :: dt
        real(rkind)                                :: t
        real(rkind), dimension(:)    , allocatable :: x_map
        real(rkind), dimension(:)    , allocatable :: y_map
        real(rkind), dimension(:,:,:), allocatable :: nodes
        integer    , dimension(:,:)  , allocatable :: grdpts_id

        integer(ikind), dimension(2)                :: x_borders
        integer(ikind), dimension(2)                :: y_borders
        integer       , dimension(:,:), allocatable :: bc_sections
        
        logical :: detailled
        logical :: test_validated

        print '(''**************************'')'
        print '(''physical_model: ns2d'')'
        print '(''bc_operators  : poinsot_xy'')'
        print '(''**************************'')'
        print '()'

        t=0
        dt=1.0

        !allocation
        allocate(x_map(nx))
        allocate(y_map(ny))
        allocate(nodes(nx,ny,ne))
        allocate(grdpts_id(nx,ny))

        !initialize the coordinate tables
        call initialize_maps(
     $       x_map,
     $       y_map)

        !initialize the nodes
        call p_model%apply_ic(
     $       nodes,
     $       x_map,
     $       y_map)

        !initialize the grdpts_id
        call initialize_grdpts_id(
     $       grdpts_id,
     $       bc_sections,
     $       x_borders,
     $       y_borders)
        
        !test compute_time_dev()
        test_validated = test_compute_time_dev(
     $       t, x_map, y_map, nodes, grdpts_id,
     $       s_op, p_model, bc_op, td_op,
     $       x_borders, y_borders,
     $       bc_sections,
     $       detailled)

        print '(''test compute_timedev()'')'
        print '(''test_validated: '',L1)', test_validated
        print '()'


        !test compute_time_integration_step()
        test_validated = test_compute_integration_step(
     $       dt, t, x_map, y_map, nodes, grdpts_id,
     $       s_op, p_model, bc_op, td_op,
     $       x_borders, y_borders,
     $       bc_sections,
     $       detailled)

        print '(''test compute_time_integration_step()'')'
        print '(''test_validated: '',L1)', test_validated
        print '()'


        contains

        subroutine initialize_maps(x_map,y_map)

          implicit none

          real(rkind), dimension(:), intent(inout) :: x_map
          real(rkind), dimension(:), intent(inout) :: y_map

          real(rkind)    :: dx,dy
          integer(ikind) :: i,j

          dx = 0.5
          dy = 0.6
          
          do i=1, nx
             x_map(i) = x_min + (i-1)*dx
          end do
          
          do j=1, ny
             y_map(j) = y_min + (j-1)*dy
          end do
          
        end subroutine initialize_maps


        subroutine initialize_grdpts_id(
     $     grdpts_id,
     $     bc_sections,
     $     x_borders,
     $     y_borders)

          implicit none

          integer       , dimension(:,:)             , intent(inout) :: grdpts_id
          integer       , dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind), dimension(2)               , intent(out)   :: x_borders
          integer(ikind), dimension(2)               , intent(out)   :: y_borders
          

          !integration borders
          x_borders = [1,nx]
          y_borders = [1,ny]

          !bc_sections
          allocate(bc_sections(4,8))

          bc_sections(1,1) = SW_corner_type
          bc_sections(2,1) = 1
          bc_sections(3,1) = 1
          
          bc_sections(1,2) = S_edge_type
          bc_sections(2,2) = bc_size+1
          bc_sections(3,2) = 1
          bc_sections(4,2) = nx-bc_size

          bc_sections(1,3) = SE_corner_type
          bc_sections(2,3) = nx-bc_size+1
          bc_sections(3,3) = 1

          bc_sections(1,4) = W_edge_type
          bc_sections(2,4) = 1
          bc_sections(3,4) = bc_size+1
          bc_sections(4,4) = ny-bc_size

          bc_sections(1,5) = E_edge_type
          bc_sections(2,5) = nx-bc_size+1
          bc_sections(3,5) = bc_size+1
          bc_sections(4,5) = ny-bc_size

          bc_sections(1,6) = NW_corner_type
          bc_sections(2,6) = 1
          bc_sections(3,6) = ny-bc_size+1
          
          bc_sections(1,7) = N_edge_type
          bc_sections(2,7) = bc_size+1
          bc_sections(3,7) = ny-bc_size+1
          bc_sections(4,7) = nx-bc_size

          bc_sections(1,8) = NE_corner_type
          bc_sections(2,8) = nx-bc_size+1
          bc_sections(3,8) = ny-bc_size+1

          !grdpts_id
          call initialize_grdpts(grdpts_id)
          
        end subroutine initialize_grdpts_id


        function test_compute_time_dev(
     $     t, x_map, y_map, nodes, grdpts_id,
     $     s_op, p_model, bc_op, td_op,
     $     x_borders, y_borders,
     $     bc_sections,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)                                , intent(in)    :: t
          real(rkind), dimension(:)                  , intent(in)    :: x_map
          real(rkind), dimension(:)                  , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)              , intent(in)    :: nodes
          integer    , dimension(:,:)                , intent(in)    :: grdpts_id
          type(sd_operators)                         , intent(in)    :: s_op
          type(pmodel_eq)                            , intent(in)    :: p_model
          type(bc_operators)                         , intent(in)    :: bc_op
          type(td_operators)                         , intent(in)    :: td_op
          integer(ikind), dimension(2)               , intent(in)    :: x_borders
          integer(ikind), dimension(2)               , intent(in)    :: y_borders
          integer       , dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(in)    :: detailled
          logical                                                    :: test_validated

          real(rkind)   , dimension(:,:,:), allocatable :: timedev
          real(rkind)   , dimension(:,:,:), allocatable :: timedev_bf
          type(bf_compute)                              :: bf_compute_used

          allocate(timedev(nx,ny,ne))

          !compute the timedev for the interior
          !using the optimized function
          timedev = td_op%compute_time_dev(
     $         t, nodes, x_map, y_map, s_op,
     $         p_model, bc_op)

          !compute the timedev for the buffer layer
          !using the non-optimized function
          call bf_compute_used%allocate_tables(
     $         size(timedev,1),
     $         size(timedev,2))

          call bf_compute_used%set_bc_sections(bc_sections)

          call bf_compute_used%compute_time_dev(
     $         td_op,
     $         t, nodes, x_map, y_map,
     $         s_op, p_model, bc_op,
     $         grdpts_id,
     $         x_borders, y_borders)

          call bf_compute_used%get_time_dev(timedev_bf)

          call bf_compute_used%deallocate_tables()        
        
          !compare the two ways time dev are computed
          test_validated = compare_timedev(
     $         timedev, timedev_bf,
     $         detailled)
          
          deallocate(timedev)

        end function test_compute_time_dev


        function test_compute_integration_step(
     $     dt, t, x_map, y_map, nodes, grdpts_id,
     $     s_op, p_model, bc_op, td_op,
     $     x_borders, y_borders,
     $     bc_sections,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)                                , intent(in)    :: dt
          real(rkind)                                , intent(in)    :: t
          real(rkind), dimension(:)                  , intent(in)    :: x_map
          real(rkind), dimension(:)                  , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)              , intent(inout) :: nodes
          integer    , dimension(:,:)                , intent(in)    :: grdpts_id
          type(sd_operators)                         , intent(in)    :: s_op
          type(pmodel_eq)                            , intent(in)    :: p_model
          type(bc_operators)                         , intent(in)    :: bc_op
          type(td_operators)                         , intent(in)    :: td_op
          integer(ikind), dimension(2)               , intent(in)    :: x_borders
          integer(ikind), dimension(2)               , intent(in)    :: y_borders
          integer       , dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(in)    :: detailled
          logical                                                    :: test_validated

          real(rkind)   , dimension(:,:,:), allocatable :: timedev
          real(rkind)   , dimension(:,:,:), allocatable :: nodes_tmp
          real(rkind)   , dimension(:,:,:), allocatable :: nodes_bf
          type(bf_compute)                              :: bf_compute_used

          integer(ikind) :: i,j
          integer        :: k

          allocate(timedev(nx,ny,ne))
          allocate(nodes_tmp(nx,ny,ne))
          allocate(nodes_bf(nx,ny,ne))

          !make a copy of nodes
          do k=1, size(nodes,3)
             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
                   nodes_bf(i,j,k) = nodes(i,j,k)
                end do
             end do
          end do

          !compute the timedev for the interior
          !using the optimized function
          timedev = td_op%compute_time_dev(
     $         t, nodes, x_map, y_map, s_op,
     $         p_model, bc_op)

          !compute the integration step
          call compute_1st_step(
     $         nodes, dt, nodes_tmp, timedev,
     $         x_borders, y_borders)

          !compute the timedev for the buffer layer
          !using the non-optimized function
          call bf_compute_used%allocate_tables(
     $         size(timedev,1),
     $         size(timedev,2))

          call bf_compute_used%set_bc_sections(bc_sections)

          call bf_compute_used%compute_time_dev(
     $         td_op,
     $         t, nodes_bf, x_map, y_map,
     $         s_op, p_model, bc_op,
     $         grdpts_id,
     $         x_borders, y_borders)

          call bf_compute_used%compute_integration_step(
     $         grdpts_id, nodes_bf, dt,
     $         x_borders, y_borders,
     $         compute_1st_step_nopt)

          call bf_compute_used%deallocate_tables()
        
          !compare the two ways time dev are computed
          test_validated = compare_timedev(
     $         nodes, nodes_bf,
     $         detailled)

          deallocate(timedev)
          deallocate(nodes_tmp)
          deallocate(nodes_bf)

        end function test_compute_integration_step


        subroutine initialize_grdpts(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          j=1
          do i=1,size(grdpts_id,1)
             grdpts_id(i,j) = bc_pt
          end do

          j=2
          grdpts_id(1,j) = bc_pt
          do i=bc_size,size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j) = bc_interior_pt
          end do
          grdpts_id(nx,j)=bc_pt

          do j=bc_size+1, size(grdpts_id,2)-bc_size

             grdpts_id(1,j) = bc_pt
             grdpts_id(2,j) = bc_interior_pt

             do i=bc_size+1, size(grdpts_id,1)-bc_size
                grdpts_id(i,j) = interior_pt
             end do

             grdpts_id(nx-bc_size+1,j) = bc_interior_pt
             grdpts_id(nx,j)           = bc_pt

          end do

          j=size(grdpts_id,2)-bc_size+1
          grdpts_id(1,j) = bc_pt
          do i=bc_size,size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j) = bc_interior_pt
          end do
          grdpts_id(nx,j)=bc_pt
          
          j=size(grdpts_id,2)
          do i=1,size(grdpts_id,1)
             grdpts_id(i,j) = bc_pt
          end do

        end subroutine initialize_grdpts


        function compare_timedev(
     $       timedev, timedev_bf,
     $       detailled)
     $       result(test_validated)

          real(rkind), dimension(:,:,:), intent(in) :: timedev
          real(rkind), dimension(:,:,:), intent(in) :: timedev_bf
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          integer(ikind) :: i,j
          integer        :: k
          logical        :: same

          test_validated = .true.

          do k=1, size(timedev,3)
             do j=1, size(timedev,2)
                do i=1, size(timedev,1)
                   
                   same = is_test_validated(
     $                  timedev(i,j,k),
     $                  timedev_bf(i,j,k),
     $                  .false.)
                   
                   if((.not.same).and.detailled) then
                      
                      print '(I3,1X,I3,1X,I3,'' timedev: '',2I10)',
     $                     i,j,k, 
     $                     int(timedev(i,j,k)*1e5),
     $                     int(timedev_bf(i,j,k)*1e5)
                      
                   end if
                   
                   test_validated = test_validated.and.same
                   
                end do
             end do
          end do

        end function compare_timedev


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1

        end function is_test_validated


      end program test_bf_compute
