      program test_comp_timedev_bf

        use bc_operators_class, only :
     $      bc_operators

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
     $      bc_pt,
     $      bc_interior_pt,
     $      interior_pt

        use parameters_input, only :
     $     x_min, x_max,
     $     y_min, y_max,
     $     nx,ny,ne,
     $     bc_size

        use parameters_kind, only :
     $     ikind,
     $     rkind

        use pmodel_eq_class, only :
     $     pmodel_eq

        use sd_operators_class, only :
     $     sd_operators

        use td_operators_class, only :
     $     td_operators

        implicit none

        type(pmodel_eq)    :: p_model
        type(sd_operators) :: s_op
        type(bc_operators) :: bc_op
        type(td_operators) :: td_op
        real(rkind)        :: dx
        real(rkind)        :: dy

        real(rkind), dimension(:,:,:), allocatable :: nodes
        real(rkind), dimension(:)    , allocatable :: x_map
        real(rkind), dimension(:)    , allocatable :: y_map

        integer    , dimension(:,:)  , allocatable :: grdpts_id
        real(rkind), dimension(:,:,:), allocatable :: timedev
        real(rkind), dimension(:,:,:), allocatable :: timedev_bf
                
        integer(ikind), dimension(2) :: x_borders
        integer(ikind), dimension(2) :: y_borders

        integer, dimension(:,:), allocatable :: bc_sections

        integer(ikind) :: i,j
        real(rkind)    :: t

        logical :: detailled
        logical :: test_validated

        print '(''**************************'')'
        print '(''physical_model: ns2d'')'
        print '(''bc_operators  : poinsot_xy'')'
        print '(''**************************'')'
        print '()'


        !allocation
        allocate(nodes(nx,ny,ne))
        allocate(x_map(nx))
        allocate(y_map(ny))

        allocate(grdpts_id(nx,ny))

        allocate(timedev(nx,ny,ne))
        allocate(timedev_bf(nx,ny,ne))


        !initialize the coordinate tables
        dx = 0.5
        dy = 0.6

        do i=1, nx
           x_map(i) = x_min + (i-1)*dx
        end do

        do j=1, ny
           y_map(j) = y_min + (j-1)*dy
        end do

        !initialize the nodes
        call p_model%apply_ic(nodes,x_map,y_map)


        !compute the timedev for the interior
        !using the optimized function
        timedev = td_op%compute_time_dev(
     $       t, nodes, x_map, y_map, s_op,
     $       p_model, bc_op)


        !initialize the grdpts_id
        call initialize_grdpts(grdpts_id)
        

        !compute the timedev for the interior
        !using the non-optimized subroutine
        x_borders = [bc_size+1,nx-bc_size]
        y_borders = [bc_size+1,ny-bc_size]

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

        call td_op%compute_time_dev_nopt(
     $       t,nodes,x_map,y_map,
     $       s_op,p_model,bc_op,
     $       timedev_bf,
     $       grdpts_id,
     $       bc_sections,
     $       x_borders, y_borders)

        !compare the two ways time dev are computed
        detailled = .true.
        test_validated = compare_timedev(
     $       timedev, timedev_bf,
     $       detailled)

        print '(''test timedev comparison'')'
        print '(''test_validated: '',L1)', test_validated
        print '()'


        contains

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


      end program test_comp_timedev_bf
