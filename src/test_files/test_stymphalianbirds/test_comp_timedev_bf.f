      program test_comp_timedev_bf

       use bc_operators_class, only :
     $     bc_operators

        use parameters_bf_layer, only :
     $     interior_pt

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
        dx = 0.1
        dy = 0.2

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
        do j=1, size(grdpts_id,2)
           do i=1, size(grdpts_id,1)
              grdpts_id(i,j) = interior_pt
           end do
        end do

        !compute the timedev for the interior
        !using the non-optimized subroutine
        x_borders = [bc_size+1,nx-bc_size]
        y_borders = [bc_size+1,ny-bc_size]

        call td_op%compute_time_dev_nopt(
     $       t,nodes,x_map,y_map,
     $       s_op,p_model,bc_op,
     $       timedev_bf,
     $       grdpts_id,
     $       bc_sections,
     $       x_borders, y_borders)

        
        !compare the two ways fluxes are computed
        test_validated = compare_timedev(
     $       timedev, timedev_bf,
     $       detailled)

        print '(''test fluxes comparison'')'
        print '(''test_validated: '',L1)', test_validated
        print '()'


        contains


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
             do j=bc_size+1, size(timedev,2)-bc_size
                do i=bc_size+1, size(timedev,1)-bc_size
                   
                   same = is_test_validated(
     $                  timedev(i,j,k),
     $                  timedev_bf(i,j,k),
     $                  .false.)
                   
                   if((.not.same).and.detailled) then
                      
                      print '(I3,1X,I3,1X,I3,'' flux_x: '',2I8)',
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
