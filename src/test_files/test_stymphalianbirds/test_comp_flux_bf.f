      program test_comp_flux_bf

        use parameters_bf_layer, only :
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

        use sd_operators_class, only :
     $       sd_operators

        implicit none

        type(pmodel_eq)    :: p_model
        type(sd_operators) :: s_op
        real(rkind)        :: dx
        real(rkind)        :: dy

        real(rkind), dimension(:,:,:), allocatable :: nodes
        real(rkind), dimension(:)    , allocatable :: x_map
        real(rkind), dimension(:)    , allocatable :: y_map

        real(rkind), dimension(:,:,:), allocatable :: flux_x
        real(rkind), dimension(:,:,:), allocatable :: flux_y

        integer    , dimension(:,:)  , allocatable :: grdpts_id
        real(rkind), dimension(:,:,:), allocatable :: flux_x_bf
        real(rkind), dimension(:,:,:), allocatable :: flux_y_bf
                
        integer(ikind), dimension(2) :: x_borders
        integer(ikind), dimension(2) :: y_borders

        integer(ikind) :: i,j

        logical        :: detailled
        logical        :: test_validated


        !allocation
        allocate(nodes(nx,ny,ne))
        allocate(x_map(nx))
        allocate(y_map(ny))

        allocate(flux_x(nx+1,ny,ne))
        allocate(flux_y(nx,ny+1,ne))

        allocate(grdpts_id(nx,ny))
        allocate(flux_x_bf(nx+1,ny,ne))
        allocate(flux_y_bf(nx,ny+1,ne))


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


        !compute the fluxes for the interior
        !using the optimized function
        flux_x = p_model%compute_flux_x(nodes,dx,dy,s_op)
        flux_y = p_model%compute_flux_y(nodes,dx,dy,s_op)
              

        !initialize the grdpts_id
        do j=1, size(grdpts_id,2)
           do i=1, size(grdpts_id,1)
              grdpts_id(i,j) = interior_pt
           end do
        end do        


        !compute the fluxes for the buffer layers
        !using the non-optimized subroutine
        x_borders = [bc_size+1,nx-bc_size]
        y_borders = [bc_size+1,ny-bc_size]

        call p_model%compute_flux_x_nopt(
     $       nodes,dx,dy,s_op,
     $       grdpts_id,
     $       flux_x_bf,
     $       x_borders, y_borders)

        call p_model%compute_flux_y_nopt(
     $       nodes,dx,dy,s_op,
     $       grdpts_id,
     $       flux_y_bf,
     $       x_borders, y_borders)

        
        !compare the two ways fluxes are computed
        test_validated = compare_fluxes(
     $       flux_x, flux_y,
     $       flux_x_bf, flux_y_bf,
     $       detailled)

        print '(''test fluxes comparison'')'
        print '(''test_validated: '',L1)', test_validated
        print '()'


        contains


        function compare_fluxes(
     $       flux_x, flux_y,
     $       flux_x_bf, flux_y_bf,
     $       detailled)
     $       result(test_validated)

          real(rkind), dimension(:,:,:), intent(in) :: flux_x
          real(rkind), dimension(:,:,:), intent(in) :: flux_y
          real(rkind), dimension(:,:,:), intent(in) :: flux_x_bf
          real(rkind), dimension(:,:,:), intent(in) :: flux_y_bf
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          integer(ikind) :: i,j
          integer        :: k
          logical        :: same

          test_validated = .true.

          do k=1, size(flux_x,3)
             do j=bc_size+1, size(flux_x,2)-bc_size
                do i=bc_size+1, size(flux_x,1)-bc_size
                   
                   same = is_test_validated(
     $                  flux_x(i,j,k),
     $                  flux_x_bf(i,j,k),
     $                  .false.)
                   
                   if((.not.same).and.detailled) then
                      
                      print '(I3,1X,I3,1X,I3,'' flux_x: '',2I8)',
     $                     i,j,k, 
     $                     int(flux_x(i,j,k)*1e5),
     $                     int(flux_x_bf(i,j,k)*1e5)
                      
                   end if
                   
                   test_validated = test_validated.and.same
                   
                end do
             end do
          end do
          
          do k=1, size(flux_y,3)
             do j=bc_size+1, size(flux_y,2)-bc_size
                do i=bc_size+1, size(flux_y,1)-bc_size

                   same = is_test_validated(
     $                  flux_y(i,j,k),
     $                  flux_y_bf(i,j,k),
     $                  .false.)

                   if((.not.same).and.detailled) then

                      print '(I3,1X,I3,1X,I3,'' flux_y: '',2I8)',
     $                     i,j,k, 
     $                     int(flux_y(i,j,k)*1e5),
     $                     int(flux_y_bf(i,j,k)*1e5)

                   end if

                   test_validated = test_validated.and.same

                end do
             end do
          end do

        end function compare_fluxes


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


      end program test_comp_flux_bf
