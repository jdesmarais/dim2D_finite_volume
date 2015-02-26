      program test_bf_layer_nf90_operators_prog

        use bf_layer_nf90_operators_module, only :
     $     print_bf_layer_on_netcdf

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_constant, only :
     $       N

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use parameters_input, only :
     $       ne

        implicit none

        integer, parameter :: n_x = 10
        integer, parameter :: n_y = 7   

        integer    , dimension(n_x,n_y)    :: grdpts_id
        real(rkind), dimension(n_x)        :: x_map
        real(rkind), dimension(n_y)        :: y_map
        real(rkind), dimension(n_x,n_y,ne) :: nodes
        real(rkind)                        :: time

        type(pmodel_eq) :: p_model

        call initialize_data(grdpts_id, x_map, y_map, nodes, time)

        call print_bf_layer_on_netcdf(
     $       'test_bf_layer.nc',
     $       p_model%get_var_name(),
     $       p_model%get_var_longname(),
     $       p_model%get_var_unit(),
     $       N, 1,
     $       grdpts_id, x_map, y_map, nodes, time)

        print '(''check file test_bf_layer.nc'')'

        
        contains

        subroutine initialize_data(grdpts_id, x_map, y_map, nodes, time)

          implicit none

          integer    , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind), dimension(:)    , intent(inout) :: x_map
          real(rkind), dimension(:)    , intent(inout) :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(out)   :: time

          integer(ikind) :: i,j
          integer        :: k

          real(rkind) :: x_start
          real(rkind) :: y_start
          real(rkind) :: dx
          real(rkind) :: dy


          !initialization of the grdpts_id
          !the diagonal is inialized with no_pt
          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)
                grdpts_id(i,j) = mod(i+j,2)+2
             end do
          end do

          do i=1, min(size(grdpts_id,1),size(grdpts_id,2))
             grdpts_id(i,i) = no_pt
          end do

          !initialization of the x_map
          x_start = 1.0
          dx = 0.5

          do i=1, size(grdpts_id,1)
             x_map(i) = x_start + (i-1)*dx
          end do

          !initialization of the y_map
          y_start = 2.0
          dy = 0.6

          do j=1, size(grdpts_id,2)
             y_map(j) = y_start + (j-1)*dy
          end do

          !initialization of the nodes
          do k=1, ne
             do j=1, size(grdpts_id,2)
                do i=1, size(grdpts_id,1)
                   nodes(i,j,k) = (i-1) + 10*(j-1) + 100*(k-1)
                end do
             end do
          end do

          !initialization of the time
          time = 1.0d0

        end subroutine initialize_data


      end program test_bf_layer_nf90_operators_prog
