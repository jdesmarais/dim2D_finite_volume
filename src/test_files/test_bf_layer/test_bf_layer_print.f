      program test_bf_layer_print

        use bf_layer_print_class, only :
     $       bf_layer_print

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


        call test_print_netcdf()

        
        contains

        subroutine test_print_netcdf()

          implicit none          
          
          type(bf_layer_print) :: bf_layer_used
          real(rkind)          :: time
          type(pmodel_eq)      :: p_model
          
          call initialize_data(bf_layer_used, time)
          
          call bf_layer_used%print_netcdf(
     $         'test_bf_layer.nc',
     $         p_model%get_var_name(),
     $         p_model%get_var_longname(),
     $         p_model%get_var_unit(),
     $         1, time)

          print '(''check file test_bf_layer.nc'')'

        end subroutine test_print_netcdf


        subroutine initialize_data(bf_layer_used,time)

          implicit none

          type(bf_layer_print), intent(inout) :: bf_layer_used
          real(rkind)         , intent(out)   :: time

          integer, parameter :: n_x = 10
          integer, parameter :: n_y = 7

          integer    , dimension(:,:)  , allocatable :: grdpts_id
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes

          integer(ikind) :: i,j
          integer        :: k

          real(rkind) :: x_start
          real(rkind) :: y_start
          real(rkind) :: dx
          real(rkind) :: dy

          allocate(grdpts_id(n_x,n_y))
          allocate(x_map(n_x))
          allocate(y_map(n_y))
          allocate(nodes(n_x,n_y,ne))


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

          !initialization of the buffer layer
          call bf_layer_used%ini(N)
          call bf_layer_used%set_grdpts_id(grdpts_id)
          call bf_layer_used%set_x_map(x_map)
          call bf_layer_used%set_y_map(y_map)
          call bf_layer_used%set_nodes(nodes)

        end subroutine initialize_data

      end program test_bf_layer_print
