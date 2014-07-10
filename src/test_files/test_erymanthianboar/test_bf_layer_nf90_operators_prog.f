      program test_bf_layer_nf90_operators_prog

        use bf_layer_nf90_operators_module, only :print_bf_layer_on_netcdf
        use dim2d_eq_class                , only : dim2d_eq
        use parameters_bf_layer           , only : no_pt
        use parameters_constant           , only : N
        use parameters_kind               , only : ikind, rkind
        use parameters_input              , only : ne

        implicit none

        integer, parameter :: n_x = 10
        integer, parameter :: n_y = 7   

        integer    , dimension(n_x,n_y)    :: grdpts_id
        real(rkind), dimension(n_x,n_y,ne) :: nodes
        real(rkind) :: time

        real(rkind) :: x_start
        real(rkind) :: y_start
        real(rkind) :: dx
        real(rkind) :: dy

        type(dim2d_eq) :: p_model

        call initialize_data(grdpts_id, nodes, time)

        x_start = 1.0
        y_start = 2.0
        dx = 0.5
        dy = 0.6

        call print_bf_layer_on_netcdf(
     $       'test_bf_layer.nc',
     $       p_model%get_var_name(),
     $       p_model%get_var_longname(),
     $       p_model%get_var_unit(),
     $       N, 1,
     $       x_start, y_start, dx, dy,
     $       grdpts_id, nodes, time)

        
        contains

        subroutine initialize_data(grdpts_id, nodes, time)

          implicit none

          integer    , dimension(n_x,n_y)   , intent(out) :: grdpts_id
          real(rkind), dimension(n_x,n_y,ne), intent(out) :: nodes
          real(rkind)                       , intent(out) :: time

          integer(ikind) :: i,j
          integer        :: k

          !initialization of the grdpts_id
          !the diagonal is inialized with no_pt
          do j=1, n_y
             do i=1, n_x
                grdpts_id(i,j) = mod(i+j,2)+2
             end do
          end do

          do i=1, min(n_x,n_y)
             grdpts_id(i,i) = no_pt
          end do

          !initialization of the nodes
          do k=1, ne
             do j=1, n_y
                do i=1, n_x
                   nodes(i,j,k) = (i-1) + 10*(j-1) + 100*(k-1)
                end do
             end do
          end do

          !initialization of the time
          time = 1.0d0

        end subroutine initialize_data


      end program test_bf_layer_nf90_operators_prog
