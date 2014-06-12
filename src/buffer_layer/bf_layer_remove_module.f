      module bf_layer_remove_module

        use bf_activation_module  , only : is_openbc_undermined
        use bf_layer_errors_module, only : error_mainlayer_id
        use parameters_bf_layer   , only : align_N,align_S,
     $                                     align_E,align_W
        use parameters_constant   , only : N,S,E,W
        use parameters_input      , only : nx,ny,ne,search_dcr
        use parameters_kind       , only : ikind, rkind
        
        implicit none

        private
        public :: check_if_can_be_removed_locally,
     $            get_check_line_param


        contains


        !< check if the grid points of the buffer layer common
        !> with the interior domain are such that the buffer layer
        !> can be removed without considering its neighbor dependencies
        !
        !> @param bf_localization
        !> cardinal coordinate identifying the position of the buffer
        !> layer
        !> 
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !>
        !> @param bf_nodes
        !> grid points of the buffer layer investigated
        !>
        !> @param interior_nodes
        !> grid points of the interior domain
        !>
        !> @param bf_removal
        !> logical identifying whether the local removal is approved or
        !> not
        !---------------------------------------------------------------
        subroutine check_if_bf_layer_remains(
     $       bf_localization, bf_alignment, bf_nodes,
     $       interior_nodes, bf_remains)

          implicit none

          integer                          , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2)   , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)     , intent(in)  :: bf_match_table
          real(rkind), dimension(:,:,:)    , intent(in)  :: bf_nodes
          real(rkind), dimension(nx,ny,ne) , intent(in)  :: interior_nodes
          logical                          , intent(out) :: bf_remains

          
          integer(ikind), dimension(2,2) :: bf_coords
          integer(ikind), dimension(2,2) :: in_coords
          

          !check if the buffer layer has grid points in common with
          !the interior domain
          if(grdpts_common_with_interior(bf_alignment)) then
             
             !depending on the buffer layer position, determine
             !the line around which the grid points will be checked
             !to see whether the region undermines the open boundary
             !conditions or not
             call get_check_line_param(
     $            bf_localization, bf_alignment, bf_match_table,
     $            bf_coords, in_coords)

             !check the neighboring points around the line
             call check_line_neighbors(
     $            bf_coords, in_coords,
     $            bf_nodes, interior_nodes,
     $            bf_remains)
             
          
          !if the buffer layer has no grid point in common with the
          !interior domain, it can be removed immediately (the
          !neighboring buffer layers depending on this buffer layer
          !are not taken into account)
          else
             bf_remains = .false.
          end if

        end subroutine check_if_can_be_removed_locally


        !< check if the alignment of the buffer layer compared to
        !> the interior is such that the buffer layer has grid points
        !> in common with the interior domain
        function grdpts_common_with_interior(bf_alignment)
     $     result(grdpts_common)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          logical                                    :: grdpts_common

          integer(ikind) :: common_layer_size_x
          integer(ikind) :: common_layer_size_y


          common_layer_size_x = min(nx-bc_size,bf_alignment(1,2))-
     $                          max(bc_size+1,bf_alignment(1,1))+
     $                          2*bc_size+1

          common_layer_size_y = min(ny-bc_size,bf_alignment(2,2))-
     $                          max(bc_size+1,bf_alignment(2,1))+
     $                          2*bc_size+1          

          grdpts_common = (common_layer_size_x.gt.0).and.
     $                    (common_layer_size_y.gt.0)


        end function grdpts_common_with_interior


        !< determine the parameters constraining the line around
        !> which the neighboring points will be checked to see
        !> whether the open boundary conditions are undermined
        !> in this region or not
        subroutine get_check_line_param(
     $     bf_localization, bf_alignment, bf_match_table,
     $     bf_coords, in_coords)


          implicit none

          integer                       , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf_match_table
          integer(ikind), dimension(2,2), intent(out) :: bf_coords
          integer(ikind), dimension(2,2), intent(out) :: in_coords
          

          select case(bf_localization)

            case(N)
               bf_coords(1,1) = max(bf_alignment(1,1) - bc_size, 1) - bf_match_table(1)
               bf_coords(2,1) = max(1, align_N - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_alignment(1,2) + bc_size, nx) - bf_match_table(1)
               bf_coords(2,2) = min(bf_size_y, align_N + search_dcr - bf_match_table(2))


               in_coords(1,1) = bf_coords(1,1)
               in_coords(2,1) = align_N - search_dcr

               in_coords(1,2) = bf_coords(1,2)
               in_coords(2,2) = bf_coords(2,1)-1 + bf_match_table(2)


            case(S)
               bf_coords(1,1) = max(bf_alignment(1,1) - bc_size, 1) - bf_match_table(1)
               bf_coords(2,1) = max(1, align_S - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_alignment(1,2) + bc_size, nx) - bf_match_table(1)
               bf_coords(2,2) = min(bf_size_y, align_S + search_dcr - bf_match_table(2))


               in_coords(1,1) = bf_coords(1,1)
               in_coords(2,1) = bf_coords(2,2)+1 + bf_match_table(2)

               in_coords(1,2) = bf_coords(1,2)
               in_coords(2,2) = align_S + search_dcr


            case(E)
               bf_coords(1,1) = max(1, align_E - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(bf_alignment(2,1) - bc_size, 1) - bf_match_table(2)

               bf_coords(1,2) = min(bf_size_x, align_E + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_alignment(2,2) + bc_size, ny) - bf_match_table(2)


               in_coords(1,1) = align_E - search_dcr
               in_coords(2,1) = bf_coords(2,1)

               in_coords(1,2) = bf_coords(1,1)-1 + bf_match_table(1)
               in_coords(2,2) = bf_coords(2,2)


            case(W)
               bf_coords(1,1) = max(1, align_W - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(bf_alignment(2,1) - bc_size, 1) - bf_match_table(2)

               bf_coords(1,2) = min(bf_size_x, align_W + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_alignment(2,2) + bc_size, ny) - bf_match_table(2)


               in_coords(1,1) = bf_coords(1,2)+1 + bf_match_table(1)
               in_coords(2,1) = bf_coords(2,1)

               in_coords(1,2) = align_W + search_dcr
               in_coords(2,2) = bf_coords(2,2)


            case default
               call error_mainlayer_id(
     $              'bf_layer_remove_module.f',
     $              'get_check_line_param',
     $              bf_localization)

          end select
          
        end subroutine get_check_line_param


        !< check the points around the line
        subroutine check_line_neighbors(
     $     bf_coords, in_coords,
     $     bf_nodes, interior_nodes,
     $     bf_remains)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_coords
          integer(ikind), dimension(2,2)     , intent(in)    :: in_coords
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          logical                            , intent(out)   :: bf_remains


          bf_remains = .false.


          !check the interior points
          call check_layer(
     $         in_coords,
     $         interior_nodes,
     $         bf_remains)
          
          
          !check the buffer layer points
          if(bf_remains) then
             call check_layer(
     $            bf_coords,
     $            bf_nodes,
     $            bf_remains)
          end if


        end subroutine check_line_neighbors


        !< check if the grid points neighboring the central point
        !> identified by cpt_coords undermine the open boundary
        !> conditions. The number of grid points to be checked
        !> is reduced knowing the previous central point investigated
        !> identified by cpt_coords_p
        subroutine check_layer(
     $     pt_coords,
     $     nodes,
     $     bf_remains)
        
          implicit none
          
          integer(ikind), dimension(2,2)     , intent(in) :: pt_coords
          real(rkind)   , dimension(:,:,:)   , intent(in) :: nodes
          logical                            , intent(out):: bf_remains
          
          
          integer(ikind) :: i,j


          do j=pt_coords(2,1), pt_coords(2,2)
             do i=pt_coords(1,1), pt_coords(1,2)

                bf_remains = is_openbc_undermined(nodes(i,j,:))

                if(bf_remains) then
                   exit
                end if

             end do

             if(bf_remains) then
                exit
             end if

          end do

        end subroutine check_layer

      end module bf_layer_remove_module
