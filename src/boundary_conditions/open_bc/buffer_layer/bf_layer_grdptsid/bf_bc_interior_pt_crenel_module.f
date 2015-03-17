      !> @file
      !> detection of bc_interior_pt crenels
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> detection of bc_interior_pt crenels
      !
      !> @date
      ! 17_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_bc_interior_pt_crenel_module

        use bf_newgrdpt_dispatch_module, only :
     $     are_grdpts_available_to_get_newgrdpt_proc

        use parameters_bf_layer, only :
     $       bc_pt

        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public ::
     $       check_if_bc_interior_pt_crenel,
     $       are_grdpts_available_to_detect_bc_interior_pt_crenel


        contains


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the bc_interior_pt leads to a bc_interior_pt
        !> crenel that should be removed (i.e. whether the suspicious
        !> bc_interior_pt has a bc_pt neighbor)
        !
        !> @date
        !> 17_03_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> configuration of the grid points around the suspicious
        !> bc_interior_pt
        !
        !> @param i
        !> x-index of the bc_interior_pt checked
        !
        !> @param j
        !> y-index of the bc_interior_pt checked
        !
        !> @return is_bc_interior_crenel
        !> logical determining whether this is a bc_crenel or not
        !--------------------------------------------------------------        
        function check_if_bc_interior_pt_crenel(grdpts_id,i_center,j_center)
     $       result(is_a_bc_interior_pt_crenel)

          implicit none

          integer       , dimension(:,:), intent(in) :: grdpts_id
          integer(ikind)                , intent(in) :: i_center
          integer(ikind)                , intent(in) :: j_center
          logical                                    :: is_a_bc_interior_pt_crenel
          

          integer(ikind) :: i,j
          

          is_a_bc_interior_pt_crenel = .true.


          do j=j_center-1,j_center+1
             do i=i_center-1, i_center+1
                if(grdpts_id(i,j).eq.bc_pt) then
                   is_a_bc_interior_pt_crenel = .false.
                   exit
                end if                
             end do
             if(.not.is_a_bc_interior_pt_crenel) exit
          end do

        end function check_if_bc_interior_pt_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are enough grid points to
        !> identify whether there is a bc_interior_pt crenel
        !> or not
        !
        !> @date
        !> 17_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the mainlayer to which the
        !> buffer layer belongs
        !
        !>@param size_x
        !> size along the x-direction of the buffer layer
        !
        !>@param size_y
        !> size along the y-direction of the buffer layer
        !
        !>@param bf_can_exchange_with_neighbor1
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor1
        !
        !>@param bf_can_exchange_with_neighbor2
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor2
        !
        !>@param bf_newgrdpt_coords
        !> local coordinates of the new grid point to be computed
        !
        !>@return grdpts_available
        !> logical determining whether there are enough grid points
        !> to identify whether there is a bc_interior_pt crenel or not
        !--------------------------------------------------------------
        function are_grdpts_available_to_detect_bc_interior_pt_crenel(
     $       bf_localization,
     $       size_x, size_y,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       bf_newgrdpt_coords)
     $       result(grdpts_available)
        
          implicit none

          integer                     , intent(in) :: bf_localization
          integer                     , intent(in) :: size_x
          integer                     , intent(in) :: size_y
          logical                     , intent(in) :: bf_can_exchange_with_neighbor1
          logical                     , intent(in) :: bf_can_exchange_with_neighbor2
          integer(ikind), dimension(2), intent(in) :: bf_newgrdpt_coords
          logical                                  :: grdpts_available

          
          grdpts_available = are_grdpts_available_to_get_newgrdpt_proc(
     $         bf_localization,
     $         size_x, size_y,
     $         bf_can_exchange_with_neighbor1,
     $         bf_can_exchange_with_neighbor2,
     $         bf_newgrdpt_coords)

        end function are_grdpts_available_to_detect_bc_interior_pt_crenel

      end module bf_bc_interior_pt_crenel_module
      
