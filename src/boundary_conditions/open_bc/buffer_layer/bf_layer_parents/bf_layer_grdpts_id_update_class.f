      !> @file
      !> bf_layer_newgrdpt augmented with procedures updating the
      !> configuration of the grdpts_id
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_newgrdpt augmented with procedures updating the
      !> configuration of the grdpts_id
      !
      !> @date
      ! 17_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_grdpts_id_update_class

        use bf_bc_interior_pt_crenel_module, only :
     $     check_if_bc_interior_pt_crenel,
     $     are_grdpts_available_to_detect_bc_interior_pt_crenel

        use bf_layer_extract_module, only :
     $     get_bf_layer_match_table,
     $     get_grdpts_id_from_interior,
     $     get_grdpts_id_from_buffer_layer

        use bf_layer_newgrdpt_class, only :
     $       bf_layer_newgrdpt

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind

        private
        public :: bf_layer_grdpts_id_update


        !> @class bf_layer_grdpts_id_update
        !> bf_layer_newgrdpt augmented with procedures updating the
        !> configuration of the grdpts_id
        !
        !> @param check_grdpts_id_pt
        !> check if the ID of the grid-point is the configuration checked
        !
        !> @param set_grdpts_id_pt
        !> set the the ID of the grid-point to the configuration given
        !
        !> @param update_grdpts_id_next_to_new_interior_pt
        !> 
        !-------------------------------------------------------------
        type, extends(bf_layer_newgrdpt) :: bf_layer_grdpts_id_update
        
          contains

          ! local check of the grdpts_id
          procedure, pass :: check_grdpts_id_pt

          ! for updating the bc_interior_pt and bc_pt
          ! around a new interior_pt
          procedure, pass :: update_grdpt_id_next_to_new_interior_pt
          procedure, pass :: finalize_update_grdpts_id_around_new_interior_pt

          ! for detecting and curbing bc_interior_pt crenels
          procedure, pass :: detect_bc_interior_pt_crenel
          


        end type bf_layer_grdpts_id_update


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the grdpts(i,j) corresponds to the configuration
        !> config_checked
        !
        !> @date
        !> 17_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param i
        !> x-index of the grdpts_id tested
        !
        !> @param j
        !> y-index of the grdpts_id tested
        !
        !> @param config_checked
        !> configuration checked
        !
        !> @return check_grdpts_id_pt
        !> logical stating whether the grdpts_id(i,j) corresponds
        !> to the configuration checked, config_checked
        !--------------------------------------------------------------
        function check_grdpts_id_pt(this, i,j, config_checked)

          implicit none

          class(bf_layer_grdpts_id_update), intent(in) :: this
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          integer                         , intent(in) :: config_checked
          logical                                      :: check_grdpts_id_pt

          check_grdpts_id_pt = this%grdpts_id(i,j).eq.config_checked

        end function check_grdpts_id_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update a grid point next to a new interior_pt such
        !> that:
        !> 
        !>   - interior_pt    -> do nothing
        !> 
        !>   - bc_interior_pt -> do_nothing
        !> 
        !>   - bc_pt          -> bc_interior_pt or bc_pt
        !>                       (depending on the distance
        !>                        to the new interior_pt)
        !>
        !>   - no_pt          -> bc_pt or bc_interior_pt +
        !>                       ask for computating the new
        !>                       grid-point
        !>                       (depending on the distance
        !>                        to the new interior_pt)
        !
        !> @date
        !> 17_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> x-index of the grid points checked
        !
        !>@param j
        !> y-index of the grid points checked
        !
        !>@param i_center
        !> x-index of the central point
        !
        !>@param j_center
        !> y-index of the central point 
        !
        !> @return compute_newgrdpt
        !> logical identifying whether the grid point should be
        !> computed
        !---------------------------------------------------------------
        function update_grdpt_id_next_to_new_interior_pt(this,i,j,i_center,j_center)
     $     result(compute_newgrdpt)

          implicit none

          class(bf_layer_grdpts_id_update), intent(inout) :: this
          integer(ikind)                  , intent(in)    :: i
          integer(ikind)                  , intent(in)    :: j
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center
          logical                                         :: compute_newgrdpt

          
          compute_newgrdpt = .false.

          if(
     $         (this%grdpts_id(i,j).eq.no_pt).or.
     $         (this%grdpts_id(i,j).eq.bc_pt)) then

             if (this%grdpts_id(i,j).eq.no_pt) then
                compute_newgrdpt = .true.
             end if

             !if the grid point is next to the new interior point,
             !it is a bc_interior_pt, otherwise, it is a bc_pt
             if((abs(i_center-i).le.1).and.(abs(j_center-j).le.1)) then
                this%grdpts_id(i,j) = bc_interior_pt
             else
                this%grdpts_id(i,j) = bc_pt
             end if

          end if        

        end function update_grdpt_id_next_to_new_interior_pt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> as the grid point is set as new interior_pt, the neighboring
        !> gridpoints should be updated: here only the identity of the
        !> neighboring grid points is updated (bc_pt-> bc_interior_pt)
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i_prev
        !> integer identifying the x-index of the
        !> previous central grid point that was
        !> updated to interior_pt
        !
        !>@param j_prev
        !> integer identifying the y-index of the
        !> previous central grid point that was
        !> updated to interior_pt
        !
        !>@param i_center
        !> integer identifying the x-index of the
        !> central grid point which is updated to
        !> interior_pt
        !
        !>@param j_center
        !> integer identifying the y-index of the
        !> central grid point which is updated to
        !> interior_pt
        !--------------------------------------------------------------
        subroutine finalize_update_grdpts_id_around_new_interior_pt(
     $     this,
     $     i_prev, j_prev,
     $     i_center, j_center)

          implicit none

          class(bf_layer_grdpts_id_update), intent(inout) :: this
          integer(ikind)                  , intent(in)    :: i_prev
          integer(ikind)                  , intent(in)    :: j_prev
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center

          
          integer(ikind) :: i,j


          do j=max(j_prev-bc_size, j_center-1), min(j_prev-bc_size, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev-bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev+bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev+bc_size, j_center-1), min(j_prev+bc_size, j_center+1)
             do i=max(i_center-1, i_prev-bc_size), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

        end subroutine finalize_update_grdpts_id_around_new_interior_pt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the bc_interior_pt leads to a bc_interior_pt
        !> crenel that should be removed
        !
        !> @date
        !> 17_03_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param i
        !> x-index of the bc_interior_pt checked
        !
        !> @param j
        !> y-index of the bc_interior_pt checked
        !
        !> @param ierror
        !> logical stating whether the detection was successful or not
        !
        !> @return is_bc_interior_crenel
        !> logical determining whether this is a bc_crenel or not
        !--------------------------------------------------------------
        function detect_bc_interior_pt_crenel(this,i,j,ierror)
     $     result(is_bc_interior_crenel)

          implicit none

          class(bf_layer_grdpts_id_update), intent(in)  :: this
          integer(ikind)                  , intent(in)  :: i
          integer(ikind)                  , intent(in)  :: j
          logical                         , intent(out) :: ierror
          logical                                       :: is_bc_interior_crenel


          logical                        :: grdpts_available
          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2)   :: gen_coords
          integer       , dimension(3,3) :: grdpts_id_tmp


          !1) check whether there are enough grid points to check whether
          !   there is a bc_interior_pt crenel or not
          !   if there are enough grid points, the presence of a bc_interior_pt
          !   crenel is directly checked on the data of the buffer layer
          if(
     $         ((i-1).le.1).and.
     $         ((i+1).le.(size(this%grdpts_id,1))).and.
     $         ((j-1).le.1).and.
     $         ((j+1).le.(size(this%grdpts_id,2)))) then

             ierror = BF_SUCCESS

             is_bc_interior_crenel = check_if_bc_interior_pt_crenel(
     $            this%grdpts_id,
     $            i,j)

          !2) otherwise, we check whether it is possible to create a temporary
          !   array gathering data from the buffer layer
          else

             grdpts_available = are_grdpts_available_to_detect_bc_interior_pt_crenel(
     $            this%localization,
     $            size(this%grdpts_id,1), size(this%grdpts_id,2),
     $            this%can_exchange_with_neighbor1(),
     $            this%can_exchange_with_neighbor2(),
     $            [i,j])

             if(grdpts_available) then

                ierror = BF_SUCCESS
                
                !> determine the extents of the grdpts_id temporary array
                !> to be created
                match_table = get_bf_layer_match_table(this%alignment)
                
                gen_coords(1,1) = i-1 + match_table(1)
                gen_coords(1,2) = i+1 + match_table(1)
                gen_coords(2,1) = j-1 + match_table(2)
                gen_coords(2,2) = j+1 + match_table(2)
                

                !> extract the grdpts_id to determine the presence
                !> of the bc_interior_pt or not

                !> extract the grdpts_id from the interior domain
                call get_grdpts_id_from_interior(
     $               grdpts_id_tmp,
     $               gen_coords)

                !> extract the grdpts_id from the current buffer layer
                call get_grdpts_id_from_bf_layer(
     $               grdpts_id_tmp,
     $               gen_coords,
     $               this%alignment,
     $               this%grdpts_id)

                !> determine whether there is a bc_interior_pt
                !> crenel or not
                is_bc_interior_crenel = check_if_bc_interior_pt_crenel(
     $               grdpts_id_tmp,
     $               2,2)

             else
                ierror = .not.BF_SUCCESS
             end if

          end if

        end function detect_bc_interior_pt_crenel

      end module bf_layer_grdpts_id_update_class
