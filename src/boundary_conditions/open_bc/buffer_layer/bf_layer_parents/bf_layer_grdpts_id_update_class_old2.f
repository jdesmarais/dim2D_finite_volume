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

        use bf_bc_crenel_module, only :
     $      detect_and_curb_bc_crenels

        use bf_layer_newgrdpt_class, only :
     $       bf_layer_newgrdpt

        use parameters_bf_layer, only :
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt,
     $       BF_SUCCESS

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        private
        public :: bf_layer_grdpts_id_update


        !> @class bf_layer_grdpts_id_update
        !> bf_layer_newgrdpt augmented with procedures updating the
        !> configuration of the grdpts_id
        !
        !> @param check_grdpts_id_pt
        !  check if the ID of the grid-point is the configuration checked
        !
        !> @param set_grdpts_id_pt
        !  set the the ID of the grid-point to the configuration given
        !
        !> @param check_neighboring_bc_interior_pts
        !> verify whether there are bc_interior_pt to be updated
        !
        !> @param check_bc_interior_pt
        !> verify if it is a bc_interior_pt
        !
        !> @param set_new_interior_pt
        !> set the ID of the grid-point as interior_pt
        !
        !> @param finalize_neighboring_grdpts_update
        !> update the grid points around the new interior point
        !
        !> @param update_neighboring_grdpt
        !> check whether the grid point asked exists or not.
        !> If it does not exist in the buffer layer, it has
        !> to be computed and updated in the gridpoint ID map
        !
        !> @param verify_if_all_bf_grdpts_exist
        !> verify if all the grdpts exists around a grdpt
        !> to turn it into an interior_pt
        !
        !> @param has_a_bc_pt_neighbor
        !> verify if there is a bc_pt in the neighboring of the
        !> grid points
        !
        !> @param detect_and_curb_bc_pt_crenels
        !> remove the crenels of grid points of type bc_pt
        !-------------------------------------------------------------
        type, extends(bf_layer_newgrdpt) :: bf_layer_grdpts_id_update
        
          contains

          ! local modification of the grdpts_id
          procedure,   pass :: check_grdpts_id_pt
          procedure,   pass :: set_grdpts_id_pt

          !for updating the bc_interior_pt and bc_pt around a new
          !interior_pt          
          procedure,   pass :: update_neighboring_grdpt
          procedure,   pass :: finalize_neighboring_grdpts_update

          ! for detecting bc_interior_pt to be turned into interior_pt
          procedure,   pass :: check_neighboring_bc_interior_pts
          procedure, nopass :: check_bc_interior_pt          

          !for suspicious bf_interior_pt checks
          procedure,   pass :: verify_if_all_bf_grdpts_exist
          procedure,   pass :: has_a_bc_pt_neighbor
                       
          !for bc_pt crenels checks
          procedure,   pass :: detect_and_curb_bc_pt_crenels

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
        !> set the configuration of the grdpts(i,j)to config_set
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
        !> @param config_set
        !> configuration set to the grdpts_id(i,j)
        !--------------------------------------------------------------
        subroutine set_grdpts_id_pt(this, i,j, config_set)

          implicit none

          class(bf_layer_grdpts_id_update), intent(inout) :: this
          integer(ikind)                  , intent(in)    :: i
          integer(ikind)                  , intent(in)    :: j
          integer                         , intent(in)    :: config_set

          this%grdpts_id(i,j) = config_set

        end subroutine set_grdpts_id_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the grid points neighboring a point identified by its
        !> general coordinates (cpt_coords) are bc_interior_pt, if so,
        !> the points are added to a list of bc_interior_pt
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i_prev
        !> integer identifying the x-index of the previous central point
        !> investigated
        !
        !>@param j_prev
        !> integer identifying the y-index of the previous central point
        !> investigated
        !
        !>@param i_center
        !> integer identifying the x-index of the central point
        !> investigated
        !
        !>@param j_central
        !> integer identifying the y-index of the central point
        !> investigated
        !
        !>@param nb_mgrdpts
        !> number of grid points in the list
        !>
        !>@param mgrdpts
        !> list of indices identifying the bc_interior_pt surrounding
        !> the central grid point
        !--------------------------------------------------------------
        subroutine check_neighboring_bc_interior_pts(
     $     this,
     $     i_prev, j_prev,
     $     i_center, j_center,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          class(bf_layer_grdpts_id_update), intent(in)    :: this
          integer(ikind)                  , intent(in)    :: i_prev
          integer(ikind)                  , intent(in)    :: j_prev
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center
          integer                         , intent(inout) :: nb_mgrdpts
          integer(ikind) , dimension(:,:) , intent(out)   :: mgrdpts


          !radius for the search of bc_interior_pt around the
          !central point identified by (i_center, j_center)
          integer, parameter :: search_r = 1

          integer(ikind), dimension(2) :: match_table
          integer(ikind) :: min_j, max_j
          integer(ikind) :: size_x, size_y
          integer(ikind) :: i,j

          !get the match table converting the general coords
          !into local coords
          match_table = this%get_general_to_local_coord_tab()

          !get the borders of the loops
          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          size_x = size(this%grdpts_id,1)
          size_y = size(this%grdpts_id,2)


          !1.
          do j=max(1,j_center-search_r-match_table(2)),
     $         min(size_y, j_center+search_r-match_table(2), j_prev-search_r-1-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x, i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do

          
          !2.
          do j=max(1,j_center-search_r-min_j-match_table(2)),
     $         min(size_y, j_center+search_r-max_j-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1), i_prev-search_r-1-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)

             end do
          end do


          !3.
          do j=max(1,j_center-search_r-min_j-match_table(2)),
     $         min(size_y,j_center+search_r-max_j-match_table(2))

             do i=max(1,i_center-search_r-match_table(1),i_prev+search_r+1-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do


          !4.
          do j=max(1,j_center-search_r-match_table(2),j_prev+search_r+1-match_table(2)),
     $         min(size_y,j_center+search_r-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do

        end subroutine check_neighboring_bc_interior_pts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid point tested is a bc_interior_pt
        !> and if so save the general coordinates of the grid point
        !> in mgrdpts
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param i
        !> integer identifying the x-index of the central point
        !> investigated
        !
        !>@param j
        !> integer identifying the y-index of the central point
        !> investigated
        !
        !>@param match_table
        !> table converting local coordinate into general coordinates
        !
        !>@param grdpts_id
        !> table where the role of the grid points is stored
        !
        !>@param nb_mgrdpts
        !> number of grid points in the list
        !>
        !>@param mgrdpts
        !> list of indices identifying the bc_interior_pt surrounding
        !> the central grid point
        !--------------------------------------------------------------
        subroutine check_bc_interior_pt(
     $     i,j,
     $     match_table,
     $     grdpts_id,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          integer(ikind)                , intent(in)    :: i
          integer(ikind)                , intent(in)    :: j
          integer(ikind), dimension(2)  , intent(in)    :: match_table
          integer       , dimension(:,:), intent(in)    :: grdpts_id
          integer                       , intent(inout) :: nb_mgrdpts
          integer(ikind), dimension(:,:), intent(out)   :: mgrdpts

          if(grdpts_id(i,j).eq.bc_interior_pt) then

             nb_mgrdpts = nb_mgrdpts+1
             mgrdpts(1,nb_mgrdpts) = i+match_table(1)
             mgrdpts(2,nb_mgrdpts) = j+match_table(2)
             
          end if

        end subroutine check_bc_interior_pt        


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> set the grid point as interior_pt
        !>
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !>
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> integer identifying the x-index of the
        !> updated grid point
        !
        !>@param j
        !> integer identifying the y-index of the
        !> updated grid point
        !--------------------------------------------------------------
        subroutine set_new_interior_pt(this,i,j)

          implicit none

          class(bf_layer_grdpts_id_update), intent(inout) :: this
          integer(ikind)                  , intent(in)    :: i
          integer(ikind)                  , intent(in)    :: j

          this%grdpts_id(i,j) = interior_pt

        end subroutine set_new_interior_pt


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
        subroutine finalize_neighboring_grdpts_update(
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

        end subroutine finalize_neighboring_grdpts_update


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid point asked exists or not.
        !> If it does not exist in the buffer layer, it has
        !> to be computed and updated in the gridpoint ID map
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
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
        function update_neighboring_grdpt(this,i,j,i_center,j_center)
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

        end function update_neighboring_grdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> verify if all the grid points around a central grid point
        !> exist in order to be turned into an interior_pt
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param 
        !> logical identifying whether all grid points around a
        !> central grid point exist in order to be turned into
        !> an interior_pt
        !--------------------------------------------------------------
        function verify_if_all_bf_grdpts_exist(
     $     this,i_c,j_c,ierror)
     $     result(all_grdpts_exist)

          implicit none

          class(bf_layer_grdpts_id_update), intent(in)  :: this
          integer(ikind)                  , intent(in)  :: i_c
          integer(ikind)                  , intent(in)  :: j_c
          logical                         , intent(out) :: ierror
          logical                                       :: all_grdpts_exist

          if(
     $         (i_c.ge.(bc_size+1)).and.
     $         (i_c.le.size(this%grdpts_id,1)-bc_size).and.
     $         (j_c.ge.(bc_size+1)).and.
     $         (j_c.le.size(this%grdpts_id,2)-bc_size)) then

             ierror = BF_SUCCESS

             all_grdpts_exist = verify_if_all_grdpts_exist(
     $            i_c, j_c, this%grdpts_id)
             
          else

             ierror = .not.BF_SUCCESS
             all_grdpts_exist = .false.

          end if

        end function verify_if_all_bf_grdpts_exist


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if there is a bc_pt in the neighboring grid
        !> points
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param local_coords 
        !> local coordinates of the grid point whose neighbors are
        !> tested
        !
        !> @return bc_pt_neighbor
        !> logical stating whether there is a neighboring grid point 
        !> which is a bc_pt
        !--------------------------------------------------------------
        function has_a_bc_pt_neighbor(
     $     this, local_coords, ierror)
     $     result(bc_pt_neighbor)

          implicit none

          class(bf_layer_grdpts_id_update), intent(in) :: this
          integer(ikind) , dimension(2)   , intent(in) :: local_coords
          logical                         , intent(out):: ierror
          logical                                      :: bc_pt_neighbor

          integer :: i,j


          if(
     $         (local_coords(1).gt.1).and.
     $         (local_coords(1).lt.size(this%grdpts_id,1)).and.
     $         (local_coords(2).gt.1).and.
     $         (local_coords(2).lt.size(this%grdpts_id,2))) then

             ierror = BF_SUCCESS
             bc_pt_neighbor = .false.
             
             do j=local_coords(2)-1,local_coords(2)+1
                do i=local_coords(1)-1,local_coords(1)+1
                   
                   if(this%grdpts_id(i,j).eq.bc_pt) then
                      bc_pt_neighbor=.true.
                      exit
                   end if
                   
                end do
             end do

          else

             ierror = .not.BF_SUCCESS
             bc_pt_neighbor = .false.

          end if

        end function has_a_bc_pt_neighbor        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if there is a bc_pt in the neighboring grid
        !> points
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param local_coords 
        !> local coordinates of the grid point whose neighbors are
        !> tested
        !
        !> @return bc_pt_neighbor
        !> logical stating whether there is a neighboring grid point 
        !> which is a bc_pt
        !--------------------------------------------------------------
        function detect_and_curb_bc_pt_crenels(
     $     this,
     $     cpt_local_coords)
     $     result(bc_pt_crenel_exists)

          implicit none

          class(bf_layer_grdpts_id_update), intent(inout) :: this
          integer(ikind), dimension(2)    , intent(in)    :: cpt_local_coords
          logical                                         :: bc_pt_crenel_exists

          
          bc_pt_crenel_exists = detect_and_curb_bc_crenels(
     $         cpt_local_coords,
     $         [size(this%grdpts_id,1),size(this%grdpts_id,2)],
     $         this%grdpts_id)


        end function detect_and_curb_bc_pt_crenels

      end module bf_layer_grdpts_id_update_class
