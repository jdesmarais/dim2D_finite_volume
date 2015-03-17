      !> @file
      !> detection of bc_pt crenels
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> detection of bc_pt crenels
      !
      !> @date
      ! 18_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_bc_pt_crenel_module
      
        use bf_newgrdpt_dispatch_module, only :
     $     are_grdpts_available_to_get_newgrdpt_proc

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt

        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public :: 
     $       remove_bc_pt_crenel,
     $       check_if_bc_pt_crenel,
     $       are_grdpts_available_to_detect_bc_pt_crenel


        contains


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> remove the bc_pt crenel
        !
        !> @date
        !> 18_03_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> configuration of the grid points around the suspicious
        !> bc_interior_pt
        !
        !> @param i
        !> x-index of the bc_pt checked
        !
        !> @param j
        !> y-index of the bc_pt checked
        !--------------------------------------------------------------        
        subroutine remove_bc_pt_crenel(grdpts_id,i,j)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: i
          integer(ikind)         , intent(in)    :: j

          logical :: remove_crenel

          
          remove_crenel = remove_corner_crenel(grdpts_id,i,j)

          if(.not.remove_crenel) then
             remove_crenel = remove_edge_crenel(grdpts_id,i,j)

             if(.not.remove_crenel) then

                print '(''bf_bc_pt_crenel_module'')'
                print '(''remove_bc_pt_crenel'')'
                print '(''configuration not recognized'')'
                print '(''grdpts_id'')'
                print '(3I2)', grdpts_id(i-1:i+1,j+1)
                print '(3I2)', grdpts_id(i-1:i+1,j)
                print '(3I2)', grdpts_id(i-1:i+1,j-1)
                print '()'
                stop ''

             end if

          end if

        end subroutine remove_bc_pt_crenel


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> remove a corner-like bc_pt crenel
        !   _____     _____
        !>  - 3 3|   |3 3 -   |-           -|
        !>    3 3|   |3 3     |3 3       3 3|
        !>      -|   |-       |3_3_-   -_3_3|
        !
        !> @date
        !> 18_03_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> configuration of the grid points around the suspicious
        !> bc_interior_pt
        !
        !> @param i
        !> x-index of the bc_pt checked
        !
        !> @param j
        !> y-index of the bc_pt checked
        !
        !> @return remove_crenel
        !> confirm that crenel is removed
        !-------------------------------------------------------------- 
        function remove_corner_crenel(grdpts_id,i,j)
     $     result(remove_crenel)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: i
          integer(ikind)         , intent(in)    :: j
          logical                                :: remove_crenel

          
          remove_crenel = .true.
          

          ! -------       ------- 
          !| - - 2 |     | - - 1 |
          !| 3 . - | ->  | 3 2 - |
          !| 3 3 - |     | 3 3 - |
          ! -------       ------- 
          if(
     $         (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $         (grdpts_id(i,j-1).eq.bc_pt).and.
     $         (grdpts_id(i-1,j).eq.bc_pt).and.
     $         (grdpts_id(i+1,j+1).eq.bc_interior_pt)) then

             grdpts_id(i,j)     = bc_interior_pt
             grdpts_id(i+1,j+1) = interior_pt

          else

          ! -------       ------- 
          !| 2 - - |     | 1 - - |
          !| - . 3 | ->  | - 2 3 |
          !| - 3 3 |     | - 3 3 |
          ! -------       ------- 
             if(
     $            (grdpts_id(i,j-1).eq.bc_pt).and.
     $            (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $            (grdpts_id(i+1,j).eq.bc_pt).and.
     $            (grdpts_id(i-1,j+1).eq.bc_interior_pt)) then
                
                grdpts_id(i,j)     = bc_interior_pt
                grdpts_id(i-1,j+1) = interior_pt
                
             else

          ! -------       ------- 
          !| 3 3 - |     | 3 3 - |
          !| 3 . - | ->  | 3 2 - |
          !| - - 2 |     | - - 1 |
          ! -------       ------- 
                if(
     $               (grdpts_id(i+1,j-1).eq.bc_interior_pt).and.
     $               (grdpts_id(i-1,j).eq.bc_pt).and.
     $               (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $               (grdpts_id(i,j+1).eq.bc_pt)) then

                   grdpts_id(i+1,j-1) = interior_pt
                   grdpts_id(i,j)     = bc_interior_pt

                else

          ! -------       ------- 
          !| - 3 3 |     | - 3 3 |
          !| - . 3 | ->  | - 2 3 |
          !| 2 - - |     | 1 - - |
          ! -------       ------- 
                   if(
     $                  (grdpts_id(i-1,j-1).eq.bc_interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.bc_pt).and.
     $                  (grdpts_id(i,j+1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j+1).eq.bc_pt)) then

                      grdpts_id(i-1,j-1) = interior_pt
                      grdpts_id(i,j)     = bc_interior_pt

                   else

                      remove_crenel = .false.

                   end if
                end if
             end if
          end if             

        end function remove_corner_crenel


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> remove a edge-like bc_pt crenel
        !
        !  |3       3|   _____     3
        !  |3 3   3 3|   3 3 3   3_3_3
        !  |3       3|     3
        !
        !
        !> @date
        !> 18_03_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> configuration of the grid points around the suspicious
        !> bc_interior_pt
        !
        !> @param i
        !> x-index of the bc_pt checked
        !
        !> @param j
        !> y-index of the bc_pt checked
        !
        !> @return remove_crenel
        !> confirm that crenel is removed
        !-------------------------------------------------------------- 
        function remove_edge_crenel(grdpts_id,i,j)
     $     result(remove_crenel)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: i
          integer(ikind)         , intent(in)    :: j
          logical                                :: remove_crenel

          
          remove_crenel = .true.

          ! -------       ------- 
          !| 3 - - |     | 3 - 1 |
          !| 3 3 2 | ->  | 3 2 1 |
          !| 3 - - |     | 3 - 1 |
          ! -------       ------- 
          if(
     $         (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $         (grdpts_id(i-1,j).eq.bc_pt).and.
     $         (grdpts_id(i+1,j).eq.bc_interior_pt).and.
     $         (grdpts_id(i-1,j+1).eq.bc_pt)) then
             
             grdpts_id(i+1,j-1) = interior_pt
             grdpts_id(i,j)     = bc_interior_pt
             grdpts_id(i+1,j)   = interior_pt
             grdpts_id(i+1,j+1) = interior_pt

          ! -------       ------- 
          !| - - 3 |     | 1 - 3 |
          !| 2 3 3 | ->  | 1 2 3 |
          !| - - 3 |     | 1 - 3 |
          ! -------       ------- 
          else
             
             if(
     $            (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $            (grdpts_id(i-1,j).eq.bc_interior_pt).and.
     $            (grdpts_id(i+1,j).eq.bc_pt).and.
     $            (grdpts_id(i+1,j+1).eq.bc_pt)) then
             
                grdpts_id(i-1,j-1) = interior_pt
                grdpts_id(i-1,j)   = interior_pt
                grdpts_id(i,j)     = bc_interior_pt
                grdpts_id(i-1,j+1) = interior_pt

          ! -------       ------- 
          !| - 2 - |     | 1 1 1 |
          !| - 3 - | ->  | - 2 - |
          !| 3 3 3 |     | 3 3 3 |
          ! -------       ------- 
             else

                if(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i,j+1).eq.bc_interior_pt)) then
             
                   grdpts_id(i,j)     = bc_interior_pt
                   grdpts_id(i-1,j+1) = interior_pt
                   grdpts_id(i,j+1)   = interior_pt
                   grdpts_id(i+1,j+1) = interior_pt

          ! -------       ------- 
          !| 3 3 3 |     | 3 3 3 |
          !| - 3 - | ->  | - 2 - |
          !| - 2 - |     | 1 1 1 |
          ! -------       ------- 
                else
                   if(
     $                  (grdpts_id(i,j-1).eq.bc_interior_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                  (grdpts_id(i,j+1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j+1).eq.bc_pt)) then
                      
                      grdpts_id(i-1,j-1) = interior_pt
                      grdpts_id(i,j-1)   = interior_pt
                      grdpts_id(i+1,j-1) = interior_pt
                      grdpts_id(i,j)     = bc_interior_pt
                      
                   else

                      remove_crenel = .false.

                   end if
                end if
             end if
          end if

        end function remove_edge_crenel

        
        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the bc_pt leads to a bc_pt crenel that
        !> should be removed (i.e. whether the suspicious bc_pt
        !> has a no_pt neighbor)
        !
        !> @date
        !> 18_03_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> configuration of the grid points around the suspicious
        !> bc_interior_pt
        !
        !> @param i_center
        !> x-index of the bc_pt checked
        !
        !> @param j_center
        !> y-index of the bc_pt checked
        !
        !> @return is_a_bc_pt_crenel
        !> logical determining whether this is a bc_crenel or not
        !--------------------------------------------------------------        
        function check_if_bc_pt_crenel(grdpts_id,i_center,j_center)
     $       result(is_a_bc_pt_crenel)

          implicit none

          integer       , dimension(:,:), intent(in) :: grdpts_id
          integer(ikind)                , intent(in) :: i_center
          integer(ikind)                , intent(in) :: j_center
          logical                                    :: is_a_bc_pt_crenel
          

          integer(ikind) :: i,j
          

          is_a_bc_pt_crenel = .true.


          do j=j_center-1,j_center+1
             do i=i_center-1, i_center+1
                if(grdpts_id(i,j).eq.no_pt) then
                   is_a_bc_pt_crenel = .false.
                   exit
                end if                
             end do
             if(.not.is_a_bc_pt_crenel) exit
          end do

        end function check_if_bc_pt_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are enough grid points to
        !> identify whether there is a bc_pt crenel
        !> or not
        !
        !> @date
        !> 18_03_2015 - initial version - J.L. Desmarais
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
        !> to identify whether there is a bc_pt crenel or not
        !--------------------------------------------------------------
        function are_grdpts_available_to_detect_bc_pt_crenel(
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

        end function are_grdpts_available_to_detect_bc_pt_crenel

      end module bf_bc_pt_crenel_module
