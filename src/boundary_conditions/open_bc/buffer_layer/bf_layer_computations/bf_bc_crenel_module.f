      !> @file
      !> module encapsulating useful functions to curb a bc_pt crenel
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating useful functions to curb a bc_pt crenel
      !
      !> @date
      ! 21_01_2015 - initial version       - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_bc_crenel_module

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use parameters_bf_layer, only : 
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind

        implicit none


        private
        public ::
     $       is_temp_array_needed_for_bc_crenel,
     $       detect_bc_double_crenel,
     $       curb_bc_double_crenel,
     $       detect_and_curb_bc_double_crenel,
     $       detect_and_curb_bc_single_crenel,
     $       detect_and_curb_bc_crenels


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether a temporary array is needed to control
        !> the existence of a bc_pt crenel
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_localization
        !> cardinal point identifying the position of the
        !> buffer layer
        !
        !>@param bf_sizes
        !> buffer layer extents
        !
        !>@param cpt_local_coords
        !> coordinates of the central bc_pt checked
        !> using the indices of the buffer layer
        !
        !>@return temp_array_needed
        !> check whether a temporary array is needed
        !--------------------------------------------------------------
        function is_temp_array_needed_for_bc_crenel(
     $       bf_localization,
     $       bf_sizes,
     $       cpt_local_coords)
     $       result(temp_array_needed)

          implicit none

          integer                     , intent(in) :: bf_localization
          integer(ikind), dimension(2), intent(in) :: bf_sizes
          integer(ikind), dimension(2), intent(in) :: cpt_local_coords
          logical                                  :: temp_array_needed


          select case(bf_localization)
            case(N)
               temp_array_needed = ((cpt_local_coords(1)-bc_size).lt.1).or.
     $                             ((cpt_local_coords(1)+bc_size).gt.bf_sizes(1)).or.
     $                             ((cpt_local_coords(2)-bc_size).lt.1)

            case(S)
               temp_array_needed = ((cpt_local_coords(1)-bc_size).lt.1).or.
     $                             ((cpt_local_coords(1)+bc_size).gt.bf_sizes(1)).or.
     $                             ((cpt_local_coords(2)+bc_size).gt.bf_sizes(2))

            case(E)
               temp_array_needed = ((cpt_local_coords(1)-bc_size).lt.1).or.
     $                             ((cpt_local_coords(2)-bc_size).lt.1).or.
     $                             ((cpt_local_coords(2)+bc_size).gt.bf_sizes(2))


            case(W)
               temp_array_needed = ((cpt_local_coords(1)+bc_size).gt.bf_sizes(1)).or.
     $                             ((cpt_local_coords(2)-bc_size).lt.1).or.
     $                             ((cpt_local_coords(2)+bc_size).gt.bf_sizes(2))


            case default
               call error_mainlayer_id(
     $              'bf_grdpts_id_bc_crenel_module',
     $              'is_temp_grdptsid_array_needed_for_bc_crenel',
     $              bf_localization)

          end select

        end function is_temp_array_needed_for_bc_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> detect the position of the bc_pt crenel if any
        !> there are four possible double crenels:
        !>
        !>   _|3     3|_        ___        _______
        !>  |3 3     3 3|     _|3 3|_     |3 3 3 3|
        !>  |3 3     3 3|    |3 3 3 3|      |3 3|
        !>    |3     3|
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param cpt_local_coords
        !> coordinates of the central bc_pt checked
        !> using the indices of the buffer layer
        !
        !>@param bf_sizes
        !> extents of the grdpts_id array
        !
        !>@param bf_grdpts_id
        !> role of grid-points
        !
        !>@param bc_pt_crenel_coords
        !> coordinates of the bc_pt_crenel
        !
        !>@return bc_pt_crenel_exists
        !> check whether a bc_pt crenel exists
        !--------------------------------------------------------------
        function detect_bc_double_crenel(
     $     cpt_local_coords,
     $     bf_sizes,
     $     bf_grdpts_id,
     $     bc_pt_crenel_coords)
     $     result(bc_pt_crenel_exists)

          implicit none

          integer(ikind), dimension(2)  , intent(in)  :: cpt_local_coords
          integer(ikind), dimension(2)  , intent(in)  :: bf_sizes
          integer       , dimension(:,:), intent(in)  :: bf_grdpts_id
          integer(ikind), dimension(2)  , intent(out) :: bc_pt_crenel_coords
          logical                                     :: bc_pt_crenel_exists


          integer(ikind) :: i,j
          logical        :: possible_to_check


          i = cpt_local_coords(1)
          j = cpt_local_coords(2)
          bc_pt_crenel_exists = .false.          


          !there are four possibilities for the bc_pt_crenel:
          !    ___
          !1) |- 0|
          !   |- -|
          !--------------------------------------------------
          possible_to_check = ((i-1).ge.1).and.((j-1).ge.1)

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then
             
             bc_pt_crenel_exists = (bf_grdpts_id(i-1,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i  ,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i-1,j  ).eq.bc_pt)

             bc_pt_crenel_coords = [i-1,j-1]

          end if

          !    ___
          !2) |0 -|
          !   |- -|
          !--------------------------------------------------
          possible_to_check = ((i+1).le.bf_sizes(1)).and.((j-1).ge.1)

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then
             
             bc_pt_crenel_exists = (bf_grdpts_id(i  ,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j  ).eq.bc_pt)

             bc_pt_crenel_coords = [i,j-1]

          end if

          !    ___
          !3) |- -|
          !   |- 0|
          !--------------------------------------------------
          possible_to_check = ((i-1).ge.1).and.((j+1).le.bf_sizes(2))

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then
             
             bc_pt_crenel_exists = (bf_grdpts_id(i-1,j  ).eq.bc_pt).and.
     $                             (bf_grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i  ,j+1).eq.bc_pt)

             bc_pt_crenel_coords = [i-1,j]

          end if

          !    ___
          !4) |- -|
          !   |0 -|
          !--------------------------------------------------
          possible_to_check = ((i+1).le.bf_sizes(1)).and.((j+1).le.bf_sizes(2))

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then
             
             bc_pt_crenel_exists = (bf_grdpts_id(i+1,j  ).eq.bc_pt).and.
     $                             (bf_grdpts_id(i  ,j+1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j+1).eq.bc_pt)

             bc_pt_crenel_coords = [i,j]

          end if

        end function detect_bc_double_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> curb the bc double crenel: from the
        !> identification, the grdpts are updated
        !
        !>   _|3     3|_        ___        _______
        !>  |3 3     3 3|     _|3 3|_     |3 3 3 3|
        !>  |3 3     3 3|    |3 3 3 3|      |3 3|
        !>    |3     3|    
        !>
        !>     |      |          |            |
        !>    \|/    \|/        \|/          \|/
        !>    
        !>    |3|    |3|      _______      _______
        !>    |3|    |3|     |3 3 3 3|    |3 3 3 3|
        !>    |3|    |3|      
        !>    |3|    |3|
        !> 
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_pt_crenel_coords
        !> coordinates of the bc_pt_crenel
        !
        !>@param bf_sizes
        !> extents of the grdpts_id array
        !
        !>@param bf_grdpts_id
        !> role of grid-points
        !
        !>@return bc_pt_crenel_dir
        !> direction of the bc_crenel
        !--------------------------------------------------------------
        subroutine curb_bc_double_crenel(
     $     bc_pt_crenel_coords,
     $     bf_sizes,
     $     bf_grdpts_id)

          implicit none

          integer(ikind), dimension(2)  , intent(in)    :: bc_pt_crenel_coords
          integer(ikind), dimension(2)  , intent(in)    :: bf_sizes
          integer       , dimension(:,:), intent(inout) :: bf_grdpts_id


          integer(ikind) :: i,j,k
          logical        :: possible_to_check
          logical        :: dir_not_found
          logical        :: dir_found

          i = bc_pt_crenel_coords(1)
          j = bc_pt_crenel_coords(2)
          dir_not_found = .true.


          !there are four possibilities:

          !      ___
          !1)  _|- -|_
          !   |-|0 -|-|
          possible_to_check = ((i-1).ge.1).and.((i+2).le.bf_sizes(1))
          if(dir_not_found.and.possible_to_check) then

             dir_found        = (bf_grdpts_id(i-1,j).eq.bc_pt).and.
     $                          (bf_grdpts_id(i+2,j).eq.bc_pt)
             dir_not_found    = .not.dir_found


             if(dir_found) then

                ! 1 1 1 1         1 1 1 1
                ! 2 2 2 2  ___\   1 1 1 1   
                ! 2 3 3 2     /   2 2 2 2
                ! 3 X 3 3         3 3 3 3
                bf_grdpts_id(i  ,j+1) = bc_interior_pt
                bf_grdpts_id(i+1,j+1) = bc_interior_pt
                do k=i-1,i+2
                   bf_grdpts_id(k,j+2) = interior_pt
                end do

             end if

          end if

          !    _ ___ _
          !2) |-|- -|-|
          !     |0 -|
          possible_to_check = ((i-1).ge.1).and.((i+2).le.bf_sizes(1))
          if(dir_not_found.and.possible_to_check) then

             dir_found        = (bf_grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                          (bf_grdpts_id(i+2,j+1).eq.bc_pt)
             dir_not_found    = .not.dir_found


             if(dir_found) then

                ! 3 3 3 3         3 3 3 3
                ! 2 X 3 2  ___\   2 2 2 2   
                ! 2 2 2 2     /   1 1 1 1
                ! 1 1 1 1         1 1 1 1
                do k=i-1,i+2
                   bf_grdpts_id(k,j-1) = interior_pt
                end do
                bf_grdpts_id(i  ,j) = bc_interior_pt
                bf_grdpts_id(i+1,j) = bc_interior_pt

             end if

          end if
           
          !      _
          !     |-|_
          !3)   |- -|
          !     |0 -|
          !     |-|
          possible_to_check = ((j-1).ge.1).and.((j+2).le.bf_sizes(2))
          if(dir_not_found.and.possible_to_check) then

             dir_found        = (bf_grdpts_id(i,j-1).eq.bc_pt).and.
     $                          (bf_grdpts_id(i,j+2).eq.bc_pt)
             dir_not_found    = .not.dir_found


             if(dir_found) then

                ! 3 2 2 1         3 2 1 1
                ! 3 3 2 1  ___\   3 2 1 1   
                ! X 3 2 1     /   X 2 1 1
                ! 3 2 2 1         3 2 1 1
                bf_grdpts_id(i+2,j-1) = interior_pt
                bf_grdpts_id(i+1,j  ) = bc_interior_pt
                bf_grdpts_id(i+2,j  ) = interior_pt
                bf_grdpts_id(i+1,j+1) = bc_interior_pt
                bf_grdpts_id(i+2,j+1) = interior_pt
                bf_grdpts_id(i+2,j+2) = interior_pt

             end if

          end if

          
          !        _
          !      _|-|
          !4)   |- -|
          !     |0 -|
          !       |-|
          possible_to_check = ((j-1).ge.1).and.((j+2).le.bf_sizes(2))
          if(dir_not_found.and.possible_to_check) then

             dir_found        = (bf_grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                          (bf_grdpts_id(i+1,j+2).eq.bc_pt)
             dir_not_found    = .not.dir_found


             if(dir_found) then

                ! 1 2 2 3         1 1 2 3
                ! 1 2 3 3  ___\   1 1 2 3   
                ! 1 2 X 3     /   1 1 2 3
                ! 1 2 2 3         1 1 2 3
                bf_grdpts_id(i-1,j-1) = interior_pt
                bf_grdpts_id(i-1,j  ) = interior_pt
                bf_grdpts_id(i  ,j  ) = bc_interior_pt
                bf_grdpts_id(i-1,j+1) = interior_pt
                bf_grdpts_id(i  ,j+1) = bc_interior_pt
                bf_grdpts_id(i-1,j+2) = interior_pt

             end if

          end if

        end subroutine curb_bc_double_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> detect and curb double crenel, i.e. crenels of
        !> the form:
        !>
        !>   _|3     3|_        ___        _______
        !>  |3 3     3 3|     _|3 3|_     |3 3 3 3|
        !>  |3 3     3 3|    |3 3 3 3|      |3 3|
        !>    |3     3|    
        !>
        !>     |      |          |            |
        !>    \|/    \|/        \|/          \|/
        !>    
        !>    |3|    |3|      _______      _______
        !>    |3|    |3|     |3 3 3 3|    |3 3 3 3|
        !>    |3|    |3|      
        !>    |3|    |3|
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param cpt_local_coords
        !> coordinates of the central bc_pt checked
        !> using the indices of the buffer layer
        !
        !>@param bf_sizes
        !> extents of the grdpts_id array
        !
        !>@param bf_grdpts_id
        !> role of grid-points
        !
        !>@return bc_pt_crenel_exists
        !> check whether a bc_pt crenel exists
        !--------------------------------------------------------------
         function detect_and_curb_bc_double_crenel(
     $     cpt_local_coords,
     $     bf_sizes,
     $     bf_grdpts_id)
     $     result(bc_pt_crenel_exists)

          implicit none

          integer(ikind), dimension(2)  , intent(in)    :: cpt_local_coords
          integer(ikind), dimension(2)  , intent(in)    :: bf_sizes
          integer       , dimension(:,:), intent(inout) :: bf_grdpts_id
          logical                                       :: bc_pt_crenel_exists


          integer(ikind), dimension(2) :: bc_pt_crenel_coords


          !1) detect whether there is a bc_pt crenel
          bc_pt_crenel_exists = detect_bc_double_crenel(
     $         cpt_local_coords,
     $         bf_sizes,
     $         bf_grdpts_id,
     $         bc_pt_crenel_coords)


          !2) curb the bc_pt crenel if there is one
          if(bc_pt_crenel_exists) then

             call curb_bc_double_crenel(
     $            bc_pt_crenel_coords,
     $            bf_sizes,
     $            bf_grdpts_id)

          end if

        end function detect_and_curb_bc_double_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> detect and curb single crenel, i.e. crenels of
        !> the form:
        !>
        !>   _|3     3|_        _          _____  
        !>  |3 3     3 3|     _|3|_       |3 3 3|  
        !>    |3     3|      |3 3 3|        |3|  
        !>
        !>     |      |         |            |
        !>    \|/    \|/       \|/          \|/
        !>
        !>    |3|    |3|      _____        _____
        !>    |3|    |3|     |3 3 3|      |3 3 3|
        !>    |3|    |3|
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param cpt_local_coords
        !> coordinates of the central bc_pt checked
        !> using the indices of the buffer layer
        !
        !>@param bf_sizes
        !> extents of the grdpts_id array
        !
        !>@param bf_grdpts_id
        !> role of grid-points
        !
        !>@return bc_pt_crenel_exists
        !> check whether a bc_pt crenel exists
        !--------------------------------------------------------------
         function detect_and_curb_bc_single_crenel(
     $     cpt_local_coords,
     $     bf_sizes,
     $     bf_grdpts_id)
     $     result(bc_pt_crenel_exists)

          implicit none

          integer(ikind), dimension(2)  , intent(in)    :: cpt_local_coords
          integer(ikind), dimension(2)  , intent(in)    :: bf_sizes
          integer       , dimension(:,:), intent(inout) :: bf_grdpts_id
          logical                                       :: bc_pt_crenel_exists


          integer(ikind) :: i,j
          logical        :: possible_to_check


          i = cpt_local_coords(1)
          j = cpt_local_coords(2)
          bc_pt_crenel_exists = .false.  


          !>1) 
          !>    _|3 
          !>   |X 3  
          !>     |3
          !>----------
          possible_to_check = ((i+1).le.bf_sizes(1)).and.
     $                        ((j-1).ge.1).and.
     $                        ((j+1).le.bf_sizes(2))

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then

             bc_pt_crenel_exists = (bf_grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j  ).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j+1).eq.bc_pt)

             if(bc_pt_crenel_exists) then

                bf_grdpts_id(i,j) = bc_interior_pt

             end if

          end if
          
          !>2)
          !>   3|_ 
          !>   3 X|
          !>   3|  
          !>----------
          possible_to_check = ((i-1).ge.1).and.
     $                        ((j-1).ge.1).and.
     $                        ((j+1).le.bf_sizes(2))

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then

             bc_pt_crenel_exists = (bf_grdpts_id(i-1,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i-1,j  ).eq.bc_pt).and.
     $                             (bf_grdpts_id(i-1,j+1).eq.bc_pt)

             if(bc_pt_crenel_exists) then

                bf_grdpts_id(i,j) = bc_interior_pt

             end if

          end if

          !>3)
          !>      _   
          !>    _|X|_ 
          !>   |3 3 3|
          !>----------
          possible_to_check = ((i-1).ge.1).and.
     $                        ((i+1).le.bf_sizes(1)).and.
     $                        ((j-1).ge.1)

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then

             bc_pt_crenel_exists = (bf_grdpts_id(i-1,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i  ,j-1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j-1).eq.bc_pt)

             if(bc_pt_crenel_exists) then

                bf_grdpts_id(i,j) = bc_interior_pt

             end if

          end if

          !>4)
          !>    _____ 
          !>   |3 3 3|
          !>     |X|
          !>----------
          possible_to_check = ((i-1).ge.1).and.
     $                        ((i+1).le.bf_sizes(1)).and.
     $                        ((j+1).le.bf_sizes(2))

          if(possible_to_check.and.(.not.bc_pt_crenel_exists)) then

             bc_pt_crenel_exists = (bf_grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i  ,j+1).eq.bc_pt).and.
     $                             (bf_grdpts_id(i+1,j+1).eq.bc_pt)

             if(bc_pt_crenel_exists) then

                bf_grdpts_id(i,j) = bc_interior_pt

             end if

          end if          

        end function detect_and_curb_bc_single_crenel


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> detect and curb single and double crenels
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param cpt_local_coords
        !> coordinates of the central bc_pt checked
        !> using the indices of the buffer layer
        !
        !>@param bf_sizes
        !> extents of the grdpts_id array
        !
        !>@param bf_grdpts_id
        !> role of grid-points
        !
        !>@return bc_pt_crenel_exists
        !> check whether a bc_pt crenel exists
        !--------------------------------------------------------------
         function detect_and_curb_bc_crenels(
     $     cpt_local_coords,
     $     bf_sizes,
     $     bf_grdpts_id)
     $     result(bc_pt_crenel_exists)

          implicit none

          integer(ikind), dimension(2)  , intent(in)    :: cpt_local_coords
          integer(ikind), dimension(2)  , intent(in)    :: bf_sizes
          integer       , dimension(:,:), intent(inout) :: bf_grdpts_id
          logical                                       :: bc_pt_crenel_exists


          logical :: bc_pt_single_crenel_exists
          logical :: bc_pt_double_crenel_exists


          bc_pt_double_crenel_exists = detect_and_curb_bc_double_crenel(
     $         cpt_local_coords,
     $         bf_sizes,
     $         bf_grdpts_id)

          bc_pt_single_crenel_exists = detect_and_curb_bc_single_crenel(
     $         cpt_local_coords,
     $         bf_sizes,
     $         bf_grdpts_id)

          bc_pt_crenel_exists = bc_pt_double_crenel_exists.or.
     $                          bc_pt_single_crenel_exists

        end function detect_and_curb_bc_crenels

      end module bf_bc_crenel_module
