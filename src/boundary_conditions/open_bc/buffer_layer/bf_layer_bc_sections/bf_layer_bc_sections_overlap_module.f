      !> @file
      !> module encapsulating the subroutines for the determination of
      !> overlapping of the bc_sections by the integration borders
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines for the determination of
      !> overlapping of the bc_sections by the integration borders
      !
      !> @date
      !> 02_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_bc_sections_overlap_module

        use bf_layer_errors_module, only :
     $       error_overlap_index,
     $       error_cpt_overlap_index

        use parameters_bf_layer, only :
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       
     $       cpt1normal_and_cpt4normal,   
     $       cpt1normal_and_cpt4not,  
     $       cpt1normal_and_cpt4overlap,
     $       cpt1not_and_cpt4normal,
     $       cpt1not_and_cpt4not,
     $       cpt1not_and_cpt4overlap,
     $       cpt1overlap_and_cpt4normal,
     $       cpt1overlap_and_cpt4not,
     $       cpt1overlap_and_cpt4overlap,
     $       
     $       cpt2normal_and_cpt3normal,
     $       cpt2normal_and_cpt3not,  
     $       cpt2normal_and_cpt3overlap,
     $       cpt2not_and_cpt3normal, 
     $       cpt2not_and_cpt3not,     
     $       cpt2not_and_cpt3overlap,
     $       cpt2overlap_and_cpt3normal,
     $       cpt2overlap_and_cpt3not,
     $       cpt2overlap_and_cpt3overlap,
     $       
     $       cptnormal_type,
     $       cptnot_type,
     $       cptoverlap_type


        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public ::
     $       overlap_square_bc_sections,
     $       overlap_corner_and_anti_corner,
     $       overlap_corners_or_anti_corners,
     $       
     $       add_compute_pt1_overlap,
     $       add_compute_pt2_overlap,
     $       add_compute_pt3_overlap,
     $       add_compute_pt4_overlap,
     $       remove_compute_pt1,
     $       remove_compute_pt2,
     $       remove_compute_pt3,
     $       remove_compute_pt4,
     $       is_an_anti_corner,
     $       is_a_corner,
     $       
     $       overlap_bc_section_by_integration_borders,
     $       overlap_N,
     $       overlap_S,
     $       overlap_E,
     $       overlap_W,
     $       
     $       determine_corner_or_anti_corner_grdpts_computed,
     $       determine_edge_grdpts_computed


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the square bc_section
        !> to remove the computation of the grid-points in common
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param square1
        !> boundary layer represented as [square_type,i_min,j_min,cpt_overlap,overlap]
        !
        !> @param square2
        !> boundary layer represented as [square_type,i_min,j_min,cpt_overlap,overlap]
        !--------------------------------------------------------------
        subroutine overlap_square_bc_sections(square1,square2)

          implicit none

          integer, dimension(5), intent(inout) :: square1
          integer, dimension(5), intent(inout) :: square2


          if(is_a_corner(square1)) then
             if(is_a_corner(square2)) then
                call overlap_corners_or_anti_corners(square1,square2)
             else
                call overlap_corner_and_anti_corner(square1,square2)
             end if

          else
             if(is_a_corner(square2)) then
                call overlap_corner_and_anti_corner(square2,square1)
             else
                call overlap_corners_or_anti_corners(square1,square2)
             end if
          end if          

        end subroutine overlap_square_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of the grid-points in common
        !> with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param corner
        !> boundary layer represented as [corner_type,i_min,j_min,cpt_overlap,overlap]
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,cpt_overlap,overlap]
        !--------------------------------------------------------------
        subroutine overlap_corner_and_anti_corner(corner,anti_corner)

          implicit none

          integer, dimension(5), intent(inout) :: corner
          integer, dimension(5), intent(inout) :: anti_corner
          

          if(corner(3).eq.anti_corner(3)) then

          !i_corner = i_anti_corner+1 AND j_corner = j_anti_corner
             if(corner(2).eq.(anti_corner(2)+1)) then
                call overlap_E(anti_corner(5))
             else
             
          !i_corner = i_anti_corner-1 AND j_corner = j_anti_corner
                if(corner(2).eq.(anti_corner(2)-1)) then
                   call overlap_W(anti_corner(5))
                end if
             end if
          end if


          if(corner(2).eq.anti_corner(2)) then

          !j_corner = j_anti_corner+1 AND i_corner = i_anti_corner
             if(corner(3).eq.(anti_corner(3)+1)) then
                call overlap_N(anti_corner(5))
             else
             
          !j_corner = j_anti_corner-1 AND i_corner = i_anti_corner
                if(corner(3).eq.(anti_corner(3)-1)) then
                   call overlap_S(anti_corner(5))
                end if
             end if
          end if

        end subroutine overlap_corner_and_anti_corner
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the square boundary layer
        !> to remove the computation of the grid-points in common
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param square1
        !> boundary layer represented as [square_type,i_min,j_min,cpt_overlap,overlap]
        !
        !> @param square2
        !> boundary layer represented as [square_type,i_min,j_min,cpt_overlap,overlap]
        !--------------------------------------------------------------
        subroutine overlap_corners_or_anti_corners(square1,square2)

          implicit none

          integer, dimension(5), intent(inout) :: square1
          integer, dimension(5), intent(inout) :: square2


          !      _____ square1
          !     |     |
          !   -----.  |
          !  |  |__|__|
          !  |     |
          !   -----  square2
          !
          if((square2(2).eq.(square1(2)-1)).and.(square2(3).eq.(square1(3)-1))) then
             
             call add_compute_pt1_overlap(square1(4))
             call remove_compute_pt4(square2(4))

          end if

          !      _____ square1
          !     |     |
          !     |  .-----
          !     |__|__|  |
          !        |     |
          !         ----- square2
          !
          if((square2(2).eq.(square1(2)+1)).and.(square2(3).eq.(square1(3)-1))) then
             
             call add_compute_pt2_overlap(square1(4))
             call remove_compute_pt3(square2(4))

          end if

          !      square2
          !   -----
          !  |   __|__ 
          !  |  |  |  |
          !   -----.  | square1
          !     |_____|
          !
          if((square2(2).eq.(square1(2)-1)).and.(square2(3).eq.(square1(3)+1))) then
             
             call add_compute_pt3_overlap(square1(4))
             call remove_compute_pt2(square2(4))

          end if

          !         -----
          !      __|__   | square2
          !     |  |  |  |
          !     |   ----- 
          !     |_____|
          !            square1
          !
          if((square2(2).eq.(square1(2)+1)).and.(square2(3).eq.(square1(3)+1))) then
             
             call add_compute_pt4_overlap(square1(4))
             call remove_compute_pt1(square2(4))

          end if

        end subroutine overlap_corners_or_anti_corners  

      
        subroutine add_compute_pt1_overlap(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt1normal_and_cpt4normal)
               cpt_overlap = cpt1overlap_and_cpt4normal

            case(cpt1normal_and_cpt4not)
               cpt_overlap = cpt1overlap_and_cpt4not

            case(cpt1normal_and_cpt4overlap)
               cpt_overlap = cpt1overlap_and_cpt4overlap

            case(cpt1not_and_cpt4normal)
               cpt_overlap = cpt1overlap_and_cpt4normal

            case(cpt1not_and_cpt4not)
               cpt_overlap = cpt1overlap_and_cpt4not

            case(cpt1not_and_cpt4overlap)
               cpt_overlap = cpt1overlap_and_cpt4overlap               
               
            case(cpt1overlap_and_cpt4normal,
     $           cpt1overlap_and_cpt4not,
     $           cpt1overlap_and_cpt4overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine add_compute_pt1_overlap


        subroutine remove_compute_pt1(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt1normal_and_cpt4normal)
               cpt_overlap = cpt1not_and_cpt4normal

            case(cpt1normal_and_cpt4not)
               cpt_overlap = cpt1not_and_cpt4not

            case(cpt1normal_and_cpt4overlap)
               cpt_overlap = cpt1not_and_cpt4overlap

            case(cpt1overlap_and_cpt4normal)
               cpt_overlap = cpt1not_and_cpt4normal

            case(cpt1overlap_and_cpt4not)
               cpt_overlap = cpt1not_and_cpt4not

            case(cpt1overlap_and_cpt4overlap)
               cpt_overlap = cpt1not_and_cpt4overlap               
               
            case(cpt1not_and_cpt4normal,
     $           cpt1not_and_cpt4not,
     $           cpt1not_and_cpt4overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine remove_compute_pt1


        subroutine add_compute_pt2_overlap(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt2normal_and_cpt3normal)
               cpt_overlap = cpt2overlap_and_cpt3normal

            case(cpt2normal_and_cpt3not)
               cpt_overlap = cpt2overlap_and_cpt3not

            case(cpt2normal_and_cpt3overlap)
               cpt_overlap = cpt2overlap_and_cpt3overlap

            case(cpt2not_and_cpt3normal)
               cpt_overlap = cpt2overlap_and_cpt3normal

            case(cpt2not_and_cpt3not)
               cpt_overlap = cpt2overlap_and_cpt3not

            case(cpt2not_and_cpt3overlap)
               cpt_overlap = cpt2overlap_and_cpt3overlap
               
            case(cpt2overlap_and_cpt3normal,
     $           cpt2overlap_and_cpt3not,
     $           cpt2overlap_and_cpt3overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine add_compute_pt2_overlap


        subroutine remove_compute_pt2(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt2normal_and_cpt3normal)
               cpt_overlap = cpt2not_and_cpt3normal

            case(cpt2normal_and_cpt3not)
               cpt_overlap = cpt2not_and_cpt3not

            case(cpt2normal_and_cpt3overlap)
               cpt_overlap = cpt2not_and_cpt3overlap

            case(cpt2overlap_and_cpt3normal)
               cpt_overlap = cpt2not_and_cpt3normal

            case(cpt2overlap_and_cpt3not)
               cpt_overlap = cpt2not_and_cpt3not

            case(cpt2overlap_and_cpt3overlap)
               cpt_overlap = cpt2not_and_cpt3overlap
               
            case(cpt2not_and_cpt3normal,
     $           cpt2not_and_cpt3not,
     $           cpt2not_and_cpt3overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine remove_compute_pt2


        subroutine add_compute_pt3_overlap(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

           select case(cpt_overlap)
            case(cpt2normal_and_cpt3normal)
               cpt_overlap = cpt2normal_and_cpt3overlap

            case(cpt2not_and_cpt3normal)
               cpt_overlap = cpt2not_and_cpt3overlap

            case(cpt2overlap_and_cpt3normal)
               cpt_overlap = cpt2overlap_and_cpt3overlap

            case(cpt2normal_and_cpt3not)
               cpt_overlap = cpt2normal_and_cpt3overlap

            case(cpt2not_and_cpt3not)
               cpt_overlap = cpt2not_and_cpt3overlap

            case(cpt2overlap_and_cpt3not)
               cpt_overlap = cpt2overlap_and_cpt3overlap
               
            case(cpt2normal_and_cpt3overlap,
     $           cpt2not_and_cpt3overlap,
     $           cpt2overlap_and_cpt3overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine add_compute_pt3_overlap


        subroutine remove_compute_pt3(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

           select case(cpt_overlap)
            case(cpt2normal_and_cpt3normal)
               cpt_overlap = cpt2normal_and_cpt3not

            case(cpt2not_and_cpt3normal)
               cpt_overlap = cpt2not_and_cpt3not

            case(cpt2overlap_and_cpt3normal)
               cpt_overlap = cpt2overlap_and_cpt3not

            case(cpt2normal_and_cpt3overlap)
               cpt_overlap = cpt2normal_and_cpt3not

            case(cpt2not_and_cpt3overlap)
               cpt_overlap = cpt2not_and_cpt3not

            case(cpt2overlap_and_cpt3overlap)
               cpt_overlap = cpt2overlap_and_cpt3not
               
            case(cpt2normal_and_cpt3not,
     $           cpt2not_and_cpt3not,
     $           cpt2overlap_and_cpt3not)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine remove_compute_pt3


        subroutine add_compute_pt4_overlap(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt1normal_and_cpt4normal)
               cpt_overlap = cpt1normal_and_cpt4overlap

            case(cpt1not_and_cpt4normal)
               cpt_overlap = cpt1not_and_cpt4overlap

            case(cpt1overlap_and_cpt4normal)
               cpt_overlap = cpt1overlap_and_cpt4overlap

            case(cpt1normal_and_cpt4not)
               cpt_overlap = cpt1normal_and_cpt4overlap

            case(cpt1not_and_cpt4not)
               cpt_overlap = cpt1not_and_cpt4overlap

            case(cpt1overlap_and_cpt4not)
               cpt_overlap = cpt1overlap_and_cpt4overlap
               
            case(cpt1normal_and_cpt4overlap,
     $           cpt1not_and_cpt4overlap,
     $           cpt1overlap_and_cpt4overlap)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine add_compute_pt4_overlap


        subroutine remove_compute_pt4(cpt_overlap)

          implicit none

          integer, intent(inout) :: cpt_overlap

          select case(cpt_overlap)
            case(cpt1normal_and_cpt4normal)
               cpt_overlap = cpt1normal_and_cpt4not

            case(cpt1not_and_cpt4normal)
               cpt_overlap = cpt1not_and_cpt4not

            case(cpt1overlap_and_cpt4normal)
               cpt_overlap = cpt1overlap_and_cpt4not

            case(cpt1normal_and_cpt4overlap)
               cpt_overlap = cpt1normal_and_cpt4not

            case(cpt1not_and_cpt4overlap)
               cpt_overlap = cpt1not_and_cpt4not

            case(cpt1overlap_and_cpt4overlap)
               cpt_overlap = cpt1overlap_and_cpt4not
               
            case(cpt1normal_and_cpt4not,
     $           cpt1not_and_cpt4not,
     $           cpt1overlap_and_cpt4not)

               cpt_overlap = cpt_overlap
               
            case default
               call error_cpt_overlap_index(
     $              'bf_layer_bc_sections_overlap_module',
     $              'add_compute_pt1_overlap',
     $              cpt_overlap)

          end select

        end subroutine remove_compute_pt4


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the boundary layer is an anti-corner
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary layer represented as [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        function is_an_anti_corner(bc_section)

          implicit none

          integer, dimension(4), intent(in) :: bc_section
          logical                           :: is_an_anti_corner

          integer :: procedure_type

          procedure_type = bc_section(1)

          is_an_anti_corner = (procedure_type.eq.NE_edge_type).or.
     $                        (procedure_type.eq.NW_edge_type).or.
     $                        (procedure_type.eq.SE_edge_type).or.
     $                        (procedure_type.eq.SW_edge_type)

        end function is_an_anti_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the boundary layer is a corner
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary layer represented as [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        function is_a_corner(bc_section)

          implicit none

          integer, dimension(5), intent(in) :: bc_section
          logical                           :: is_a_corner

          integer :: procedure_type

          procedure_type = bc_section(1)

          is_a_corner = (procedure_type.eq.NE_corner_type).or.
     $                  (procedure_type.eq.NW_corner_type).or.
     $                  (procedure_type.eq.SE_corner_type).or.
     $                  (procedure_type.eq.SW_corner_type)

        end function is_a_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> overlap the bc_section depnding whether it fits inside 
        !> the integration borders
        !
        !> @date
        !> 02_03_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section_type
        !> type of bc_section (N_edge_type,S_edge_type,...)
        !
        !> @param bc_section_borders
        !> borders of the bc_sections ([i_min,i_max]x[j_min,j_max])
        !
        !> @param bc_section_overlap_type
        !> overlap of the bc_section (no_overlap,...)
        !
        !> @param gen_borders
        !> integration borders ([i_min,i_max]x[j_min,j_max])
        !--------------------------------------------------------------
        subroutine overlap_bc_section_by_integration_borders(
     $       bc_section_type,
     $       bc_section_borders,
     $       bc_section_overlap_type,
     $       gen_borders)

          implicit none

          integer                       , intent(in)    :: bc_section_type
          integer(ikind), dimension(2,2), intent(inout) :: bc_section_borders
          integer                       , intent(inout) :: bc_section_overlap_type
          integer(ikind), dimension(2,2), intent(in)    :: gen_borders


          integer                        :: k
          integer(ikind), dimension(2,2) :: overlap_borders


          do k=1,2
             overlap_borders(k,1) = max(bc_section_borders(k,1),gen_borders(k,1))
             overlap_borders(k,2) = min(bc_section_borders(k,2),gen_borders(k,2))
          end do


          select case(bc_section_type)

            case(E_edge_type,W_edge_type)

               if(overlap_borders(1,2).le.bc_section_borders(1,1)) then
                  if(overlap_borders(1,2).eq.bc_section_borders(1,1)) then
                     call overlap_E(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = EW_overlap
                  end if
               end if

               if(overlap_borders(1,1).ge.bc_section_borders(1,2)) then
                  if(overlap_borders(1,1).eq.bc_section_borders(1,2)) then
                     call overlap_W(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = EW_overlap
                  end if
               end if

               bc_section_borders(2,1) = overlap_borders(2,1)
               bc_section_borders(2,2) = overlap_borders(2,2)
               

            case(N_edge_type,S_edge_type)

               if(overlap_borders(2,2).le.bc_section_borders(2,1)) then
                  if(overlap_borders(2,2).eq.bc_section_borders(2,1)) then
                     call overlap_N(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = NS_overlap
                  end if
               end if

               if(overlap_borders(2,1).ge.bc_section_borders(2,2)) then
                  if(overlap_borders(2,1).eq.bc_section_borders(2,2)) then
                     call overlap_S(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = NS_overlap
                  end if
               end if

               bc_section_borders(1,1) = overlap_borders(1,1)
               bc_section_borders(1,2) = overlap_borders(1,2)

            case(   NE_corner_type, NW_corner_type,
     $              SE_corner_type, SW_corner_type,
     $              NE_edge_type, NW_edge_type,
     $              SE_edge_type, SW_edge_type)

               if(overlap_borders(1,2).le.bc_section_borders(1,1)) then
                  if(overlap_borders(1,2).eq.bc_section_borders(1,1)) then
                     call overlap_E(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = EW_overlap
                  end if
               end if

               if(overlap_borders(1,1).ge.bc_section_borders(1,2)) then
                  if(overlap_borders(1,1).eq.bc_section_borders(1,2)) then
                     call overlap_W(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = EW_overlap
                  end if
               end if

               if(overlap_borders(2,1).ge.bc_section_borders(2,2)) then
                  if(overlap_borders(2,1).eq.bc_section_borders(2,2)) then
                     call overlap_S(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = NS_overlap
                  end if
               end if

               if(overlap_borders(2,2).le.bc_section_borders(2,1)) then
                  if(overlap_borders(2,2).eq.bc_section_borders(2,1)) then
                     call overlap_N(bc_section_overlap_type)
                  else
                     bc_section_overlap_type = NS_overlap
                  end if
               end if              
            
            case default
               print '(''bf_bc_sections_overlap_module'')'
               print '(''overlap_bc_section_by_integration_borders'')'
               print '(''bc_section_type not recognized: '',I2)', bc_section_type
               stop ''
          end select

        end subroutine overlap_bc_section_by_integration_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the bc_section to remove
        !> the computation of North the grid-points
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> overlapping of the bc_section
        !--------------------------------------------------------------
        subroutine overlap_N(overlap_type)

          implicit none

          integer, intent(inout) :: overlap_type

          
          select case(overlap_type)

            case(no_overlap)
               overlap_type = N_overlap

            case(N_overlap,NE_overlap,NW_overlap,NS_overlap,EW_overlap)
               overlap_type = overlap_type

            case(S_overlap,SE_overlap,SW_overlap)
               overlap_type = NS_overlap

            case(E_overlap)
               overlap_type = NE_overlap

            case(W_overlap)
               overlap_type = NW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_N',
     $              overlap_type)

          end select

        end subroutine overlap_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the bc_section to remove
        !> the computation of South the grid-points
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> overlapping of the bc_section
        !--------------------------------------------------------------
        subroutine overlap_S(overlap_type)

          implicit none

          integer, intent(inout) :: overlap_type

          
          select case(overlap_type)

            case(no_overlap)
               overlap_type = S_overlap

            case(S_overlap,SE_overlap,SW_overlap,NS_overlap,EW_overlap)
               overlap_type = overlap_type

            case(N_overlap,NE_overlap,NW_overlap)
               overlap_type = NS_overlap

            case(E_overlap)
               overlap_type = SE_overlap

            case(W_overlap)
               overlap_type = SW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_S',
     $              overlap_type)

          end select

        end subroutine overlap_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the bc_section to remove
        !> the computation of East the grid-points
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> overlapping of the bc_section
        !--------------------------------------------------------------
        subroutine overlap_E(overlap_type)

          implicit none

          integer, intent(inout) :: overlap_type

          
          select case(overlap_type)

            case(no_overlap)
               overlap_type = E_overlap

            case(E_overlap,SE_overlap,NE_overlap,NS_overlap,EW_overlap)
               overlap_type = overlap_type
               
            case(W_overlap,NW_overlap,SW_overlap)
               overlap_type = EW_overlap

            case(N_overlap)
               overlap_type = NE_overlap

            case(S_overlap)
               overlap_type = SE_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_E',
     $              overlap_type)

          end select

        end subroutine overlap_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the bc_section to remove
        !> the computation of West the grid-points
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> overlapping of the bc_section
        !--------------------------------------------------------------
        subroutine overlap_W(overlap_type)

          implicit none

          integer, intent(inout) :: overlap_type

          
          select case(overlap_type)

            case(no_overlap)
               overlap_type = W_overlap

            case(W_overlap,SW_overlap,NW_overlap,NS_overlap,EW_overlap)
               overlap_type = overlap_type
               
            case(E_overlap,SE_overlap,NE_overlap)
               overlap_type = EW_overlap

            case(N_overlap)
               overlap_type = NW_overlap

            case(S_overlap)
               overlap_type = SW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_W',
     $              overlap_type)

          end select

        end subroutine overlap_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine which grid-points are computed when
        !> the corner or anti-corner bc_section is overlap
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> index identifying how the anti-corner bc_section
        !> is overlaped by a neighboring corner bc_section
        !
        !> @param compute_point1
        !> logical indicating whether the SW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point2
        !> logical indicating whether the SE grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point3
        !> logical indicating whether the NW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point4
        !> logical indicating whether the NE grid-point of
        !> the bc_section should be computed
        !--------------------------------------------------------------
        subroutine determine_corner_or_anti_corner_grdpts_computed(
     $     cpt_overlap_type,
     $     overlap_type,
     $     compute_point)

          implicit none

          integer              , intent(in)  :: cpt_overlap_type
          integer              , intent(in)  :: overlap_type
          integer, dimension(4), intent(out) :: compute_point


          logical :: compute_pt
          

          !1st point
          !--------------------------------------------------
          compute_pt = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          if(compute_pt) then
             
             select case(cpt_overlap_type)

               case(
     $            cpt1normal_and_cpt4normal,
     $            cpt1normal_and_cpt4not,
     $            cpt1normal_and_cpt4overlap)

                 compute_point(1) = cptnormal_type

               case(
     $                cpt1not_and_cpt4normal,
     $                cpt1not_and_cpt4not,
     $                cpt1not_and_cpt4overlap)

                 compute_point(1) = cptnot_type

               case(
     $                cpt1overlap_and_cpt4normal,
     $                cpt1overlap_and_cpt4not,
     $                cpt1overlap_and_cpt4overlap)

                 compute_point(1) = cptoverlap_type

               case(
     $                cpt2normal_and_cpt3not,
     $                cpt2normal_and_cpt3overlap,
     $                cpt2not_and_cpt3normal,
     $                cpt2not_and_cpt3not,    
     $                cpt2not_and_cpt3overlap,
     $                cpt2overlap_and_cpt3normal,
     $                cpt2overlap_and_cpt3not,    
     $                cpt2overlap_and_cpt3overlap)

                 compute_point(1) = cptnormal_type

               case default

                  call error_cpt_overlap_index(
     $                 'bf_layer_bc_scetions_overlap_module',
     $                 'determine_corner_or_anti_corner_grdpts_computed',
     $                 cpt_overlap_type)
                  
             end select

          else
             compute_point(1) = cptnot_type
          end if


          !2nd point
          !--------------------------------------------------
          compute_pt = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          if(compute_pt) then
             
             select case(cpt_overlap_type)

               case(
     $            cpt2normal_and_cpt3normal,
     $            cpt2normal_and_cpt3not,
     $            cpt2normal_and_cpt3overlap)

                 compute_point(2) = cptnormal_type

               case(
     $                cpt2not_and_cpt3normal,
     $                cpt2not_and_cpt3not,
     $                cpt2not_and_cpt3overlap)

                 compute_point(2) = cptnot_type

               case(
     $                cpt2overlap_and_cpt3normal,
     $                cpt2overlap_and_cpt3not,
     $                cpt2overlap_and_cpt3overlap)

                 compute_point(2) = cptoverlap_type

               case(
     $                cpt1normal_and_cpt4not,
     $                cpt1normal_and_cpt4overlap,
     $                cpt1not_and_cpt4normal,
     $                cpt1not_and_cpt4not,    
     $                cpt1not_and_cpt4overlap,
     $                cpt1overlap_and_cpt4normal, 
     $                cpt1overlap_and_cpt4not,    
     $                cpt1overlap_and_cpt4overlap)

                 compute_point(2) = cptnormal_type

               case default

                  call error_cpt_overlap_index(
     $                 'bf_layer_bc_scetions_overlap_module',
     $                 'determine_corner_or_anti_corner_grdpts_computed',
     $                 cpt_overlap_type)
                  
             end select

          else
             compute_point(2) = cptnot_type
          end if


          !3rd point
          !--------------------------------------------------
          compute_pt = .not.(
     $         (overlap_type.eq.N_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          if(compute_pt) then
             
             select case(cpt_overlap_type)

               case(
     $            cpt2normal_and_cpt3normal,
     $            cpt2not_and_cpt3normal,
     $            cpt2overlap_and_cpt3normal)

                 compute_point(3) = cptnormal_type

               case(
     $                cpt2normal_and_cpt3not,
     $                cpt2not_and_cpt3not,
     $                cpt2overlap_and_cpt3not)

                 compute_point(3) = cptnot_type

               case(
     $                cpt2normal_and_cpt3overlap,
     $                cpt2not_and_cpt3overlap,
     $                cpt2overlap_and_cpt3overlap)

                 compute_point(3) = cptoverlap_type

               case(
     $                cpt1normal_and_cpt4not,
     $                cpt1normal_and_cpt4overlap,
     $                cpt1not_and_cpt4normal,
     $                cpt1not_and_cpt4not,    
     $                cpt1not_and_cpt4overlap,
     $                cpt1overlap_and_cpt4normal, 
     $                cpt1overlap_and_cpt4not,    
     $                cpt1overlap_and_cpt4overlap)

                 compute_point(3) = cptnormal_type

               case default

                  call error_cpt_overlap_index(
     $                 'bf_layer_bc_scetions_overlap_module',
     $                 'determine_corner_or_anti_corner_grdpts_computed',
     $                 cpt_overlap_type)
                  
             end select

          else
             compute_point(3) = cptnot_type
          end if


          !4th point
          !--------------------------------------------------
          compute_pt = .not.(
     $         (overlap_type.eq.N_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          if(compute_pt) then
             
             select case(cpt_overlap_type)

               case(
     $            cpt1normal_and_cpt4normal,
     $            cpt1not_and_cpt4normal,
     $            cpt1overlap_and_cpt4normal)

                 compute_point(4) = cptnormal_type

               case(
     $                cpt1normal_and_cpt4not,
     $                cpt1not_and_cpt4not,
     $                cpt1overlap_and_cpt4not)

                 compute_point(4) = cptnot_type

               case(
     $                cpt1normal_and_cpt4overlap,
     $                cpt1not_and_cpt4overlap,
     $                cpt1overlap_and_cpt4overlap)

                 compute_point(4) = cptoverlap_type

               case(
     $                cpt2normal_and_cpt3not,
     $                cpt2normal_and_cpt3overlap,
     $                cpt2not_and_cpt3normal,
     $                cpt2not_and_cpt3not,    
     $                cpt2not_and_cpt3overlap,
     $                cpt2overlap_and_cpt3normal, 
     $                cpt2overlap_and_cpt3not,    
     $                cpt2overlap_and_cpt3overlap)

                 compute_point(4) = cptnormal_type

               case default

                  call error_cpt_overlap_index(
     $                 'bf_layer_bc_scetions_overlap_module',
     $                 'determine_corner_or_anti_corner_grdpts_computed',
     $                 cpt_overlap_type)
                  
             end select

          else
             compute_point(4) = cptnot_type
          end if

        end subroutine determine_corner_or_anti_corner_grdpts_computed


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine which grid-points are computed when
        !> the corner or anti-corner bc_section is overlap
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> index identifying how the anti-corner bc_section
        !> is overlaped by a neighboring corner bc_section
        !
        !> @param compute_point1
        !> logical indicating whether the SW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point2
        !> logical indicating whether the SE grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point3
        !> logical indicating whether the NW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point4
        !> logical indicating whether the NE grid-point of
        !> the bc_section should be computed
        !--------------------------------------------------------------
        subroutine determine_edge_grdpts_computed(
     $     overlap_type,
     $     compute_edge)

          implicit none

          integer              , intent(in)  :: overlap_type
          logical, dimension(2), intent(out) :: compute_edge

          select case(overlap_type)

            case(no_overlap)
               compute_edge(1) = .true.
               compute_edge(2) = .true.

            case(NS_overlap,EW_overlap)
               compute_edge(1) = .false.
               compute_edge(2) = .false.
               
            case(N_overlap,E_overlap)
               compute_edge(1) = .true.
               compute_edge(2) = .false.

            case(S_overlap,W_overlap)
               compute_edge(1) = .false.
               compute_edge(2) = .true.

            case default
               print '(''bf_layer_bc_sections_overlap_module'')'
               print '(''determine_edge_grdpts_computed'')'
               print '(''overlap type not compatible'')'
               
          end select

        end subroutine determine_edge_grdpts_computed

      end module bf_layer_bc_sections_overlap_module
