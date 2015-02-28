      module bf_bc_sections_overlap_module

        use bf_layer_errors_module, only :
     $       error_overlap_index

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
     $       SW_edge_type

        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public ::
     $       overlap_bc_section_by_integration_borders,
     $       overlap_N,
     $       overlap_S,
     $       overlap_E,
     $       overlap_W,
     $       determine_corner_or_anti_corner_grdpts_computed


        contains


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
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
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
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
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
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
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
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
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
     $     overlap_type,
     $     compute_point1,
     $     compute_point2,
     $     compute_point3,
     $     compute_point4)

          implicit none

          integer, intent(in)  :: overlap_type
          logical, intent(out) :: compute_point1
          logical, intent(out) :: compute_point2
          logical, intent(out) :: compute_point3
          logical, intent(out) :: compute_point4


          compute_point1 = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          compute_point2 = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          compute_point3 = .not.(
     $         (overlap_type.eq.N_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.SW_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

          compute_point4 = .not.(
     $         (overlap_type.eq.N_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.NE_overlap).or.
     $         (overlap_type.eq.NW_overlap).or.
     $         (overlap_type.eq.SE_overlap).or.
     $         (overlap_type.eq.NS_overlap).or.
     $         (overlap_type.eq.EW_overlap))

        end subroutine determine_corner_or_anti_corner_grdpts_computed

      end module bf_bc_sections_overlap_module
