      !> @file
      !> subroutines for the identification of the crenels
      !> to the inside of the computational domain that should
      !> be activated to remove the crenels
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines for the identification of the crenels
      !> to the inside of the computational domain that should
      !> be activated to remove the crenels
      !
      !> @date
      ! 28_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_bc_sections_icr_module

        use parameters_bf_layer, only :
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt
      
        use parameters_constant, only :
     $       left, right,
     $       
     $       no_bc_procedure_type,
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
     $       get_edge_crenel_id_param


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the grid-points that need to be tested
        !> to identify whether there is a crenel to the inside
        !> of the computational domain
        !
        !> @date
        !> 28_04_2015 - initial version - J.L. Desmarais
        !
        !@param[in] edge_bc_section
        ! edge bc_section whose overlap with neighboring
        ! anti-corner is investigated
        !
        !@param[in] side
        ! either left or right, indicating which side of the
        ! bc_section is studied
        !
        !@param[out] grdpts_ex_borders
        ! [i_min,i_max]x[j_min,j_max] of the grdpts_id to be extracted
        ! for the tests
        !
        !@param[out] test1_loc_borders
        ! which grdpts from the test1_array are tested
        !
        !@param[out] test1_array
        ! the pattern the grdpts_id should match for the test1
        !
        !@param[out] test2_loc_borders
        ! which grdpts from the test2_array are tested
        !
        !@param[out] test2_array
        ! the pattern the grdpts_id should match for the test2
        !------------------------------------------------------------
        subroutine get_edge_crenel_id_param(
     $       edge_bc_section,
     $       side,
     $       grdpts_ex_borders,
     $       test1_loc_borders,
     $       test1_array,
     $       test2_loc_borders,
     $       test2_array)

          implicit none

          integer(ikind), dimension(5)  , intent(in)  :: edge_bc_section
          logical                       , intent(in)  :: side
          integer(ikind), dimension(2,2), intent(out) :: grdpts_ex_borders
          integer(ikind), dimension(2,2), intent(out) :: test1_loc_borders
          integer(ikind), dimension(3,3), intent(out) :: test1_array
          integer(ikind), dimension(2,2), intent(out) :: test2_loc_borders
          integer(ikind), dimension(3,3), intent(out) :: test2_array
          
          integer(ikind) :: i_min
          integer(ikind) :: i_max
          integer(ikind) :: j_min
          integer(ikind) :: j_max


          select case(edge_bc_section(1))

            case(N_edge_type)

               i_min = edge_bc_section(2)
               j_min = edge_bc_section(3)
               i_max = edge_bc_section(4)

               test1_loc_borders = reshape((/1,1,2,3/),(/2,2/))
               test2_loc_borders = reshape((/1,1,2,3/),(/2,2/))

               if(side.eq.left) then
                  grdpts_ex_borders = reshape((/
     $                 i_min-2, j_min,
     $                 i_min-1, j_min+2/),
     $                 (/2,2/))                  

                  test1_array(1:2,1:3) = reshape((/
     $                 bc_interior_pt, bc_interior_pt,
     $                 bc_interior_pt, bc_pt         ,
     $                 bc_interior_pt, bc_pt         /),
     $                 (/2,3/))

                  test2_array(1:2,1:3) = reshape((/
     $                 bc_interior_pt, bc_interior_pt,
     $                 bc_interior_pt, bc_pt         ,
     $                 bc_pt         , bc_pt         /),
     $                 (/2,3/))

               else
                  grdpts_ex_borders = reshape((/
     $                 i_max+1, j_min,
     $                 i_max+2, j_min+2/),
     $                 (/2,2/))

                  test1_array(1:2,1:3) = reshape((/
     $                 bc_interior_pt, bc_interior_pt,
     $                 bc_pt         , bc_interior_pt,
     $                 bc_pt         , bc_interior_pt/),
     $                 (/2,3/))

                  test2_array(1:2,1:3) = reshape((/
     $                 bc_interior_pt, bc_interior_pt,
     $                 bc_pt         , bc_interior_pt,
     $                 bc_pt         , bc_pt         /),
     $                 (/2,3/))

               end if


            case(S_edge_type)

               i_min = edge_bc_section(2)
               j_min = edge_bc_section(3)
               i_max = edge_bc_section(4)

               test1_loc_borders = reshape((/1,1,2,3/),(/2,2/))
               test2_loc_borders  = reshape((/1,1,2,3/),(/2,2/))

               if(side.eqv.left) then
                  grdpts_ex_borders = reshape((/
     $                 i_min-2,j_min-1,
     $                 i_min-1,j_min+1/),
     $                 (/2,2/))

                  test1_array(1:2,1:3) = reshape((/
     $                 bc_interior_pt, bc_pt,
     $                 bc_interior_pt, bc_pt,
     $                 bc_interior_pt, bc_interior_pt/),
     $                 (/2,3/))

                  test2_array(1:2,1:3) = reshape((/
     $                 bc_pt         , bc_pt,
     $                 bc_interior_pt, bc_pt,
     $                 bc_interior_pt, bc_interior_pt/),
     $                 (/2,3/))

               else
                  grdpts_ex_borders = reshape((/
     $                 i_max+1,j_min-1,
     $                 i_max+2,j_min+1/),
     $                 (/2,2/))

                  test1_array(1:2,1:3) = reshape((/
     $                 bc_pt         , bc_interior_pt,
     $                 bc_pt         , bc_interior_pt,
     $                 bc_interior_pt, bc_interior_pt/),
     $                 (/2,3/))

                  test2_array(1:2,1:3) = reshape((/
     $                 bc_pt         , bc_pt,
     $                 bc_pt         , bc_interior_pt,
     $                 bc_interior_pt, bc_interior_pt/),
     $                 (/2,3/))

               end if
               

            case(E_edge_type)

               i_min = edge_bc_section(2)
               j_min = edge_bc_section(3)
               j_max = edge_bc_section(4)

               test1_loc_borders = reshape((/1,1,3,2/),(/2,2/))
               test2_loc_borders  = reshape((/1,1,3,2/),(/2,2/))

               if(side.eqv.left) then
                  grdpts_ex_borders = reshape((/
     $                 i_min  ,j_min-2,
     $                 i_min+2,j_min-1/),
     $                 (/2,2/))

                  test1_array(1:3,1:2) = reshape((/
     $                 bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $                 bc_interior_pt, bc_pt         , bc_pt         /),
     $                 (/3,2/))

                  test2_array(1:3,1:2) = reshape((/
     $                 bc_interior_pt, bc_interior_pt, bc_pt,
     $                 bc_interior_pt, bc_pt         , bc_pt/),
     $                 (/3,2/))

               else
                  grdpts_ex_borders = reshape((/
     $                 i_min  ,j_max+1,
     $                 i_min+2,j_max+2/),
     $                 (/2,2/))

                  test1_array(1:3,1:2) = reshape((/
     $                 bc_interior_pt, bc_pt         , bc_pt         ,
     $                 bc_interior_pt, bc_interior_pt, bc_interior_pt/),
     $                 (/3,2/))

                  test2_array(1:3,1:2) = reshape((/
     $                 bc_interior_pt, bc_pt         , bc_pt,
     $                 bc_interior_pt, bc_interior_pt, bc_pt/),
     $                 (/3,2/))

               end if
               

            case(W_edge_type)

               i_min = edge_bc_section(2)
               j_min = edge_bc_section(3)
               j_max = edge_bc_section(4)

               test1_loc_borders = reshape((/1,1,3,2/),(/2,2/))
               test2_loc_borders = reshape((/1,1,3,2/),(/2,2/))

               if(side.eqv.left) then
                  grdpts_ex_borders = reshape((/
     $                 i_min-1,j_min-2,
     $                 i_min+1,j_min-1/),
     $                 (/2,2/))

                  test1_array(1:3,1:2) = reshape((/
     $                 bc_interior_pt, bc_interior_pt, bc_interior_pt, 
     $                 bc_pt         , bc_pt         , bc_interior_pt/),
     $                 (/3,2/))

                  test2_array(1:3,1:2) = reshape((/
     $                 bc_pt, bc_interior_pt, bc_interior_pt, 
     $                 bc_pt, bc_pt         , bc_interior_pt/),
     $                 (/3,2/))

               else
                  grdpts_ex_borders = reshape((/
     $                 i_min-1,j_max+1,
     $                 i_min+1,j_max+2/),
     $                 (/2,2/))

                  test1_array(1:3,1:2) = reshape((/
     $                 bc_pt         , bc_pt         , bc_interior_pt,
     $                 bc_interior_pt, bc_interior_pt, bc_interior_pt/),
     $                 (/3,2/))

                  test2_array(1:3,1:2) = reshape((/
     $                 bc_pt, bc_pt         , bc_interior_pt,
     $                 bc_pt, bc_interior_pt, bc_interior_pt/), 
     $                 (/3,2/))

               end if

            case default
               print '(''bf_layer_bc_sections_icr_module'')'
               print '(''get_edge_crenel_id_param'')'
               print '(''edge type not recognized'')'
               print '(''edge_type: '',I3)', edge_bc_section(1)

          end select

        end subroutine get_edge_crenel_id_param


      end module bf_layer_bc_sections_icr_module
