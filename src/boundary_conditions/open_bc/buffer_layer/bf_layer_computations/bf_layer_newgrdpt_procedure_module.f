      !> @file
      !> identification of the procedure needed to compute the new
      !> grid points
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module with the subroutines identifying the procedure
      !> to compute the new grid points
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_layer_newgrdpt_procedure_module

        use bf_layer_bc_procedure_module, only : 
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use bf_layer_errors_module, only : 
     $       error_mainlayer_id

        use bf_layer_sync_module, only :
     $       get_sync_indices_to_extract_interior_data

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: 
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type,
     $       get_newgrdpt_procedure,
     $       error_gradient_type,
     $       get_grdpts_id_from_interior,
     $       get_interior_data_for_newgrdpt,
     $       are_intermediate_newgrdpt_data_needed,
     $       get_x_map_for_newgrdpt,
     $       get_y_map_for_newgrdpt


        integer, parameter :: no_gradient_type=0
        integer, parameter :: gradient_I_type=1
        integer, parameter :: gradient_L0_type=2
        integer, parameter :: gradient_R0_type=3
        integer, parameter :: gradient_xLR0_yI_type=4
        integer, parameter :: gradient_xI_yLR0_type=5
        integer, parameter :: gradient_xLR0_yLR0_type=6


        contains

        !i: x-index identifying the position of the no_pt in grdpts_id
        !j: y-index identifying the position of the no_pt in grdpts_id
        !procedure_type : type of procedure needed to compute the new grdpt
        !gradient_type  : number of grdpts in the surroundings to compute the new grdpt
        subroutine get_newgrdpt_procedure(
     $       i,j,
     $       grdpts_id,
     $       procedure_type,
     $       gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type


          !  -------
          ! |       |
          ! |   0*  |
          ! | 3     |
          !  -------
          !============================================================
          !procedure 1
          !============================================================
          if(grdpts_id(i-1,j-1).eq.bc_pt) then


          !  -------
          ! |       |
          ! |   0*  |
          ! | 3 3   |
          !  -------
             if(grdpts_id(i,j-1).eq.bc_pt) then

                call get_newgrdpt_procedure_1_1(
     $               i,j,grdpts_id,
     $               procedure_type,gradient_type)

             else

          !  -------
          ! |       |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
                if(grdpts_id(i-1,j).eq.bc_pt) then

                   call get_newgrdpt_procedure_1_2(
     $                  i,j,grdpts_id,
     $                  procedure_type,gradient_type)

          !  -------
          ! |       |
          ! | x 0*  |
          ! | 3 x   |
          !  -------
                else

                   procedure_type = NE_corner_type
                   gradient_type  = no_gradient_type

                end if

             end if

          else

          !  -------
          ! |       |
          ! |   0*  |
          ! | x 3   |
          !  -------
          !============================================================
          !procedure 2
          !============================================================
             if(grdpts_id(i,j-1).eq.bc_pt) then

          !  -------
          ! |       |
          ! |   0*  |
          ! | x 3 3 |
          !  -------
                if(grdpts_id(i+1,j-1).eq.bc_pt) then

                   call get_newgrdpt_procedure_2_1(
     $                  i,j,grdpts_id,
     $                  procedure_type,gradient_type)

                else

                   call error_newgrdpt_procedure(i,j,grdpts_id)

                end if

             else

          !  -------
          ! |       |
          ! |   0*  |
          ! | x x 3 |
          !  -------
          !============================================================
          !procedure 3
          !============================================================
                if(grdpts_id(i+1,j-1).eq.bc_pt) then

          !  -------
          ! |       |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
                   if(grdpts_id(i+1,j).eq.bc_pt) then

                      call get_newgrdpt_procedure_3_1(
     $                     i,j,grdpts_id,
     $                     procedure_type,gradient_type)

          !  -------
          ! |       |
          ! |   0*x |
          ! | x x 3 |
          !  -------
                   else

                      procedure_type = NW_corner_type
                      gradient_type  = no_gradient_type
                      

                   end if

                else

          !  -------
          ! |       |
          ! | 3 0*  |
          ! | x x x |
          !  -------
          !============================================================
          !procedure 4
          !============================================================
                   if(grdpts_id(i-1,j).eq.bc_pt) then
                      
          !  -------
          ! | - - 3 |
          ! | 3 0*  |
          ! | x x x |
          !  -------
          ! WARNING: assumptions:
          !   - grdpts_id(i-1,j+1) = bc_pt
          !   - grdpts_id(i  ,j+1) = bc_pt
          !-------------------------------
                      if(grdpts_id(i+1,j+1).eq.bc_pt) then

                         procedure_type = SE_edge_type
                         gradient_type  = gradient_xI_yLR0_type

                      else

          !  -------
          ! | - 3 x |
          ! | 3 0*  |
          ! | x x x |
          !  -------
          ! WARNING: assumptions:
          !   - grdpts_id(i-1,j+1) = bc_pt
          !-------------------------------
                         if(grdpts_id(i,j+1).eq.bc_pt) then

                            procedure_type = SE_edge_type
                            gradient_type  = gradient_xLR0_yLR0_type

                         else

          !  -------
          ! | 3 x x |
          ! | 3 0*  |
          ! | x x x |
          !  -------
                            if(grdpts_id(i-1,j+1).eq.bc_pt) then
                               
                               procedure_type = E_edge_type
                               gradient_type  = gradient_L0_type

          !  -------
          ! | x x x |
          ! | 3 0*  |
          ! | x x x |
          !  -------
                            else

                               call error_newgrdpt_procedure(i,j,grdpts_id)

                            end if

                         end if

                      end if

                   else

          !  -------
          ! |       |
          ! | x 0*3 |
          ! | x x x |
          !  -------
          !============================================================
          !procedure 5
          !============================================================
                      if(grdpts_id(i+1,j).eq.bc_pt) then

          !  -------
          ! | 3 - - |
          ! | x 0*3 |
          ! | x x x |
          !  -------
          ! WARNING: assumptions:
          !   - grdpts_id(i,j+1)   = bc_pt
          !   - grdpts_id(i+1,j+1) = bc_pt
          !-------------------------------
                         if(grdpts_id(i-1,j+1).eq.bc_pt) then

                            procedure_type = SW_edge_type
                            gradient_type  = gradient_xI_yLR0_type


                         else

          !  -------
          ! | x 3 - |
          ! | x 0*3 |
          ! | x x x |
          !  -------
          ! WARNING: assumptions:
          !   - grdpts_id(i+1,j+1) = bc_pt
          !-------------------------------
                            if(grdpts_id(i,j+1).eq.bc_pt) then

                               procedure_type = SW_edge_type
                               gradient_type  = gradient_xLR0_yLR0_type

                            else
          !  -------
          ! | x x 3 |
          ! | x 0*3 |
          ! | x x x |
          !  -------
                               if(grdpts_id(i+1,j+1).eq.bc_pt) then

                                  procedure_type = W_edge_type
                                  gradient_type  = gradient_L0_type
          !  -------
          ! | x x x |
          ! | x 0*3 |
          ! | x x x |
          !  -------
                               else
                                  
                                  call error_newgrdpt_procedure(i,j,grdpts_id)

                               end if
                            end if
                         end if                               


                      else
          !  -------
          ! | 3     |
          ! | x 0*x |
          ! | x x x |
          !  -------
          !============================================================
          !procedure 6
          !============================================================
                         if(grdpts_id(i-1,j+1).eq.bc_pt) then

          !  -------
          ! | 3 - 3 |
          ! | x 0*x |
          ! | x x x |
          !  -------
          ! WARNING: assumptions:
          !   - grdpts_id(i,j+1)   = bc_pt
          !-------------------------------
                            if(grdpts_id(i+1,j+1).eq.bc_pt) then

                               procedure_type = S_edge_type
                               gradient_type  = gradient_I_type

                            else
          !  -------
          ! | 3 3 x |
          ! | x 0*x |
          ! | x x x |
          !  -------
                               if(grdpts_id(i,j+1).eq.bc_pt) then
                                  
                                  procedure_type = S_edge_type
                                  gradient_type  = gradient_R0_type
          !  -------
          ! | 3 x   |
          ! | x 0*x |
          ! | x x x |
          !  -------
                               else
                                  procedure_type = SE_corner_type
                                  gradient_type  = no_gradient_type
                               end if
                            end if
                         else
          !  -------
          ! | x 3   |
          ! | x 0*x |
          ! | x x x |
          !  -------
          !============================================================
          !procedure 7
          !============================================================
                            if(grdpts_id(i,j+1).eq.bc_pt) then

          !  -------
          ! | x 3 3 |
          ! | x 0*x |
          ! | x x x |
          !  -------
                               if(grdpts_id(i+1,j+1).eq.bc_pt) then

                                  procedure_type = S_edge_type
                                  gradient_type  = gradient_L0_type
          !  -------
          ! | x 3 x |
          ! | x 0*x |
          ! | x x x |
          !  -------
                               else

                                  call error_newgrdpt_procedure(i,j,grdpts_id)

                               end if

          !  -------
          ! | x x   |
          ! | x 0*x |
          ! | x x x |
          !  -------                              
                            else

          !  -------
          ! | x x 3 |
          ! | x 0*x |
          ! | x x x |
          !  ------- 
                               if(grdpts_id(i+1,j+1).eq.bc_pt) then

                                  procedure_type = SW_corner_type
                                  gradient_type  = no_gradient_type

                               else

                                  call error_newgrdpt_procedure(i,j,grdpts_id)

                               end if

                            end if

                         end if

                      end if

                   end if

                end if

             end if

          end if

        end subroutine get_newgrdpt_procedure


        subroutine get_newgrdpt_procedure_1_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type

          
          !  -------
          ! |       |
          ! |   0*  |
          ! | 3 3 3 |
          !  -------
          !============================================================
          !procedure 1.1.1
          !============================================================
          if(grdpts_id(i+1,j-1).eq.bc_pt) then

          !  -------
          ! |     3 |
          ! |   0*- |
          ! | 3 3 3 |
          !  -------
          !WARNING: assumptions:
          ! - grdpts_id(i+1,j) = bc_pt
             if(grdpts_id(i+1,j+1).eq.bc_pt) then

                procedure_type = NW_edge_type
                gradient_type  = gradient_I_type

             else

          !  -------
          ! |     x |
          ! |   0*3 |
          ! | 3 3 3 |
          !  -------
                if(grdpts_id(i+1,j).eq.bc_pt) then

                   procedure_type = NW_edge_type
                   gradient_type  = gradient_xI_yLR0_type

                else

          !  -------
          ! | 3   x |
          ! | - 0*x |
          ! | 3 3 3 |
          !  -------
          !WARNING: assumptions:
          ! - grdpts_id(i-1,j) = bc_pt
                   if(grdpts_id(i-1,j+1).eq.bc_pt) then

                      procedure_type = NE_edge_type
                      gradient_type  = gradient_I_type

                   else
          !  -------
          ! | x   x |
          ! | 3 0*x |
          ! | 3 3 3 |
          !  -------
                      if(grdpts_id(i-1,j).eq.bc_pt) then
                         
                         procedure_type = NE_edge_type
                         gradient_type  = gradient_xI_yLR0_type

                      else
          !  -------
          ! | x   x |
          ! | x 0*x |
          ! | 3 3 3 |
          !  -------
                         procedure_type = N_edge_type
                         gradient_type  = gradient_I_type

                      end if
                   end if
                end if
             end if

          !  -------
          ! |       |
          ! |   0*  |
          ! | 3 3 x |
          !  -------
          !============================================================
          !procedure 1.1.2
          !============================================================
          else

          !  -------
          ! | 3     |
          ! | - 0*  |
          ! | 3 3 x |
          !  -------
             if(grdpts_id(i-1,j+1).eq.bc_pt) then
                
                procedure_type = NE_edge_type
                gradient_type  = gradient_xLR0_yI_type

             else
          !  -------
          ! | x     |
          ! | 3 0*  |
          ! | 3 3 x |
          !  -------
                if(grdpts_id(i-1,j).eq.bc_pt) then

                   procedure_type = NE_edge_type
                   gradient_type  = gradient_xLR0_yLR0_type

          !  -------
          ! | x     |
          ! | x 0*  |
          ! | 3 3 x |
          !  -------
                else
                   
                   procedure_type = N_edge_type
                   gradient_type  = gradient_R0_type

                end if
             end if
          end if          

        end subroutine get_newgrdpt_procedure_1_1


        subroutine get_newgrdpt_procedure_1_2(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type


          !  -------
          ! | - - 3 |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
          !WARNING: assumptions:
          ! - grdpts_id(i-1,j) = bc_pt
          ! - grdpts_id(i,j)   = bc_pt
          if(grdpts_id(i+1,j+1).eq.bc_pt) then

             procedure_type = SE_edge_type
             gradient_type  = gradient_I_type

          else

          !  -------
          ! | - 3 x |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
             if(grdpts_id(i,j+1).eq.bc_pt) then

                procedure_type = SE_edge_type
                gradient_type  = gradient_xLR0_yI_type

             else

          !  -------
          ! | 3 x x |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
                if(grdpts_id(i-1,j+1).eq.bc_pt) then

                   procedure_type = E_edge_type
                   gradient_type  = gradient_I_type

          !  -------
          ! | x x x |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
                else

                   procedure_type = E_edge_type
                   gradient_type  = gradient_R0_type

                end if
             end if
          end if

        end subroutine get_newgrdpt_procedure_1_2


        subroutine get_newgrdpt_procedure_2_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type


          !  -------
          ! |     3 |
          ! |   0*- |
          ! | x 3 3 |
          !  -------
          !WARNING : assumptions
          !  - grdpts(i+1,j) = bc_pt
          if(grdpts_id(i+1,j+1).eq.bc_pt) then

             procedure_type = NW_edge_type
             gradient_type  = gradient_xLR0_yI_type

          else
          !  -------
          ! |     x |
          ! |   0*3 |
          ! | x 3 3 |
          !  -------
             if(grdpts_id(i+1,j).eq.bc_pt) then

                procedure_type = NW_edge_type
                gradient_type  = gradient_xLR0_yLR0_type

             else
          !  -------
          ! |     x |
          ! |   0*x |
          ! | x 3 3 |
          !  -------
                procedure_type = N_edge_type
                gradient_type  = gradient_L0_type

             end if

          end if

        end subroutine get_newgrdpt_procedure_2_1


        subroutine get_newgrdpt_procedure_3_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type


          !  -------
          ! | 3 - - |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
          !WARNING : assumptions
          ! - grdpts_id(i  ,j+1) = bc_pt
          ! - grdpts_id(i+1,j+1) = bc_pt
          if(grdpts_id(i-1,j+1).eq.bc_pt) then

             procedure_type = SW_edge_type
             gradient_type  = gradient_I_type

          else
          !  -------
          ! | x 3 - |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
          !WARNING : assumptions
          ! - grdpts_id(i+1,j+1) = bc_pt
             if(grdpts_id(i,j+1).eq.bc_pt) then

                procedure_type = SW_edge_type
                gradient_type  = gradient_xLR0_yI_type

             else
          !  -------
          ! | x x 3 |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
                if(grdpts_id(i+1,j+1).eq.bc_pt) then
                   
                   procedure_type = W_edge_type
                   gradient_type  = gradient_I_type

          !  -------
          ! | x x x |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
                else

                   procedure_type = W_edge_type
                   gradient_type  = gradient_R0_type

                end if
             end if
          end if

        end subroutine get_newgrdpt_procedure_3_1


        subroutine error_newgrdpt_procedure(i,j,grdpts_id)

          implicit none

          integer(ikind)         , intent(in) :: i
          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(in) :: grdpts_id

          integer :: k

          do k=j-1,j+1
             print '(3I2)', grdpts_id(i-1:i+1,(j+1)-k+1)
          end do

        end subroutine error_newgrdpt_procedure


        subroutine error_gradient_type(
     $     module_name,
     $     subroutine_name,
     $     gradient_type)

          implicit none

          character(*), intent(in) :: module_name
          character(*), intent(in) :: subroutine_name
          integer     , intent(in) :: gradient_type

          print *, module_name
          print *, subroutine_name
          print '(''gradient_type not recognized: '',I2)', gradient_type
          stop ''

        end subroutine error_gradient_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id, the coordinate maps and the
        !> nodes at t-dt and t corresponding to the general
        !> coordinates gen_coords
        !    ___________________
        !   |                  _|_________
        !   |    interior     |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_nodes0
        !> interior nodes at t=t-dt
        !
        !>@param interior_nodes1
        !> interior nodes at t=t
        !
        !>@param tmp_grdpts_id0
        !> array with the grdpts_id data
        !
        !>@param tmp_nodes0
        !> array with the grid points data at t-dt
        !
        !>@param tmp_nodes1
        !> array with the grid points data at t
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> tmp arrays computed
        !--------------------------------------------------------------
        subroutine get_interior_data_for_newgrdpt(
     $     interior_nodes0,
     $     interior_nodes1,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_coords)

          implicit none

          real(rkind)    , dimension(nx,ny,ne)                          , intent(in)    :: interior_nodes0
          real(rkind)    , dimension(nx,ny,ne)                          , intent(in)    :: interior_nodes1
          integer        , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1)   , intent(inout) :: tmp_grdpts_id0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes1
          integer(ikind) , dimension(2,2)                               , intent(in)    :: gen_coords


          integer(ikind)                 :: size_x,size_y
          integer(ikind)                 :: i_recv,i_send,j_recv,j_send
          integer(ikind)                 :: i,j
          integer                        :: k


          !synchronize the nodes
          call get_sync_indices_to_extract_interior_data(
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)


          do k=1,ne
             do j=1, size_y
                do i=1, size_x

                   tmp_nodes0(i_recv+i-1,j_recv+j-1,k) =
     $                  interior_nodes0(i_send+i-1,j_send+j-1,k)

                   tmp_nodes1(i_recv+i-1,j_recv+j-1,k) =
     $                  interior_nodes1(i_send+i-1,j_send+j-1,k)

                end do
             end do
          end do


          !synchronize the grdpts_id
          !initialize with no_pt
          call get_grdpts_id_from_interior(
     $         tmp_grdpts_id0,
     $         gen_coords)          

        end subroutine get_interior_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id corresponding to the general
        !> coordinates gen_coords for the interior
        !    ___________________
        !   |                  _|_________
        !   |    interior     |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> tmp arrays computed
        !--------------------------------------------------------------
        subroutine get_grdpts_id_from_interior(
     $     tmp_grdpts_id,
     $     gen_coords)
        
          implicit none
          
          integer       , dimension(:,:), intent(out) :: tmp_grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: gen_coords

          integer                        :: i,j
          integer(ikind), dimension(2,2) :: grdpts_id_coords

          do j=1, size(tmp_grdpts_id,2)
             do i=1, size(tmp_grdpts_id,1)
                tmp_grdpts_id(i,j) = no_pt
             end do
          end do

          !overlap with the bc_pt if any
          grdpts_id_coords(1,1) = 1
          grdpts_id_coords(1,2) = nx
          grdpts_id_coords(2,1) = 1
          grdpts_id_coords(2,2) = ny

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         bc_pt)

          !overlap with the bc_interior_pt if any
          grdpts_id_coords(1,1) = 2
          grdpts_id_coords(1,2) = nx-1
          grdpts_id_coords(2,1) = 2
          grdpts_id_coords(2,2) = ny-1

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         bc_interior_pt)

          !overlap with the interior_pt if any
          grdpts_id_coords(1,1) = bc_size+1
          grdpts_id_coords(1,2) = nx-bc_size
          grdpts_id_coords(2,1) = bc_size+1
          grdpts_id_coords(2,2) = ny-bc_size

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         interior_pt)

        end subroutine get_grdpts_id_from_interior

      
        subroutine set_grdptid(
     $     grdpts_id,
     $     grdpts_id_coords,
     $     gen_coords,
     $     pt_type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)    :: grdpts_id_coords
          integer(ikind), dimension(2,2), intent(in)    :: gen_coords
          integer                       , intent(in)    :: pt_type


          integer(ikind) :: i_min,i_max,j_min,j_max
          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,j_recv
          integer(ikind) :: i,j


          i_min = max(grdpts_id_coords(1,1), gen_coords(1,1))
          i_max = min(grdpts_id_coords(1,2), gen_coords(1,2))
          j_min = max(grdpts_id_coords(2,1), gen_coords(2,1))
          j_max = min(grdpts_id_coords(2,2), gen_coords(2,2))

          size_x = i_max-i_min+1
          size_y = j_max-j_min+1 

          i_recv = i_min-gen_coords(1,1)+1
          j_recv = j_min-gen_coords(2,1)+1

          
          do j=1,size_y
             do i=1,size_x
                grdpts_id(i_recv+i-1,j_recv+j-1) = pt_type
             end do
          end do

        end subroutine set_grdptid


        function are_intermediate_newgrdpt_data_needed(
     $     mainlayer_id, gen_coords)
     $     result(needed)

          implicit none

          integer                       , intent(in) :: mainlayer_id
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          logical                                    :: needed
          

          select case(mainlayer_id)
            case(N)
               needed = .not.(gen_coords(2,1).gt.ny)

            case(S)
               needed = .not.(gen_coords(2,2).lt.1)
               
            case(E)
               needed = .not.(
     $              (gen_coords(1,1).gt.nx).and.
     $              (gen_coords(2,1).gt.bc_size).and.
     $              (gen_coords(2,2).lt.(ny-bc_size+1)))
               
            case(W)
               needed = .not.(
     $              (gen_coords(1,2).lt.1).and.
     $              (gen_coords(2,1).gt.bc_size).and.
     $              (gen_coords(2,2).lt.(ny-bc_size+1)))
               
            case default
               call error_mainlayer_id(
     $              'bf_layer_newgrdpt_procedure_module',
     $              'are_intermediate_newgrdpt_data_needed',
     $              mainlayer_id)

          end select

        end function are_intermediate_newgrdpt_data_needed


        function get_x_map_for_newgrdpt(
     $     interior_x_map, gen_coords)
     $     result(tmp_x_map)

          implicit none

          real(rkind)   , dimension(nx) , intent(in) :: interior_x_map
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_x_map

          real(rkind) :: dx
          integer     :: i
          integer     :: size

          ! --[-----]---|-------|--------
          if(gen_coords(1,2).le.0) then

             dx = interior_x_map(2)-interior_x_map(1)

             do i=gen_coords(1,1),gen_coords(1,2)
                
                tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(1) + (i-1)*dx

             end do

          else
          ! --------[--|--]-----|--------
             if(gen_coords(1,1).le.0) then

                dx = interior_x_map(2)-interior_x_map(1)

                do i=gen_coords(1,1),0
                   tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(1) + (i-1)*dx
                end do

                size = -gen_coords(1,1) + 1

                do i=1,gen_coords(1,2)
                   tmp_x_map(size+i) = interior_x_map(i)
                end do

          ! -----------|-[-----]|--------
             else
                if(gen_coords(1,2).le.nx) then

                   do i=gen_coords(1,1),gen_coords(1,2)
                      tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(i)
                   end do

          ! -----------|-----[--|--]-----
                else

                   dx = interior_x_map(nx) - interior_x_map(nx-1)

                   if(gen_coords(1,1).le.nx) then
                      
                      do i=gen_coords(1,1),nx
                         tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(i)
                      end do

                      size = -gen_coords(1,1) + 1

                      do i=nx+1,gen_coords(1,2)
                         tmp_x_map(size+i) = interior_x_map(nx) + (i-nx)*dx
                      end do

           ! -----------|--------|-[-----]
                   else

                      do i=gen_coords(1,1),gen_coords(1,2)
                         tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(nx) + (i-nx)*dx
                      end do

                   end if
                end if
             end if
          end if

        end function get_x_map_for_newgrdpt


        function get_y_map_for_newgrdpt(
     $     interior_y_map, gen_coords)
     $     result(tmp_y_map)

          implicit none

          real(rkind)   , dimension(ny) , intent(in) :: interior_y_map
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_y_map

          real(rkind) :: dy
          integer     :: i
          integer     :: size

          ! --[-----]---|-------|--------
          if(gen_coords(2,2).le.0) then

             dy = interior_y_map(2)-interior_y_map(1)

             do i=gen_coords(2,1),gen_coords(2,2)
                
                tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(1) + (i-1)*dy

             end do

          else
          ! --------[--|--]-----|--------
             if(gen_coords(2,1).le.0) then

                dy = interior_y_map(2)-interior_y_map(1)

                do i=gen_coords(2,1),0
                   tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(1) + (i-1)*dy
                end do

                size = -gen_coords(2,1) + 1

                do i=1,gen_coords(2,2)
                   tmp_y_map(size+i) = interior_y_map(i)
                end do

          ! -----------|-[-----]|--------
             else
                if(gen_coords(2,2).le.ny) then

                   do i=gen_coords(2,1),gen_coords(2,2)
                      tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(i)
                   end do

          ! -----------|-----[--|--]-----
                else

                   dy = interior_y_map(ny) - interior_y_map(ny-1)

                   if(gen_coords(2,1).le.ny) then
                      
                      do i=gen_coords(2,1),ny
                         tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(i)
                      end do

                      size = -gen_coords(2,1) + 1

                      do i=ny+1,gen_coords(2,2)
                         tmp_y_map(size+i) = interior_y_map(ny) + (i-ny)*dy
                      end do

           ! -----------|--------|-[-----]
                   else

                      do i=gen_coords(2,1),gen_coords(2,2)
                         tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(ny) + (i-ny)*dy
                      end do

                   end if
                end if
             end if
          end if

        end function get_y_map_for_newgrdpt

      end module bf_layer_newgrdpt_procedure_module
