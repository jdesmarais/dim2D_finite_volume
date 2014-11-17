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

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt

        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public :: 
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       get_newgrdpt_procedure,
     $       error_gradient_type


        integer, parameter :: no_gradient_type=0
        integer, parameter :: gradient_I_type=1
        integer, parameter :: gradient_L0_type=2
        integer, parameter :: gradient_R0_type=3


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
          if((i.gt.1).and.(j.gt.1).and.(grdpts_id(i-1,j-1).eq.bc_pt)) then


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
             if((j.gt.1).and.(grdpts_id(i,j-1).eq.bc_pt)) then

          !  -------
          ! |       |
          ! |   0*  |
          ! | x 3 3 |
          !  -------
                if((i.lt.size(grdpts_id,1)).and.(grdpts_id(i+1,j-1).eq.bc_pt)) then

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
                if((i.lt.size(grdpts_id,1)).and.(j.gt.1).and.(grdpts_id(i+1,j-1).eq.bc_pt)) then

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
                   if((i.gt.1).and.(grdpts_id(i-1,j).eq.bc_pt)) then
                      
          !  -------
          ! | 3     |
          ! | 3 0*  |
          ! | x x x |
          !  -------
                      if((j.le.size(grdpts_id,2)).and.(grdpts_id(i-1,j+1).eq.bc_pt)) then

                         call get_newgrdpt_procedure_4_1(
     $                     i,j,grdpts_id,
     $                     procedure_type,gradient_type)

                      else

                         call error_newgrdpt_procedure(i,j,grdpts_id)

                      end if                      

                   else

          !  -------
          ! |       |
          ! | x 0*3 |
          ! | x x x |
          !  -------
                      if((i.lt.size(grdpts_id,1)).and.(grdpts_id(i+1,j).eq.bc_pt)) then

                         call get_newgrdpt_procedure_5_1(
     $                        i,j,grdpts_id,
     $                        procedure_type,gradient_type)

                      else
          !  -------
          ! | 3     |
          ! | x 0*x |
          ! | x x x |
          !  -------
                         if((i.gt.1).and.(j.lt.size(grdpts_id,2)).and.(grdpts_id(i-1,j+1).eq.bc_pt)) then

          !  -------
          ! | 3 3   |
          ! | x 0*x |
          ! | x x x |
          !  -------
                            if(grdpts_id(i,j+1).eq.bc_pt) then

                               call get_newgrdpt_procedure_6_1(
     $                              i,j,grdpts_id,
     $                              procedure_type,gradient_type)

          !  -------
          ! | 3 x   |
          ! | x 0*x |
          ! | x x x |
          !  -------
                            else
                               
                               procedure_type = SE_corner_type
                               gradient_type  = no_gradient_type

                            end if


                         else

          !  -------
          ! | x 3   |
          ! | x 0*x |
          ! | x x x |
          !  -------
                            if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i,j+1).eq.bc_pt)) then

          !  -------
          ! | x 3 3 |
          ! | x 0*x |
          ! | x x x |
          !  -------
                               if((i.lt.size(grdpts_id,1)).and.(grdpts_id(i+1,j+1).eq.bc_pt)) then

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
                               if((i.lt.size(grdpts_id,1)).and.(j.lt.size(grdpts_id,2)).and.(grdpts_id(i+1,j+1).eq.bc_pt)) then

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
          if((i.lt.size(grdpts_id,1)).and.(grdpts_id(i+1,j-1).eq.bc_pt)) then

             procedure_type = N_edge_type
             gradient_type  = gradient_I_type


          else

          !  -------
          ! |       |
          ! | 3 0*  |
          ! | 3 3 x |
          !  -------
             if(grdpts_id(i-1,j).eq.bc_pt) then
                
                procedure_type = NE_edge_type
                gradient_type  = no_gradient_type

          !  -------
          ! |       |
          ! | x 0*  |
          ! | 3 3 x |
          !  -------
             else
                
                procedure_type = N_edge_type
                gradient_type  = gradient_R0_type

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
          ! | 3     |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
          if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i-1,j+1).eq.bc_pt)) then

             procedure_type = E_edge_type
             gradient_type  = gradient_I_type

          !  -------
          ! | x     |
          ! | 3 0*  |
          ! | 3 x   |
          !  -------
          else

             procedure_type = E_edge_type
             gradient_type  = gradient_R0_type

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
          ! |       |
          ! |   0*3 |
          ! | x 3 3 |
          !  -------
          if(grdpts_id(i+1,j).eq.bc_pt) then

             procedure_type = NW_edge_type
             gradient_type  = no_gradient_type

          !  -------
          ! |       |
          ! |   0*x |
          ! | x 3 3 |
          !  -------
          else

             procedure_type = N_edge_type
             gradient_type  = gradient_L0_type

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
          ! |   3   |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
          if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i,j+1).eq.bc_pt)) then

             procedure_type = SW_edge_type
             gradient_type  = no_gradient_type

          else
          !  -------
          ! |   x 3 |
          ! |   0*3 |
          ! | x x 3 |
          !  -------
             if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i+1,j+1).eq.bc_pt)) then

                procedure_type = W_edge_type
                gradient_type  = gradient_I_type

          !  -------
          ! |   x x |
          ! |   0*3 |
          ! | x x 3 |
          !  -------                
             else

                procedure_type = W_edge_type
                gradient_type  = gradient_R0_type

             end if

          end if

        end subroutine get_newgrdpt_procedure_3_1


        subroutine get_newgrdpt_procedure_4_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type

          
          !  -------
          ! | 3 3   |
          ! | 3 0*  |
          ! | x x x |
          !  -------
          if(grdpts_id(i,j+1).eq.bc_pt) then

             procedure_type = SE_edge_type
             gradient_type  = no_gradient_type

          !  -------
          ! | 3 x   |
          ! | 3 0*  |
          ! | x x x |
          !  -------
          else

             procedure_type = E_edge_type
             gradient_type  = gradient_L0_type

          end if

        end subroutine get_newgrdpt_procedure_4_1


        subroutine get_newgrdpt_procedure_5_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type

          
          !  -------
          ! |   3   |
          ! | x 0*3 |
          ! | x x x |
          !  -------
          if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i,j+1).eq.bc_pt)) then

             procedure_type = SW_edge_type
             gradient_type  = no_gradient_type

          else

          !  -------
          ! |   x 3 |
          ! | x 0*3 |
          ! | x x x |
          !  -------
             if((j.lt.size(grdpts_id,2)).and.(grdpts_id(i+1,j+1).eq.bc_pt)) then

                procedure_type = W_edge_type
                gradient_type  = gradient_L0_type

          !  -------
          ! |   x x |
          ! | x 0*3 |
          ! | x x x |
          !  -------
             else

                call error_newgrdpt_procedure(i,j,grdpts_id)

             end if

          end if

        end subroutine get_newgrdpt_procedure_5_1


        subroutine get_newgrdpt_procedure_6_1(
     $     i,j,grdpts_id,
     $     procedure_type,gradient_type)

          implicit none

          integer(ikind)         , intent(in)  :: i
          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: gradient_type

          
          !  -------
          ! | 3 3 3 |
          ! | x 0*x |
          ! | x x x |
          !  -------
          if((i.lt.size(grdpts_id,1)).and.(grdpts_id(i+1,j+1).eq.bc_pt)) then

             procedure_type = S_edge_type
             gradient_type  = gradient_I_type

          !  -------
          ! | 3 3 x |
          ! | x 0*x |
          ! | x x x |
          !  -------
          else

             procedure_type = S_edge_type
             gradient_type  = gradient_R0_type

          end if

        end subroutine get_newgrdpt_procedure_6_1


        subroutine error_newgrdpt_procedure(i,j,grdpts_id)

          implicit none

          integer(ikind)         , intent(in) :: i
          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(in) :: grdpts_id

          integer :: min_i
          integer :: max_i
          integer :: min_j
          integer :: max_j

          character(len=3) :: format_grdpt

          integer :: k

          min_i = max(1,i-1)
          max_i = min(size(grdpts_id,1),i+1)
          min_j = max(1,j-1)
          max_j = min(size(grdpts_id,2),j+1)

          write(format_grdpt,'(I1,''I2'')'), max_i-min_i+1

          do k=min_j,max_j

             print (format_grdpt), grdpts_id(min_i:max_i,max_j-(k-min_j))

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

      end module bf_layer_newgrdpt_procedure_module
