      !module containing subroutines determing the 
      !type of procedure used at a boundary gridpoint
      module bf_layer_bc_procedure_module
      
        use parameters_bf_layer, only :
     $     bc_interior_pt,
     $     bc_pt

        implicit none

        integer, parameter :: SW_corner_type=0
        integer, parameter :: SE_corner_type=1
        integer, parameter :: NW_corner_type=2
        integer, parameter :: NE_corner_type=3
        integer, parameter :: S_edge_type=4
        integer, parameter :: E_edge_type=5
        integer, parameter :: W_edge_type=6
        integer, parameter :: N_edge_type=7

        private
        public ::
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       N_edge_type,
     $       get_bc_interior_pt_procedure

        contains

        !the grid point (i,j) has been identified as a bc_interior_pt
        !we want to decide which open boundary procedure should be
        !applied at this grid point
        !1) is it a corner or an edge procedure ?
        !2.1) if it is a corner procedure, what is the orientation of
        !     the corner to know in which directions are the waves
        !     incoming ?
        !2.2) if it is an edge procedure, in which direction are the
        !     time derivatives evaluated using the fluxes and in the
        !     transverse direction, in which direction are the waves
        !     incoming ?
        function get_bc_interior_pt_procedure(i,j,grdpts_id)
     $       result(procedure_type)

          implicit none

          integer                , intent(in) :: i
          integer                , intent(in) :: j
          integer, dimension(:,:), intent(in) :: grdpts_id
          integer                             :: procedure_type


          if(grdpts_id(i,j-1).eq.bc_interior_pt) then
           
             !  -------
             ! |       |
             ! | 0 0   |
             ! |   0   |
             !  -------
             if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                !SW corner
                if(grdpts_id(i-1,j-1).eq.bc_pt) then
                   procedure_type=SW_corner_type

                !NE corner
                else
                   procedure_type=NE_corner_type
                end if

             else

             !  -------
             ! |       |
             ! | X 0 0 |
             ! |   0   |
             !  -------
                if(grdpts_id(i+1,j).eq.bc_interior_pt) then
                   
                   !SE corner
                   if(grdpts_id(i+1,j-1).eq.bc_pt) then
                      procedure_type=SE_corner_type
                      
                   !NW corner
                   else
                      procedure_type=NW_corner_type
                   end if

             
                else

             !  -------
             ! |   0   |
             ! | X 0 X |
             ! |   0   |
             !  -------
                   if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                      !E edge
                      if(grdpts_id(i+1,j).eq.bc_pt) then
                         procedure_type=E_edge_type

                      !W edge
                      else
                         procedure_type=W_edge_type
                      end if

             !  -------
             ! |   X   |
             ! | X 0 X |
             ! |   0   |
             !  -------
                   else

                      call error_bc_interior_pt_procedure(i,j,grdpts_id)

                   end if
                end if
             end if

          else

             if(grdpts_id(i-1,j).eq.bc_interior_pt) then

             !  -------
             ! |       |
             ! | 0 0 0 |
             ! |   X   |
             !  -------
                if(grdpts_id(i+1,j).eq.bc_interior_pt) then
                   
                  !S edge
                  if(grdpts_id(i,j-1).eq.bc_pt) then
                     procedure_type=S_edge_type

                  !N edge
                  else
                     procedure_type=N_edge_type
                  end if

                else

             !  -------
             ! |   0   |
             ! | 0 0 X |
             ! |   X   |
             !  -------                      
                   if(grdpts_id(i,j+1).eq.bc_interior_pt) then
                      
                      !NW corner
                      if(grdpts_id(i-1,j+1).eq.bc_pt) then
                         procedure_type=NW_corner_type
                      !SE corner
                      else
                         procedure_type=SE_corner_type
                      end if

             !  -------
             ! |   X   |
             ! | 0 0 X |
             ! |   X   |
             !  -------
                   else
                      
                      call error_bc_interior_pt_procedure(i,j,grdpts_id)

                   end if

                end if

             else

                if(grdpts_id(i+1,j).eq.bc_interior_pt) then
                   
             !  -------
             ! |   0   |
             ! | X 0 0 |
             ! |   X   |
             !  -------
                   if(grdpts_id(i,j+1).eq.bc_interior_pt) then
                      
                      !SW corner
                      if(grdpts_id(i-1,j-1).eq.bc_pt) then
                         procedure_type=SW_corner_type
                      !NE corner
                      else
                         procedure_type=NE_corner_type
                      end if

                   else

                      call error_bc_interior_pt_procedure(i,j,grdpts_id)

                   end if
                end if                
             end if
          end if

        end function get_bc_interior_pt_procedure


        subroutine error_bc_interior_pt_procedure(i,j,grdpts_id)

          implicit none

          integer                , intent(in) :: i
          integer                , intent(in) :: j
          integer, dimension(:,:), intent(in) :: grdpts_id

          print '(''bf_bc_procedure_module'')'
          print '(''get_bc_interior_pt_procedure'')'
          print '(''grdpts(i-1:i+1,j-1:j+1)'')'
          print '(3I2)', grdpts_id(i-1:i+1,j+1)
          print '(3I2)', grdpts_id(i-1:i+1,j)
          print '(3I2)', grdpts_id(i-1:i+1,j-1)
          stop 'case not recognized'

        end subroutine error_bc_interior_pt_procedure

      end module bf_layer_bc_procedure_module
