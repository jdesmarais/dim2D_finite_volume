      !module containing subroutines determing the 
      !type of procedure used at a boundary gridpoint
      module bf_layer_bc_procedure_module
      
c$$$        use bc_operators_class, only :
c$$$     $       bc_operators

        use parameters_bf_layer, only :
     $     bc_interior_pt,
     $     bc_pt

        implicit none

        integer, parameter :: no_bc_procedure_type=0
        integer, parameter :: SW_corner_type=1
        integer, parameter :: SE_corner_type=2
        integer, parameter :: NW_corner_type=3
        integer, parameter :: NE_corner_type=4
        integer, parameter :: S_edge_type=5
        integer, parameter :: E_edge_type=6
        integer, parameter :: W_edge_type=7
        integer, parameter :: N_edge_type=8
        integer, parameter :: SE_edge_type=9
        integer, parameter :: SW_edge_type=10
        integer, parameter :: NE_edge_type=11
        integer, parameter :: NW_edge_type=12

        integer, parameter :: test1=1
        integer, parameter :: test2=2
        integer, parameter :: test3=3
        integer, parameter :: test4=4
        integer, parameter :: test51=51
        integer, parameter :: test52=52
        integer, parameter :: test61=61
        integer, parameter :: test62=62        

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
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       get_bc_interior_pt_procedure,
     $       get_bc_pt_procedure,
     $       is_bc_pt_procedure_of_edge_type,
     $       does_corner_procedure_1_match,
     $       does_corner_procedure_2_match,
     $       does_corner_procedure_3_match,
     $       does_corner_procedure_4_match,
     $       does_corner_procedure_5_1_match,
     $       does_corner_procedure_5_2_match,
     $       does_corner_procedure_6_1_match,
     $       does_corner_procedure_6_2_match

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
                   procedure_type=SW_edge_type

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
                      procedure_type=SE_edge_type
                      
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


        !error pattern not found for bc_interior_pt
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


        subroutine get_bc_pt_procedure(
     $     i,j,
     $     grdpts_id,
     $     procedure_type,
     $     nbpts_x,
     $     nbpts_y)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nbpts_x
          integer                , intent(out) :: nbpts_y

          logical :: edge_type

          integer, dimension(2), parameter :: tests11 = [test4,test62]
          integer, dimension(4), parameter :: tests12 = [test2,test4,test52,test62]
          integer, dimension(3), parameter :: tests13 = [test2,test51,test52]

          integer, dimension(4), parameter :: tests21 = [test3,test4,test61,test62]
          integer, dimension(8), parameter :: tests22 = [test1,test2,test3,test4,test51,test52,test61,test62]
          integer, dimension(4), parameter :: tests23 = [test1,test2,test51,test52]

          integer, dimension(2), parameter :: tests31 = [test3,test61]
          integer, dimension(4), parameter :: tests32 = [test1,test3,test51,test61]
          integer, dimension(2), parameter :: tests33 = [test1,test51]

          logical :: match
          integer :: k

          !determine whether the procedure is of edge type or not
          edge_type = is_bc_pt_procedure_of_edge_type(
     $         i,j,
     $         grdpts_id,
     $         procedure_type)

          !if the procedure is not of edge type, we try to find
          !a match for a procedure of corner type
          if(.not.edge_type) then

          !i=1
             if(i.eq.1) then
                
                !j=1
                if(j.eq.1) then
                   
                   do k=1, size(tests11,1)
                      match = does_corner_procedure_match(
     $                     tests11(k),
     $                     i,j,
     $                     grdpts_id,
     $                     procedure_type,
     $                     nbpts_x,
     $                     nbpts_y)
                      if(match) then
                         exit
                      end if
                   end do

                else

                !1<j<ny
                   if(j.le.size(grdpts_id,2)-1) then

                      do k=1, size(tests12,1)
                         match = does_corner_procedure_match(
     $                        tests12(k),
     $                        i,j,
     $                        grdpts_id,
     $                        procedure_type,
     $                        nbpts_x,
     $                        nbpts_y)
                         if(match) then
                            exit
                         end if
                      end do

                !j=ny
                   else

                      do k=1, size(tests13,1)
                         match = does_corner_procedure_match(
     $                        tests13(k),
     $                        i,j,
     $                        grdpts_id,
     $                        procedure_type,
     $                        nbpts_x,
     $                        nbpts_y)
                         if(match) then
                            exit
                         end if
                      end do                      

                   end if
                end if

             else

          !1<i<nx
                if(i.le.size(grdpts_id,1)-1) then

                !j=1
                   if(j.eq.1) then

                      do k=1, size(tests21,1)
                         match = does_corner_procedure_match(
     $                        tests21(k),
     $                        i,j,
     $                        grdpts_id,
     $                        procedure_type,
     $                        nbpts_x,
     $                        nbpts_y)
                         if(match) then
                            exit
                         end if
                      end do

                   else

                !1<j<ny
                      if(j.lt.size(grdpts_id,2)) then

                         do k=1, size(tests22,1)
                            match = does_corner_procedure_match(
     $                           tests22(k),
     $                           i,j,
     $                           grdpts_id,
     $                           procedure_type,
     $                           nbpts_x,
     $                           nbpts_y)
                            if(match) then
                               exit
                            end if
                         end do

                !j=ny
                      else

                         do k=1, size(tests23,1)
                            match = does_corner_procedure_match(
     $                           tests23(k),
     $                           i,j,
     $                           grdpts_id,
     $                           procedure_type,
     $                           nbpts_x,
     $                           nbpts_y)
                            if(match) then
                               exit
                            end if
                         end do

                      end if

                   end if

          !i=nx
                else

                   !j=1
                   if(j.eq.1) then

                      do k=1, size(tests31,1)
                         match = does_corner_procedure_match(
     $                        tests31(k),
     $                        i,j,
     $                        grdpts_id,
     $                        procedure_type,
     $                        nbpts_x,
     $                        nbpts_y)
                         if(match) then
                            exit
                         end if
                      end do

                   else

                !1<j<ny
                      if(j.lt.size(grdpts_id,2)) then

                         do k=1, size(tests32,1)
                            match = does_corner_procedure_match(
     $                           tests32(k),
     $                           i,j,
     $                           grdpts_id,
     $                           procedure_type,
     $                           nbpts_x,
     $                           nbpts_y)
                            if(match) then
                               exit
                            end if
                         end do

                !j=ny
                      else

                         do k=1, size(tests33,1)
                            match = does_corner_procedure_match(
     $                           tests33(k),
     $                           i,j,
     $                           grdpts_id,
     $                           procedure_type,
     $                           nbpts_x,
     $                           nbpts_y)
                            if(match) then
                               exit
                            end if
                         end do

                      end if

                   end if
                   
                end if
                
             end if
             
          end if

        end subroutine get_bc_pt_procedure


        function is_bc_pt_procedure_of_edge_type(
     $     i,j,grdpts_id,procedure_type)
     $     result(edge_type)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          logical                              :: edge_type

          integer :: k


          !first determine whether the procedure can be of
          !N_edge or of S_edge type
          !for this, there should be 2 pts on both sides of
          !the current grid point in the x-direction of type
          !bc_pt
          if((i.ge.3).and.(i.le.(size(grdpts_id,1)-2))) then
             
             edge_type = .true.

             do k=-2,-1
                edge_type = edge_type.and.(grdpts_id(k+i,j).eq.bc_pt)
             end do

             if(edge_type) then

                do k=1,2
                   edge_type = edge_type.and.(grdpts_id(i+k,j).eq.bc_pt)
                end do

                if(edge_type) then

                   if(j.eq.1) then
                      procedure_type = S_edge_type
                   else
                      if(grdpts_id(i,j-1).eq.bc_interior_pt) then
                         procedure_type = N_edge_type
                      else
                         procedure_type = S_edge_type
                      end if
                   end if

                end if
             end if

          else

             edge_type=.false.

          end if


          !if the procedure is not of edge_type in the x-direction
          !determine whether the procedure can be of E_edge or of
          !W_edge type
          !for this, there should be 2 pts on both sides of
          !the current grid point in the y-direction of type
          !bc_pt
          if(.not.edge_type) then

             if((j.ge.3).and.(j.le.(size(grdpts_id,2)-2))) then

                edge_type = .true.

                do k=-2,-1
                   edge_type = edge_type.and.(grdpts_id(i,k+j).eq.bc_pt)
                end do

                if(edge_type) then

                   do k=1,2
                      edge_type = edge_type.and.(grdpts_id(i,j+k).eq.bc_pt)
                   end do

                   if(edge_type) then

                      if(i.eq.1) then
                         procedure_type = W_edge_type
                      else

                         if(grdpts_id(i-1,j).eq.bc_interior_pt) then
                            procedure_type = E_edge_type
                         else
                            procedure_type = W_edge_type
                         end if

                      end if

                   end if
                end if
             end if
          end if

        end function is_bc_pt_procedure_of_edge_type


        function does_corner_procedure_match(
     $     test_id, i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)
     $     result(match)

          implicit none
          
          integer                , intent(in)  :: test_id
          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match

          select case(test_id)
            case(test1)
               match = does_corner_procedure_1_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test2)
               match = does_corner_procedure_2_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test3)
               match = does_corner_procedure_3_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test4)
               match = does_corner_procedure_4_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test51)
               match = does_corner_procedure_5_1_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test52)
               match = does_corner_procedure_5_2_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test61)
               match = does_corner_procedure_6_1_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case(test62)
               match = does_corner_procedure_6_2_match(
     $              i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

            case default
               print '(''bf_layer_bc_procedure'')'
               print '(''does_corner_procedure_match'')'
               print '(''case : '',I2)', test_id
               stop 'case not recognized'

          end select

        end function does_corner_procedure_match




        !  -------       -------        ------- 
        ! |       |     |       |      |       |
        ! |   1*  |  -> |   1*  |  or  | 1 1*  |
        ! | 1     |     | 1 1   |      | 1     |
        !  -------       -------        ------- 
        !-----------------------------------------
        function does_corner_procedure_1_match(
     $     i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if(grdpts_id(i-1,j-1).eq.bc_pt) then

             if(grdpts_id(i,j-1).eq.bc_pt) then

                match = does_corner_procedure_1_1_match(
     $               i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

             else

                match = does_corner_procedure_1_2_match(
     $               i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)

             end if

          else
             match = .false.
          end if

        end function does_corner_procedure_1_match


        !  -------       ------- 
        ! |       |     |       |
        ! |   1*  |  -> |   1*  |
        ! | 1     |     | 1 1   |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_1_1_match(
     $     i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          !  ------- 
          ! |       |
          ! | 0 1*  | SE_corner, [0,1]
          ! | 1 1   |
          !  ------- 
          !----------------------
          if(grdpts_id(i-1,j).eq.bc_interior_pt) then
             procedure_type = SE_corner_type
             nb_pts_x       = 0
             nb_pts_y       = 1
             match          = .true.
             
          else
             if(i.lt.size(grdpts_id,1)) then

          !  ------- 
          ! |       |
          ! |   1*1 | NW_corner, [0,0]
          ! | 1 1   |
          !  ------- 
          !----------------------
                if(grdpts_id(i+1,j).eq.bc_pt) then
                   procedure_type = NW_corner_type
                   nb_pts_x       = 0
                   nb_pts_y       = 0
                   match          = .true.

                else
                   
                   if(j.lt.size(grdpts_id,2)) then
                      
          !  ------- 
          ! |   1 1 |
          ! |   1*  | NW_corner, [0,0]
          ! | 1 1   |
          !  ------- 
          !----------------------                            
                      if(grdpts_id(i+1,j+1).eq.bc_pt) then
                         procedure_type = NW_corner_type
                         nb_pts_x       = 0
                         nb_pts_y       = 1
                         match          = .true.
          !  _______
          ! |   1   |
          ! |   1   |
          ! |   1*  | edge_W [0]
          ! | 1 1   |
          !  ------- 
          !----------------------                     
                      else
                         procedure_type = W_edge_type
                         nb_pts_x       = 0
                         nb_pts_y       = 0
                         match          = .true.
                         
                      end if

                   else
                      call error_bc_pt_procedure(grdpts_id(i-1:i+1,j-1:j))
                   end if
                   
                end if
                
             else
                call error_bc_pt_procedure(grdpts_id(i-1:i,j-1:j))
             end if
          end if

        end function does_corner_procedure_1_1_match


        !  -------       ------- 
        ! |       |     |       |
        ! |   1*  |  -> | 1 1*  |
        ! | 1     |     | 1 x   |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_1_2_match(
     $     i, j, grdpts_id, procedure_type, nb_pts_x, nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if(grdpts_id(i-1,j).eq.bc_pt) then

          !  ------- 
          ! |       |
          ! | 1 1*  | NW_corner, [0,1]
          ! | 1 0 x |
          !  ------- 
          !----------------------
             if(grdpts_id(i,j-1).eq.bc_interior_pt) then
                procedure_type = NW_corner_type
                nb_pts_x       = 1
                nb_pts_y       = 0
                match          = .true.
                
             else
                if(j.lt.size(grdpts_id,2)) then

          !  ------- 
          ! |   1   |
          ! | 1 1*  | SE_corner, [0,0]
          ! | 1   x |
          !  ------- 
          !----------------------
                   if(grdpts_id(i,j+1).eq.bc_pt) then
                      procedure_type = SE_corner_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
                      
          !  ------- 
          ! |     1 |
          ! | 1 1*1 | SE_corner, [1,0]
          ! | 1   x |
          !  ------- 
          !----------------------
                   else
                      if(grdpts_id(i+1,j+1).eq.bc_pt) then
                         procedure_type = SE_corner_type
                         nb_pts_x       = 1
                         nb_pts_y       = 0
                         match          = .true.
          !  ---------
          ! |         |
          ! | 1 1*1 1 | S_edge, [0]
          ! | 1   x   |
          !  ---------
          !----------------------
                      else
                         procedure_type = S_edge_type
                         nb_pts_x       = 0
                         nb_pts_y       = 0 
                         match          = .true.                              
                         
                      end if
                            
                   end if

                else
                   call error_bc_pt_procedure(grdpts_id(i-1:i,j-1:j))
                end if
                
             end if
          else
             call error_bc_pt_procedure(grdpts_id(i-1:i,j-1:j))
          end if
          
        end function does_corner_procedure_1_2_match


        !  -------       -------        ------- 
        ! |       |     |       |      |       |
        ! |   1*  |  -> |   1*  |  or  |   1*1 |
        ! | x   1 |     | x 1 1 |      | x   1 |
        !  -------       -------        ------- 
        !-----------------------------------------
        function does_corner_procedure_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          
          if(grdpts_id(i+1,j-1).eq.bc_pt) then

             if(grdpts_id(i,j-1).eq.bc_pt) then

                match = does_corner_procedure_2_1_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             else

                match = does_corner_procedure_2_2_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             end if

          else
             match = .false.
          end if

        end function does_corner_procedure_2_match

      
        !  -------       ------- 
        ! |       |     |       |
        ! |   1*  |  -> |   1*  |
        ! | x   1 |     | x 1 1 |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_2_1_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          !  ------- 
          ! |       |
          ! |   1*0 | SW_corner, [0,1]
          ! | x 1 1 |
          !  ------- 
          !----------------------
          if(grdpts_id(i+1,j).eq.bc_interior_pt) then
             procedure_type = SW_corner_type
             nb_pts_x       = 0
             nb_pts_y       = 1
             match          = .true.

          else
                   
             if(i.gt.1) then
                      
          !  ------- 
          ! |       |
          ! | 1 1*  | NE_corner, [0,0]
          ! | x 1 1 |
          !  ------- 
          !----------------------
                if(grdpts_id(i-1,j).eq.bc_pt) then
                   
                   procedure_type = NE_corner_type
                   nb_pts_x       = 0
                   nb_pts_y       = 0
                   match          = .true.
                   
                else
                   if(j.lt.size(grdpts_id,2)) then
                      
          !  ------- 
          ! | 1 1   |
          ! | x 1*  | NE_corner, [0,1]
          ! | x 1 1 |
          !  ------- 
          !----------------------
                      if(grdpts_id(i-1,j+1).eq.bc_pt) then
                         
                         procedure_type = NE_corner_type
                         nb_pts_x       = 0
                         nb_pts_y       = 1
                         match          = .true.
                         
          !  ------- 
          ! |   1   |
          ! | x 1   |
          ! | x 1*  | E_edge, [0]
          ! | x 1 1 |
          !  ------- 
          !----------------------
                      else
                         
                         procedure_type = E_edge_type
                         nb_pts_x       = 0
                         nb_pts_y       = 0
                         match          = .true.
                         
                      end if
                      
                   else
                      call error_bc_pt_procedure(grdpts_id(i-1:i+1,j-1:j))                            
                   end if
                end if

             else
                call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j))
             end if
             
          end if

        end function does_corner_procedure_2_1_match


        !  -------       ------- 
        ! |       |     |       |
        ! |   1*  |  -> |   1*1 |
        ! | x   1 |     | x x 1 |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_2_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match
      
          
          if(grdpts_id(i+1,j).eq.bc_pt) then

          !  ------- 
          ! |       |
          ! |   1*1 | NE_corner, [1,0]
          ! | x 0 1 |
          !  ------- 
          !----------------------
             if(grdpts_id(i,j-1).eq.bc_interior_pt) then
                procedure_type = NE_corner_type
                nb_pts_x       = 1
                nb_pts_y       = 0
                match          = .true.
                
             else
                
                if(j.lt.size(grdpts_id,2)) then

          !  ------- 
          ! |   1   |
          ! |   1*1 | SW_corner, [0,0]
          ! | x   1 |
          !  ------- 
          !----------------------
                   if(grdpts_id(i,j+1).eq.bc_pt) then
                      procedure_type = SW_corner_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
                      
                   else
                      
                      if(i.gt.1) then
                         
          !  ------- 
          ! | 1     |
          ! | 1 1*1 | SW_corner, [1,0]
          ! | x   1 |
          !  ------- 
          !----------------------
                         if(grdpts_id(i-1,j+1).eq.bc_pt) then
                            procedure_type = SW_corner_type
                            nb_pts_x       = 1
                            nb_pts_y       = 0
                            match          = .true.
                            
          !  -------- 
          ! |  x     |
          ! |1 1 1*1 | S_edge, [0]
          ! |  x   1 |
          !  ------- 
          !----------------------
                         else
                            procedure_type = S_edge_type
                            nb_pts_x       = 0
                            nb_pts_y       = 0
                            match          = .true.
                            
                         end if
                         
                      else
                         call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j+1))
                      end if
                   end if
                   
                else
                   call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j))
                end if
             end if
          else
             call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j))
          end if

        end function does_corner_procedure_2_2_match
             


        !  -------       -------        ------- 
        ! | 1     |     | 1     |      | 1 1   |
        ! |   1*  |  -> | 1 1*  |  or  |   1*  |
        ! | x   x |     | x   x |      | x   x |
        !  -------       -------        ------- 
        !-----------------------------------------
        function does_corner_procedure_3_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match

          if(grdpts_id(i-1,j+1).eq.bc_pt) then

             if(grdpts_id(i-1,j).eq.bc_pt) then

                match = does_corner_procedure_3_1_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             else

                match = does_corner_procedure_3_2_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             end if

          else
             
             match = .false.

          end if

        end function does_corner_procedure_3_match


        !  -------       ------- 
        ! | 1     |     | 1     |
        ! |   1*  |  -> | 1 1*  |
        ! | x   x |     | x   x |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_3_1_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          !  -------
          ! | 1 0   |
          ! | 1 1*  | SW_corner, [1,0]
          ! | x   x |
          !  ------- 
          !----------------------
          if(grdpts_id(i,j+1).eq.bc_interior_pt) then
             procedure_type = SW_corner_type
             nb_pts_x       = 1
             nb_pts_y       = 0
             match          = .true.
             
          else

             if(j.gt.1) then

          !  -------
          ! | 1     |
          ! | 1 1*  | NE_corner, [0,0]
          ! | x 1 x |
          !  ------- 
          !----------------------
                if(grdpts_id(i,j-1).eq.bc_pt) then
                   procedure_type = NE_corner_type
                   nb_pts_x       = 0
                   nb_pts_y       = 0
                   match          = .true.
                      
          !  -------
          ! | 1     |
          ! | 1 1*1 | N_edge, [0]
          ! | x x x |
          !  ------- 
          !----------------------
                else
                   
                   if(i.lt.size(grdpts_id,1)) then
                      procedure_type = N_edge_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
                      
                   else
                      call error_bc_pt_procedure(grdpts_id(i-1:i,j-1:j+1))
                   end if
                   
                end if                
             else
                call error_bc_pt_procedure(grdpts_id(i-1:i,j:j+1))
             end if
          end if
             
        end function does_corner_procedure_3_1_match


        !  -------       ------- 
        ! | 1     |     | 1 1   |
        ! |   1*  |  -> |   1*  |
        ! | x   x |     | x   x |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_3_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if(grdpts_id(i,j+1).eq.bc_pt) then

          !  -------
          ! | 1 1   |
          ! |   1*  | NE_corner, [0,1]
          ! | x   x |
          !  ------- 
          !----------------------
             if(grdpts_id(i-1,j).eq.bc_interior_pt) then
                procedure_type = NE_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 1
                match          = .true.

             else
                if(i.lt.size(grdpts_id,1)) then
                   
          !  -------
          ! | 1 1   |
          ! |   1*1 | SW_corner, [0,0]
          ! | x   x |
          !  ------- 
          !----------------------
                   if(grdpts_id(i+1,j).eq.bc_pt) then
                      procedure_type = SW_corner_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
          !  -------
          ! | 1 1   |
          ! |   1*  | W_edge, [0]
          ! | x 1 x |
          !  ------- 
          !----------------------
                   else
                      if(j.gt.1) then
                         procedure_type = W_edge_type
                         nb_pts_x       = 0
                         nb_pts_y       = 0
                         match          = .true.
                         
                      else
                         call error_bc_pt_procedure(grdpts_id(i-1:i+1,j:j+1))
                      end if
                      
                   end if
                   
                else
                   call error_bc_pt_procedure(grdpts_id(i-1:i,j:j+1))
                end if
             end if
          else
             call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j))
          end if
          
        end function does_corner_procedure_3_2_match


        !  -------       -------        ------- 
        ! | x   1 |     | x   1 |      |   1 1 |
        ! |   1*  |  -> |   1*1 |  or  |   1*  |
        ! | x   x |     | x   x |      | x   x |
        !  -------       -------        ------- 
        !-----------------------------------------
        function does_corner_procedure_4_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match

          if(grdpts_id(i+1,j+1).eq.bc_pt) then

             if(grdpts_id(i+1,j).eq.bc_pt) then

                match = does_corner_procedure_4_1_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             else

                match = does_corner_procedure_4_2_match(
     $               i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)

             end if

          else

             match = .false.

          end if

        end function does_corner_procedure_4_match


        !  -------       ------- 
        ! | x   1 |     | x   1 |
        ! |   1*  |  -> |   1*1 |
        ! | x   x |     | x   x |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_4_1_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          !  -------
          ! | x 0 1 |
          ! |   1*1 | SE_corner, [1,0]
          ! | x   x |
          !  ------- 
          !----------------------
          if(grdpts_id(i,j+1).eq.bc_interior_pt) then
             procedure_type = SE_corner_type
             nb_pts_x       = 1
             nb_pts_y       = 0
             match          = .true.
             
          else
             
             if(j.gt.1) then
                
          !  -------
          ! | x   1 |
          ! |   1*1 | NW_corner, [0,0]
          ! | x 1 x |
          !  ------- 
          !----------------------
                if(grdpts_id(i,j-1).eq.bc_pt) then
                   procedure_type = NW_corner_type
                   nb_pts_x       = 0
                   nb_pts_y       = 0
                   match          = .true.
                   
                else
                   
          !  -------
          ! | x   1 |
          ! | 1 1*1 | N_edge, [0]
          ! | x x x |
          !  ------- 
          !----------------------
                   if(i.gt.1) then
                      procedure_type = N_edge_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
                      
                   else
                      call error_bc_pt_procedure(grdpts_id(i:i+1,j-1:j+1))
                   end if
                   
                end if
                
             else
                call error_bc_pt_procedure(grdpts_id(i:i+1,j:j+1))
             end if
          end if

        end function does_corner_procedure_4_1_match


        !  -------       ------- 
        ! | x   1 |     | x 1 1 |
        ! |   1*  |  -> |   1*x |
        ! | x   x |     | x   x |
        !  -------       ------- 
        !-----------------------------------------
        function does_corner_procedure_4_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match

          
          if(grdpts_id(i,j+1).eq.bc_pt) then

          !  -------
          ! | x 1 1 |
          ! |   1*0 | NW_corner, [0,1]
          ! | x   x |
          !  ------- 
          !----------------------
             if(grdpts_id(i+1,j).eq.bc_interior_pt) then
                procedure_type = NW_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 1
                match          = .true.
                
             else
                      
                if(i.gt.1) then
          !  -------
          ! | x 1 1 |
          ! | 1 1*x | SE_corner, [0,0]
          ! | x   x |
          !  ------- 
          !----------------------
                   if(grdpts_id(i-1,j).eq.bc_pt) then
                      procedure_type = SE_corner_type
                      nb_pts_x       = 0
                      nb_pts_y       = 0
                      match          = .true.
                      
                   else
                      
          !  -------
          ! | x 1 1 |
          ! | x 1*x | E_edge, [0]
          ! | x 1 x |
          !  ------- 
          !----------------------
                      if(j.gt.1) then
                         procedure_type = E_edge_type
                         nb_pts_x       = 0
                         nb_pts_y       = 0
                         match          = .true.
                         
                      else
                         call error_bc_pt_procedure(grdpts_id(i-1:i+1,j:j+1))
                      end if
                      
                   end if
                else
                   call error_bc_pt_procedure(grdpts_id(i:i+1,j:j+1))
                end if
             end if
             
          else
             call error_bc_pt_procedure(grdpts_id(i:i+1,j:j+1))
          end if
       
        end function does_corner_procedure_4_2_match




        !  -------       ------- 
        ! | x   x |     | x   x |
        ! |   1*  |  -> | 1 1*  |
        ! | x 1 x |     | x 1 x |
        !  -------       ------- 
        !------------------------
        function does_corner_procedure_5_1_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if((grdpts_id(i,j-1).eq.bc_pt).and.(grdpts_id(i-1,j).eq.bc_pt)) then

          !  -------
          ! | x   x |
          ! | 1 1*  | NE_corner, [0,0]
          ! | 0 1 x |
          !  ------- 
          !----------------------
             if(grdpts_id(i-1,j-1).eq.bc_interior_pt) then
                procedure_type = NE_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.
                
             else

          !  -------
          ! | x   x |
          ! | 1 1*  | SW_edge, [0,0]
          ! |   1 x |
          !  ------- 
          !----------------------
                procedure_type = SW_edge_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.

             end if
             
          else
             match = .false.
          end if

        end function does_corner_procedure_5_1_match


        !  -------       ------- 
        ! | x   x |     | x   x |
        ! |   1*  |  -> |   1*1 |
        ! | x 1 x |     | x 1 x |
        !  -------       ------- 
        !------------------------
        function does_corner_procedure_5_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if((grdpts_id(i,j-1).eq.bc_pt).and.(grdpts_id(i+1,j).eq.bc_pt)) then

          !  -------
          ! | x   x |
          ! |   1*1 | NW_corner, [0,0]
          ! | x 1 x |
          !  ------- 
          !----------------------
             if(grdpts_id(i+1,j-1).eq.bc_interior_pt) then
                procedure_type = NW_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.
                
             else

          !  -------
          ! | x   x |
          ! |   1*1 | SE_edge, [0]
          ! | x 1 x |
          !  ------- 
          !----------------------
                procedure_type = SE_edge_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.

             end if
             
          else
             match = .false.
          end if

        end function does_corner_procedure_5_2_match


        !  -------       ------- 
        ! | x 1 x |     | x 1 x |
        ! |   1*  |  -> | 1 1*  |
        ! | x x x |     | x x x |
        !  -------       ------- 
        !------------------------
        function does_corner_procedure_6_1_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if((grdpts_id(i,j+1).eq.bc_pt).and.(grdpts_id(i-1,j).eq.bc_pt)) then

          !  -------
          ! | 0 1 x |
          ! | 1 1*  | SE_corner, [0,0]
          ! | x x x |
          !  ------- 
          !----------------------
             if(grdpts_id(i-1,j+1).eq.bc_interior_pt) then
                procedure_type = SE_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.
                
             else

          !  -------
          ! | x 1 x |
          ! | 1 1*  | NW_edge, [0,0]
          ! | x x x |
          !  ------- 
          !----------------------
                procedure_type = NW_edge_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.

             end if
             
          else
             match = .false.
          end if

        end function does_corner_procedure_6_1_match


        !  -------       ------- 
        ! | x 1 x |     | x 1 x |
        ! |   1*  |  -> | x 1*1 |
        ! | x x x |     | x x x |
        !  -------       ------- 
        !------------------------
        function does_corner_procedure_6_2_match(
     $     i,j,grdpts_id,procedure_type,nb_pts_x,nb_pts_y)
     $     result(match)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: nb_pts_x
          integer                , intent(out) :: nb_pts_y
          logical                              :: match


          if((grdpts_id(i,j+1).eq.bc_pt).and.(grdpts_id(i+1,j).eq.bc_pt)) then

          !  -------
          ! | x 1 0 |
          ! | x 1*1 | SW_corner, [0,0]
          ! | x x x |
          !  ------- 
          !----------------------
             if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then
                procedure_type = SW_corner_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.
                
             else

          !  -------
          ! | x 1 x |
          ! | x 1*1 | NE_edge, [0,0]
          ! | x x x |
          !  ------- 
          !----------------------
                procedure_type = NE_edge_type
                nb_pts_x       = 0
                nb_pts_y       = 0
                match          = .true.

             end if
             
          else
             match = .false.
          end if

        end function does_corner_procedure_6_2_match


        !error pattern not found for bc_interior_pt
        subroutine error_bc_pt_procedure(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(in) :: grdpts_id

          print '(''bf_bc_procedure_module'')'
          print '(''get_bc_interior_pt_procedure'')'
          print '(''grdpts_id'')'
          print '(2I2)', grdpts_id(:,2)
          print '(2I2)', grdpts_id(:,1)
          stop 'case not recognized'

        end subroutine error_bc_pt_procedure

      end module bf_layer_bc_procedure_module
