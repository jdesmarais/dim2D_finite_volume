      !> @file
      !> subroutines to determine the type of boundary procedure
      !> that should be used at the edge of the computational domain
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines to determine the type of boundary procedure
      !> that should be used at the edge of the computational domain
      !
      !> @date
      ! 28_02_2015 - documentation update  - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_bc_procedure_module
      
        use parameters_bf_layer, only :
     $     
     $     interior_pt,
     $     bc_interior_pt,
     $     bc_pt,
     $     
     $     BF_SUCCESS,
     $     
     $     SW_corner_type,
     $     SE_corner_type,
     $     NW_corner_type,
     $     NE_corner_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     N_edge_type,
     $     SE_edge_type,
     $     SW_edge_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     
     $     bc_procedure_extra_checks

        implicit none

        private
        public ::
     $       get_bc_interior_pt_procedure

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the type of boudnary procedure that should be
        !> applied at the edge of the computational domain
        !
        !> the grid point (i,j) has been identified as a bc_interior_pt
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
        !
        !> @date
        !> 28_02_2015 - initial version - J.L. Desmarais
        !
        !> @param i
        !> x-index identifying the position of the bc_interior_pt in
        !> grdpts_id
        !
        !> @param j
        !> y-index identifying the position of the bc_interior_pt in
        !> grdpts_id
        !
        !> @param grdpts_id
        !> role of the grid points
        !
        !> @param procedure_type
        !> type of procedure identifying the boundary
        !
        !> @param i_proc
        !> x-index identifying where the boundary is located
        !
        !> @param j_proc
        !> y-index identifying where the boundary is located
        !
        !> @param ierror
        !> integer identifying whether the identification of the
        !> boundary procedure was successful or not
        !--------------------------------------------------------------
        subroutine get_bc_interior_pt_procedure(
     $       i,
     $       j,
     $       grdpts_id,
     $       procedure_type,
     $       i_proc,
     $       j_proc,
     $       ierror)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer                , intent(out) :: procedure_type
          integer                , intent(out) :: i_proc
          integer                , intent(out) :: j_proc
          logical                , intent(out) :: ierror


          ierror = BF_SUCCESS


          !==============================
          !procedure 1
          !==============================
          !  -------
          ! |       |
          ! |   2*  |
          ! | 2     |
          !  -------
          if(grdpts_id(i-1,j-1).eq.bc_interior_pt) then
           
          !------------------------------
          !procedure 1.1
          !------------------------------
          !  -------
          ! |       |
          ! |   2*  |
          ! | 2 2   |
          !  -------
             if(grdpts_id(i,j-1).eq.bc_interior_pt) then

                call get_bc_interior_pt_procedure_1_1(
     $               i,j,grdpts_id,
     $               procedure_type,i_proc,j_proc,
     $               ierror)

             else

          !------------------------------
          !procedure 1.2
          !------------------------------
          !  -------
          ! |       |
          ! | 2 2*  |
          ! | 2 x   |
          !  -------
                if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                   call get_bc_interior_pt_procedure_1_2(
     $                  i,j,grdpts_id,
     $                  procedure_type,i_proc,j_proc,
     $                  ierror)
             
                else

          !------------------------------
          !procedure 1.3
          !------------------------------
          !  -------
          ! |       |
          ! | x 2*  |
          ! | 2 x   |
          !  -------
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                      
                end if
             end if


          !==============================
          !procedure 2
          !==============================
          !  -------
          ! |       |
          ! |   2*  |
          ! | x 2   |
          !  -------
          else

             if(grdpts_id(i,j-1).eq.bc_interior_pt) then

          !------------------------------
          !procedure 2.1
          !------------------------------
          !  -------
          ! |       |
          ! |   2*  |
          ! | x 2 2 |
          !  -------
                if(grdpts_id(i+1,j-1).eq.bc_interior_pt) then

                   call get_bc_interior_pt_procedure_2_1(
     $                  i,j,grdpts_id,
     $                  procedure_type,i_proc,j_proc,
     $                  ierror)

                else

          !------------------------------
          !procedure 2.2
          !------------------------------
          !  -------
          ! |       |
          ! | 2 2*  |
          ! | x 2 x |
          !  -------
                   if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                      call get_bc_interior_pt_procedure_2_2(
     $                     i,j,grdpts_id,
     $                     procedure_type,i_proc,j_proc,
     $                     ierror)

          !------------------------------
          !procedure 2.3
          !------------------------------
          !  -------
          ! |       |
          ! | x 2*2 |
          ! | x 2 x |
          !  -------
                   else

                      if(grdpts_id(i+1,j).eq.bc_interior_pt) then

                         call get_bc_interior_pt_procedure_2_3(
     $                        i,j,grdpts_id,
     $                        procedure_type,i_proc,j_proc,
     $                        ierror)
                         
                      else

          !------------------------------
          !procedure 2.4
          !------------------------------
          !  -------
          ! | 2     |
          ! | x 2*x |
          ! | x 2 x |
          !  -------
                         if(grdpts_id(i-1,j+1).eq.bc_interior_pt) then

                            !  -------
                            ! | 2 2   |
                            ! | x 2*x |
                            ! | x 2 x |
                            !  -------
                            if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                               call get_bc_interior_pt_procedure_2_4(
     $                              i,j,grdpts_id,
     $                              procedure_type,i_proc,j_proc,
     $                              ierror)

                            else

                               call error_bc_interior_pt_procedure(
     $                              i,j,grdpts_id,ierror)

                            end if
                               
                         else
          !------------------------------
          !procedure 2.5
          !------------------------------
          !  -------
          ! | x 2   |
          ! | x 2*x |
          ! | x 2 x |
          !  -------
                            if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                               call get_bc_interior_pt_procedure_2_5(
     $                              i,j,grdpts_id,
     $                              procedure_type,i_proc,j_proc,
     $                              ierror)
                               
                            else
                               
                               call error_bc_interior_pt_procedure(
     $                              i,j,grdpts_id,ierror)
                               
                            end if

                         end if
                         
                      end if
                      
                   end if

                end if

             else
                
          !==============================
          !procedure 3
          !==============================
          !  -------
          ! |       |
          ! |   2*  |
          ! | x x 2 |
          !  -------

                if(grdpts_id(i+1,j-1).eq.bc_interior_pt) then

          !------------------------------
          !procedure 3.1
          !------------------------------
          !  -------
          ! |       |
          ! |   2*2 |
          ! | x x 2 |
          !  -------
                   if(grdpts_id(i+1,j).eq.bc_interior_pt) then

                      call get_bc_interior_pt_procedure_3_1(
     $                     i,j,grdpts_id,
     $                     procedure_type,i_proc,j_proc,
     $                     ierror)

          !------------------------------
          !procedure 3.2
          !------------------------------
          !  -------
          ! |       |
          ! |   2*x |
          ! | x x 2 |
          !  -------
                   else

                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if

                else

          !==============================
          !procedure 4
          !==============================
          !  -------
          ! |       |
          ! | 2 2*  |
          ! | x x x |
          !  -------
                   if(grdpts_id(i-1,j).eq.bc_interior_pt) then

          !------------------------------
          !procedure 4.1
          !------------------------------
          !  -------
          ! |       |
          ! | 2 2*2 |
          ! | x x x |
          !  -------
                      if(grdpts_id(i+1,j).eq.bc_interior_pt) then

                         call get_bc_interior_pt_procedure_4_1(
     $                        i,j,grdpts_id,
     $                        procedure_type,i_proc,j_proc,
     $                        ierror)

          
          
                      else
          !------------------------------
          !procedure 4.2
          !------------------------------
          !  -------
          ! |   2   |
          ! | 2 2*x |
          ! | x x x |
          !  -------
                         if(grdpts_id(i,j+1).eq.bc_interior_pt) then
                            
                            call get_bc_interior_pt_procedure_4_2(
     $                           i,j,grdpts_id,
     $                           procedure_type,i_proc,j_proc,
     $                           ierror)

          !------------------------------
          !procedure 4.3
          !------------------------------
          !  -------
          ! |   x   |
          ! | 2 2*x |
          ! | x x x |
          !  -------
                         else

                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)

                         end if

                      end if

                   else                

          !==============================
          !procedure 5
          !==============================
          !  -------
          ! |       |
          ! | x 2*2 |
          ! | x x x |
          !  -------
                      if(grdpts_id(i+1,j).eq.bc_interior_pt) then

          !------------------------------
          !procedure 5.1
          !------------------------------
          !  -------
          ! |   2   |
          ! | x 2*2 |
          ! | x x x |
          !  -------
                         if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                            call get_bc_interior_pt_procedure_5_1(
     $                           i,j,grdpts_id,
     $                           procedure_type,i_proc,j_proc,
     $                           ierror)

          !------------------------------
          !procedure 5.2
          !------------------------------
          !  -------
          ! |   x   |
          ! | x 2*2 |
          ! | x x x |
          !  -------                            
                         else

                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror) 
                            
                         end if

          !==============================
          !procedure 6
          !==============================
          !  -------
          ! |       |
          ! | x 2*x |
          ! | x x x |
          !  -------
                      else

                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)

                      end if
                   end if
                end if
             end if
          end if

        end subroutine get_bc_interior_pt_procedure


        subroutine get_bc_interior_pt_procedure_1_1(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror

          !  -------
          ! |       |
          ! | 3 2*  |
          ! | 2 2   |
          !  -------
          if(grdpts_id(i-1,j).eq.bc_pt) then

          !  -------
          ! |   2   |
          ! | 3 2*  |
          ! | 2 2   |
          !  -------             
             if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                if(bc_procedure_extra_checks) then
                   !  --------
                   ! | 3 2 1/2|
                   ! | 3 2*1  |
                   ! | 2 2 1  |
                   !  --------
                   if(.not.(
     $                  (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.interior_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.bc_pt).and.(
     $                  (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $                  (grdpts_id(i+1,j+1).eq.interior_pt)))) then
                      
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                      
                   end if

                end if
               
                procedure_type = NW_edge_type
                i_proc         = i-1
                j_proc         = j-1
                
             else

                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i,j+1).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! |   3   |
          ! | 3 2*2 |
          ! | 2 2   |
          !  ------- 
                if(grdpts_id(i+1,j).eq.bc_interior_pt) then

                   if(bc_procedure_extra_checks) then
                   !  --------
                   ! | 3 3 3/2|
                   ! | 3 2*2  |
                   ! | 2 2 1  |
                   !  --------
                      if(.not.(
     $                     (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $                     (grdpts_id(i-1,j+1).eq.bc_pt).and.(
     $                     (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $                     (grdpts_id(i+1,j+1).eq.bc_pt)))) then
                         
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)
                         
                      end if
                   end if
                      
                   procedure_type = NW_corner_type
                   i_proc         = i-1
                   j_proc         = j
                      
                else

          !  -------
          ! |   x   |
          ! | 3 2*x |
          ! | 2 2   |
          !  ------- 
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

          else

             if(bc_procedure_extra_checks) then
                if(grdpts_id(i-1,j).ne.interior_pt) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
             
          !  -------
          ! |       |
          ! | 1 2*2 |
          ! | 2 2 3 |
          !  -------
             if(grdpts_id(i+1,j).eq.bc_interior_pt) then

                !  --------
                ! | 1 1 1/2|
                ! | 1 2*2  |
                ! | 2 2 3  |
                !  ---------
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.interior_pt).and.
     $                  (grdpts_id(i,j+1).eq.interior_pt).and.(
     $                  (grdpts_id(i+1,j+1).eq.interior_pt).or.
     $                  (grdpts_id(i+1,j+1).eq.bc_interior_pt)))) then
                
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if                      
                end if              

                procedure_type = SE_edge_type
                i_proc         = i
                j_proc         = j-1

             else

                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i+1,j).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! |   2   |
          ! | 1 2*3 |
          ! | 2 2 3 |
          !  -------
                if(grdpts_id(i,j+1).eq.bc_interior_pt) then

          !  -------
          ! |   2 2 |
          ! | 1 2*3 |
          ! | 2 2 3 |
          !  -------
                   if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then

                      !  -------
                      ! | 1 2 2 |
                      ! | 1 2*3 |
                      ! | 2 2 3 |
                      !  -------
                      if(bc_procedure_extra_checks) then
                         if(grdpts_id(i-1,j+1).ne.interior_pt) then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = SE_edge_type
                      i_proc         = i
                      j_proc         = j

          !  -------
          ! |   2 3 |
          ! | 1 2*3 |
          ! | 2 2 3 |
          !  -------
                   else

          !  ---------
          ! | 1/2 2 3 |
          ! | 1   2*3 |
          ! | 2   2 3 |
          !  ---------
                      if(bc_procedure_extra_checks) then
                         if(.not.(
     $                        ((grdpts_id(i-1,j+1).eq.interior_pt).or.
     $                         (grdpts_id(i-1,j+1).eq.bc_interior_pt)).and.
     $                        (grdpts_id(i+1,j+1).eq.bc_pt))) then

                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)

                         end if
                      end if

                      procedure_type = E_edge_type
                      i_proc         = i
                      j_proc         = j

                   end if
          !  -------
          ! |   x   |
          ! | 1 2*3 |
          ! | 2 2 3 |
          !  -------
                else
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
          end if

        end subroutine get_bc_interior_pt_procedure_1_1


        subroutine get_bc_interior_pt_procedure_1_2(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror


          !  -------
          ! | 3     |
          ! | 2 2*  |
          ! | 2 1   |
          !  -------
          if(grdpts_id(i,j-1).eq.interior_pt) then

             if(bc_procedure_extra_checks) then
                if(grdpts_id(i-1,j+1).ne.bc_pt) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

          !  -------
          ! | 3 2   |
          ! | 2 2*  |
          ! | 2 1   |
          !  -------
             if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                !  --------
                ! | 3 2 1/2|
                ! | 2 2*1  |
                ! | 2 1 1  |
                !  --------
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.interior_pt).and.(
     $                  (grdpts_id(i+1,j+1).eq.interior_pt).or.
     $                  (grdpts_id(i+1,j+1).eq.bc_interior_pt)))) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if

                end if

                procedure_type = NW_edge_type
                i_proc         = i-1
                j_proc         = j

             else

                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i,j+1).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! | 3 3   |
          ! | 2 2*2 |
          ! | 2 1   |
          !  -------
                if(grdpts_id(i+1,j).eq.bc_interior_pt) then
                   
          !  -------
          ! | 3 3 2 |
          ! | 2 2*2 |
          ! | 2 1   |
          !  -------
                   if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then

          !  -------
          ! | 3 3 2 |
          ! | 2 2*2 |
          ! | 2 1 1 |
          !  -------
                      if(bc_procedure_extra_checks) then
                         if(grdpts_id(i+1,j-1).ne.interior_pt) then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = NW_edge_type
                      i_proc         = i
                      j_proc         = j

          !  --------
          ! | 3 3 3  |
          ! | 2 2*2  |
          ! | 2 1 1/2|
          !  --------                  
                   else

                      if(bc_procedure_extra_checks) then
                         if(.not.(
     $                        (grdpts_id(i+1,j-1).ne.interior_pt).or.
     $                        (grdpts_id(i+1,j-1).ne.bc_interior_pt)
     $                        ))then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = N_edge_type
                      i_proc         = i
                      j_proc         = j                      

                   end if

          !  -------
          ! | 3 3   |
          ! | 2 2*x |
          ! | 2 1   |
          !  -------  
                else

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if

             end if
                   
          !  -------
          ! | 1     |
          ! | 2 2*  |
          ! | 2 3   |
          !  -------
          else

             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i,j-1).eq.bc_pt).and.
     $               (grdpts_id(i-1,j+1).eq.interior_pt))) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

          !  -------
          ! | 1     |
          ! | 2 2*2 |
          ! | 2 3   |
          !  -------             
             if(grdpts_id(i+1,j).eq.bc_interior_pt) then

          !  -------
          ! | 1     |
          ! | 2 2*2 |
          ! | 2 3 2 |
          !  -------
                if(grdpts_id(i+1,j-1).eq.bc_interior_pt) then

                   print '(''unable to differentiate:'')'
                   print '(''SW_edge and SE_edge'')'
                   print '()'

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

          !  --------
          ! | 1 1 1/2|
          ! | 2 2*2  |
          ! | 2 3 3  |
          !  --------
                else

                   if(bc_procedure_extra_checks) then

                      if(.not.(
     $                     (grdpts_id(i,j+1).eq.interior_pt).and.(
     $                     (grdpts_id(i+1,j+1).eq.interior_pt).or.
     $                     (grdpts_id(i+1,j+1).eq.bc_interior_pt))
     $                     )) then
                         
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)

                      end if

                   end if

                   procedure_type = SE_edge_type
                   i_proc         = i-1
                   j_proc         = j-1
                   
                end if

             else

          !  ---------
          ! | 1 2 2/3|
          ! | 2 2*3  |
          ! | 2 3 3  |
          !  --------
                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i+1,j).ne.bc_pt) then
                      
                      call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)

                   end if
                end if

                if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                   if(bc_procedure_extra_checks) then
                      if(.not.(
     $                     (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $                     (grdpts_id(i+1,j+1).eq.bc_pt))) then
                      
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)

                      end if
                   end if

                   procedure_type = SE_corner_type
                   i_proc         = i
                   j_proc         = j-1

                else
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if

             end if
             
          end if

        end subroutine get_bc_interior_pt_procedure_1_2


        subroutine get_bc_interior_pt_procedure_2_1(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc 
          logical                , intent(inout) :: ierror

          !  -------
          ! |       |
          ! |   2*3 |
          ! | 1 2 2 |
          !  -------
          if(grdpts_id(i+1,j).eq.bc_pt) then

             if(bc_procedure_extra_checks) then
                if(grdpts_id(i-1,j-1).ne.interior_pt) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

          !  ---------
          ! | 2/3 3 3 |
          ! | 2   2*3 |
          ! | 1   2 2 |
          !  ---------
             if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                if(bc_procedure_extra_checks) then

                   if(.not.(
     $                  ((grdpts_id(i-1,j+1).eq.bc_interior_pt).or.
     $                   (grdpts_id(i-1,j+1).eq.bc_pt)).and.
     $                   (grdpts_id(i,j+1).eq.bc_pt).and.
     $                   (grdpts_id(i+1,j+1).eq.bc_pt))) then

                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if

                end if

                procedure_type = NE_corner_type
                i_proc         = i
                j_proc         = j

             else
          !  -------
          ! |   2   |
          ! | 1 2*3 |
          ! | 1 2 2 |
          !  -------
                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i-1,j).ne.interior_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

                if(grdpts_id(i,j+1).eq.bc_interior_pt) then

          !  -------
          ! |   2 2 |
          ! | 1 2*3 |
          ! | 1 2 2 |
          !  -------                   
                   if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then

                      print '(''unable to differentiate:'')'
                      print '(''SE_edge and NE_edge'')'

                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

          !  ---------
          ! | 1/2 2 3 |
          ! | 1   2*3 |
          ! | 1   2 2 |
          !  ---------
                   else

                      if(bc_procedure_extra_checks) then
                         if(.not.(
     $                        ((grdpts_id(i-1,j+1).eq.interior_pt).or.
     $                         (grdpts_id(i-1,j+1).eq.bc_interior_pt)).and.
     $                        (grdpts_id(i+1,j+1).eq.bc_pt))) then
     $                        
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if                      

                      procedure_type = NE_edge_type
                      i_proc         = i
                      j_proc         = j-1

                   end if

          !  -------
          ! |   x 3 |
          ! | x 2*3 |
          ! | x 2 2 |
          !  -------
                else
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                   
                end if
             end if

          !  -------
          ! |       |
          ! |   2*1 |
          ! | 3 2 2 |
          !  -------
          else

             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j).eq.interior_pt)
     $               ))then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
             
          !  ---------
          ! | 2/1 1 1 |
          ! | 2   2*1 |
          ! | 3   2 2 |
          !  ---------
             if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                if(bc_procedure_extra_checks) then
                   
                   if(.not.(
     $                  ((grdpts_id(i-1,j+1).eq.interior_pt).or.
     $                   (grdpts_id(i-1,j+1).eq.bc_interior_pt)).and.
     $                  (grdpts_id(i,j+1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j+1).eq.interior_pt))) then

                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if

                end if

                procedure_type = SW_edge_type
                i_proc         = i-1
                j_proc         = j-1

             else

          !  -------
          ! |   2   |
          ! | 3 2*1 |
          ! | 3 2 2 |
          !  -------
                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i-1,j).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

                if(grdpts_id(i,j+1).eq.bc_interior_pt) then

          !  -------
          ! | 2 2 1 |
          ! | 3 2*1 |
          ! | 3 2 2 |
          !  -------
                   if(grdpts_id(i-1,j+1).eq.bc_interior_pt) then

                      if(bc_procedure_extra_checks) then
                         if(grdpts_id(i+1,j+1).ne.interior_pt) then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = SW_edge_type
                      i_proc         = i-1
                      j_proc         = j

          !  --------
          ! | 3 2 2/1|
          ! | 3 2*1  |
          ! | 3 2 2  |
          !  --------
                   else

                      if(bc_procedure_extra_checks) then
                         if(.not.(
     $                        (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                        ((grdpts_id(i+1,j+1).eq.interior_pt).or.
     $                         (grdpts_id(i+1,j+1).eq.bc_interior_pt))))then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = W_edge_type
                      i_proc         = i
                      j_proc         = j

                   end if

          !  -------
          ! |   x   |
          ! | 3 2*1 |
          ! | 3 2 2 |
          !  -------
                else
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if

             end if
          end if

        end subroutine get_bc_interior_pt_procedure_2_1


        subroutine get_bc_interior_pt_procedure_2_2(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror

          !  -------
          ! |   3 3 |
          ! | 2 2*3 |
          ! | 1 2 x |
          !  -------
          if(grdpts_id(i-1,j-1).eq.interior_pt) then

             !  -------
             ! | 2/3 3 3 |
             ! | 2   2*3 |
             ! | 1   2 3 |
             !  ---------
             if(bc_procedure_extra_checks) then
                
                if(.not.(
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j).eq.bc_pt).and.(
     $               (grdpts_id(i-1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i-1,j+1).eq.bc_pt)).and.
     $               (grdpts_id(i+1,j+1).eq.bc_pt))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if

             end if

             procedure_type = NE_corner_type
             i_proc         = i
             j_proc         = j

          !  -------
          ! |   1 1 |
          ! | 2 2*1 |
          ! | 3 2 x |
          !  -------
          else

             !  ---------
             ! | 2/1 1 1 |
             ! | 2   2*1 |
             ! | 3   2 1 |
             !  ---------
             if(bc_procedure_extra_checks) then
                
                if(.not.(
     $               (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j).eq.interior_pt).and.(
     $               (grdpts_id(i-1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i-1,j+1).eq.interior_pt)).and.
     $               (grdpts_id(i+1,j+1).eq.interior_pt))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if

             end if

             procedure_type = SW_edge_type
             i_proc         = i-1
             j_proc         = j-1
             
          end if

        end subroutine get_bc_interior_pt_procedure_2_2


        subroutine get_bc_interior_pt_procedure_2_3(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror
 
          !  -------
          ! |       |
          ! | x 2*2 |
          ! | x 2 3 |
          !  -------
          if(grdpts_id(i+1,j-1).eq.bc_pt) then

             !  --------
             ! | 1 1 2/1|
             ! | 1 2*2  |
             ! | 1 2 3  |
             !  --------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i-1,j).eq.interior_pt).and.
     $               (grdpts_id(i-1,j+1).eq.interior_pt).and.
     $               (grdpts_id(i,j+1).eq.interior_pt).and.(
     $               (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i+1,j+1).eq.interior_pt)))) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

             procedure_type = SE_edge_type
             i_proc         = i
             j_proc         = j-1

          !  -------
          ! |       |
          ! | x 2*2 |
          ! | x 2 1 |
          !  -------
          else

             !  --------
             ! | 3 3 2/3|
             ! | 3 2*2  |
             ! | 3 2 1  |
             !  --------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i-1,j).eq.bc_pt).and.
     $               (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $               (grdpts_id(i,j+1).eq.bc_pt).and.(
     $               (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i+1,j+1).eq.bc_pt))
     $               )) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

             procedure_type = NW_corner_type
             i_proc         = i-1
             j_proc         = j

          end if

        end subroutine get_bc_interior_pt_procedure_2_3


        subroutine get_bc_interior_pt_procedure_2_4(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror
 
          !  -------
          ! | 2 2 3 |
          ! | x 2*3 |
          ! | x 2 3 |
          !  -------
          if(grdpts_id(i+1,j).eq.bc_pt) then

          !  -------
          ! | 2 2 3 |
          ! | 1 2*3 |
          ! | 1 2 3 |
          !  -------
             if(bc_procedure_extra_checks) then

                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i-1,j).eq.interior_pt).and.
     $               (grdpts_id(i+1,j+1).eq.bc_pt))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

             procedure_type = E_edge_type
             i_proc         = i
             j_proc         = j

          !  -------
          ! | 2 2 1 |
          ! | 3 2*1 |
          ! | 3 2 1 |
          !  -------
          else

             if(bc_procedure_extra_checks) then

                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i-1,j).eq.bc_pt).and.
     $               (grdpts_id(i+1,j+1).eq.interior_pt))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

             procedure_type = SW_edge_type
             i_proc         = i-1
             j_proc         = j

          end if

        end subroutine get_bc_interior_pt_procedure_2_4


        subroutine get_bc_interior_pt_procedure_2_5(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror
 
          !  -------
          ! | x 2 2 |
          ! | x 2*x |
          ! | x 2 x |
          !  -------
          if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then

          !  -------
          ! | x 2 2 |
          ! | 3 2*x |
          ! | x 2 x |
          !  -------
             if(grdpts_id(i-1,j).eq.bc_pt) then

                !  -------
                ! | 3 2 2 |
                ! | 3 2*1 |
                ! | 3 2 1 |
                !  -------
                if(bc_procedure_extra_checks) then

                   if(.not.(
     $                  (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.interior_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.bc_pt))) then
                      
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if
                end if

                procedure_type = W_edge_type
                i_proc         = i
                j_proc         = j

          !  -------
          ! | x 2 2 |
          ! | 1 2*x |
          ! | x 2 x |
          !  -------
             else

                !  -------
                ! | 1 2 2 |
                ! | 1 2*3 |
                ! | 1 2 3 |
                !  -------
                if(bc_procedure_extra_checks) then

                   if(.not.(
     $                  (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j).eq.bc_pt).and.
     $                  (grdpts_id(i-1,j).eq.interior_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.interior_pt))) then
                      
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if
                end if

                procedure_type = SE_edge_type
                i_proc         = i
                j_proc         = j

             end if

          !  -------
          ! | x 2 x |
          ! | x 2*x |
          ! | x 2 x |
          !  -------
          else

          !  -------
          ! | 3 2 x |
          ! | 3 2*x |
          ! | 3 2 x |
          !  -------
             if(grdpts_id(i-1,j).eq.bc_pt) then

                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.interior_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                  (grdpts_id(i+1,j+1).eq.interior_pt))) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

                procedure_type = W_edge_type
                i_proc         = i
                j_proc         = j

          !  -------
          ! | 1 2 3 |
          ! | 1 2*3 |
          ! | 1 2 3 |
          !  -------
             else
                
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                  (grdpts_id(i-1,j).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j).eq.bc_pt).and.
     $                  (grdpts_id(i-1,j+1).eq.interior_pt).and.
     $                  (grdpts_id(i+1,j+1).eq.bc_pt))) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

                procedure_type = E_edge_type
                i_proc         = i
                j_proc         = j

             end if

          end if

        end subroutine get_bc_interior_pt_procedure_2_5


        subroutine get_bc_interior_pt_procedure_3_1(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror
 
          !  -------
          ! |     x |
          ! |   2*2 |
          ! | x 3 2 |
          !  -------
          if(grdpts_id(i,j-1).eq.bc_pt) then

          !  -------
          ! |     1 |
          ! |   2*2 |
          ! | 3 3 2 |
          !  -------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j+1).eq.interior_pt))) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
             
          !  -------
          ! |   2 1 |
          ! |   2*2 |
          ! | 3 3 2 |
          !  -------             
             if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                !  ---------
                ! | 3/2 2 1 |
                ! | 3   2*2 |
                ! | 3   3 2 |
                !  ---------
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i-1,j).eq.bc_pt).and.
     $                  ((grdpts_id(i-1,j+1).eq.bc_pt).or.
     $                   (grdpts_id(i-1,j+1).eq.bc_interior_pt)))) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

                procedure_type = SW_corner_type
                i_proc         = i-1
                j_proc         = j-1

             else
          !  -------
          ! |   x 1 |
          ! | 2 2*2 |
          ! | 3 3 2 |
          !  -------             
                if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                   !  ---------
                   ! | 2/1 1 1 |
                   ! | 2   2*2 |
                   ! | 3   3 2 |
                   !  ---------
                   if(bc_procedure_extra_checks) then

                      if(.not.(
     $                     ((grdpts_id(i-1,j+1).eq.bc_interior_pt).or.
     $                     (grdpts_id(i-1,j+1).eq.interior_pt)).and.
     $                     (grdpts_id(i,j+1).eq.interior_pt))) then
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)
                      end if
                      
                   end if

                   procedure_type = SW_edge_type
                   i_proc         = i
                   j_proc         = j-1

                else
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                   
                end if
             end if
          
          !  -------
          ! |     3 |
          ! |   2*2 |
          ! | 1 1 2 |
          !  -------
          else

             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j+1).eq.bc_pt))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

          !  -------
          ! |   2 3 |
          ! |   2*2 |
          ! | 1 1 2 |
          !  -------
             if(grdpts_id(i,j+1).eq.bc_interior_pt) then

                !  ---------
                ! | 1/2 2 3 |
                ! | 1   2*2 |
                ! | 1   1 2 |
                !  ---------
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i-1,j).eq.interior_pt).and.(
     $                  (grdpts_id(i-1,j+1).eq.interior_pt).or.
     $                  (grdpts_id(i-1,j+1).eq.bc_interior_pt)))) then

                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                      
                   end if
                end if

                procedure_type = NE_edge_type
                i_proc         = i
                j_proc         = j

          !  -------
          ! |   3 3 |
          ! |   2*2 |
          ! | 1 1 2 |
          !  -------
             else

                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i,j+1).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! | 2 3 3 |
          ! |   2*2 |
          ! | 1 1 2 |
          !  -------
                if(grdpts_id(i-1,j+1).eq.bc_interior_pt) then

          !  -------
          ! | 2 3 3 |
          ! | 2 2*2 |
          ! | 1 1 2 |
          !  -------
                   if(bc_procedure_extra_checks) then
                      if(grdpts_id(i-1,j).ne.bc_interior_pt) then
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)
                      end if
                   end if

                   procedure_type = NE_edge_type
                   i_proc         = i-1
                   j_proc         = j

          !  -------
          ! | x 3 3 |
          ! | 2 2*2 |
          ! | 1 1 2 |
          !  -------
                else

                   if(grdpts_id(i-1,j).eq.bc_interior_pt) then

                      if(bc_procedure_extra_checks) then
                         if(grdpts_id(i-1,j+1).ne.bc_pt) then
                            call error_bc_interior_pt_procedure(
     $                           i,j,grdpts_id,ierror)
                         end if
                      end if

                      procedure_type = N_edge_type
                      i_proc         = i
                      j_proc         = j

          !  -------
          ! | x x 3 |
          ! | x 2*2 |
          ! | x x 2 |
          !  -------
                   else
                      
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)

                   end if
                end if
             end if
          end if

        end subroutine get_bc_interior_pt_procedure_3_1


        subroutine get_bc_interior_pt_procedure_4_1(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror

          
          !  -------
          ! |       |
          ! | 2 2*2 |
          ! |   3   |
          !  -------
          if(grdpts_id(i,j-1).eq.bc_pt) then

          !  -----------
          ! | 2/1 1 2/1 |
          ! | 2   2*2   |
          ! | 3   3 3   |
          !  -----------
             if(bc_procedure_extra_checks) then

                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.(
     $               (grdpts_id(i-1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i-1,j+1).eq.interior_pt)).and.
     $               (grdpts_id(i,j+1).eq.interior_pt).and.(
     $               (grdpts_id(i+1,j+1).eq.bc_interior_pt).or.
     $               (grdpts_id(i+1,j+1).eq.interior_pt))
     $               )) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if

             end if

             procedure_type = S_edge_type
             i_proc         = i
             j_proc         = j

          !  -------
          ! |       |
          ! | 2 2*2 |
          ! |   1   |
          !  -------
          else
             
          !  -------
          ! |       |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j-1).eq.interior_pt)
     $               )) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
             
          !  -------
          ! | 2     |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
             if(grdpts_id(i-1,j+1).eq.bc_interior_pt) then
                
          !  -------
          ! | 2 3   |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                if(bc_procedure_extra_checks) then
                   if(grdpts_id(i,j+1).ne.bc_pt) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! | 2 3 2 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then
                   print '(''unable to differentiate:'')'
                   print '(''NE_edge or NW_edge'')'
                   print '()'
                   
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

          !  -------
          ! | 2 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                else

                   if(bc_procedure_extra_checks) then
                      if(grdpts_id(i+1,j+1).ne.bc_pt) then
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)
                      end if
                   end if
                
                   procedure_type = NE_edge_type
                   i_proc         = i-1
                   j_proc         = j

                end if

             else
          !  -------
          ! | 3 3   |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                if(bc_procedure_extra_checks) then
                   if(.not.(
     $                  (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                  (grdpts_id(i,j+1).eq.bc_pt))) then
                      call error_bc_interior_pt_procedure(
     $                     i,j,grdpts_id,ierror)
                   end if
                end if

          !  -------
          ! | 3 3 2 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                if(grdpts_id(i+1,j+1).eq.bc_interior_pt) then

                   procedure_type = NW_edge_type
                   i_proc         = i
                   j_proc         = j

          !  ------- 
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
                else

                   if(bc_procedure_extra_checks) then
                      if(grdpts_id(i+1,j+1).ne.bc_pt) then
                         call error_bc_interior_pt_procedure(
     $                        i,j,grdpts_id,ierror)
                      end if
                   end if

                   procedure_type = N_edge_type
                   i_proc         = i
                   j_proc         = j

                end if
             end if
          end if

        end subroutine get_bc_interior_pt_procedure_4_1


        subroutine get_bc_interior_pt_procedure_4_2(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror

          !  -------
          ! | x 2   |
          ! | 2 2*3 |
          ! | x x x |
          !  -------
          if(grdpts_id(i+1,j).eq.bc_pt) then

          !  --------
          ! | 1 2 3/2|
          ! | 2 2*3  |
          ! | 3 3 3  |
          !  --------  
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i-1,j+1).eq.interior_pt).and.(
     $               (grdpts_id(i+1,j+1).eq.bc_pt).or.
     $               (grdpts_id(i+1,j+1).eq.bc_interior_pt)))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

             procedure_type = SE_corner_type
             i_proc         = i
             j_proc         = j-1

          !  -------
          ! | x 2   |
          ! | 2 2*x |
          ! | x x x |
          !  -------
          else

          !  --------
          ! | 3 2 1/2|
          ! | 2 2*1  |
          ! | 1 1 1  |
          !  -------- 
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i-1,j+1).eq.bc_pt).and.(
     $               (grdpts_id(i+1,j+1).eq.interior_pt).or.
     $               (grdpts_id(i+1,j+1).eq.bc_interior_pt)))) then

                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)

                end if
             end if

             procedure_type = NW_edge_type
             i_proc         = i-1
             j_proc         = j

          end if

        end subroutine get_bc_interior_pt_procedure_4_2


        subroutine get_bc_interior_pt_procedure_5_1(
     $     i,j,grdpts_id,
     $     procedure_type,i_proc,j_proc,
     $     ierror)

          implicit none
          
          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer                , intent(out)   :: procedure_type
          integer                , intent(out)   :: i_proc
          integer                , intent(out)   :: j_proc
          logical                , intent(inout) :: ierror

          !  -------
          ! |   2 3 |
          ! | x 2*2 |
          ! | x x x |
          !  -------
          if(grdpts_id(i+1,j+1).eq.bc_pt) then

          !  ---------
          ! | 1/2 2 3 |
          ! | 1   2*2 |
          ! | 1   1 1 |
          !  ---------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i,j-1).eq.interior_pt).and.
     $               (grdpts_id(i+1,j-1).eq.interior_pt).and.
     $               (grdpts_id(i-1,j).eq.interior_pt).and.(
     $               (grdpts_id(i-1,j+1).eq.interior_pt).or.
     $               (grdpts_id(i-1,j+1).eq.bc_interior_pt)))) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if
             
             procedure_type = NE_edge_type
             i_proc         = i
             j_proc         = j

          !  -------
          ! |   2 3 |
          ! | x 2*2 |
          ! | x x x |
          !  -------
          else

          !  ---------
          ! | 3/2 2 1 |
          ! | 3   2*2 |
          ! | 3   3 3 |
          !  ---------
             if(bc_procedure_extra_checks) then
                if(.not.(
     $               (grdpts_id(i-1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i,j-1).eq.bc_pt).and.
     $               (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $               (grdpts_id(i-1,j).eq.bc_pt).and.(
     $               (grdpts_id(i-1,j+1).eq.bc_pt).or.
     $               (grdpts_id(i-1,j+1).eq.bc_interior_pt)))) then
                   call error_bc_interior_pt_procedure(
     $                  i,j,grdpts_id,ierror)
                end if
             end if

             procedure_type = SW_corner_type
             i_proc         = i-1
             j_proc         = j-1
             
          end if

        end subroutine get_bc_interior_pt_procedure_5_1


        !error pattern not found for bc_interior_pt
        subroutine error_bc_interior_pt_procedure(i,j,grdpts_id,ierror)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          logical                , intent(out) :: ierror

          print '(''bf_bc_procedure_module'')'
          print '(''get_bc_interior_pt_procedure'')'
          print '(''[i,j]: '', 2I4)', i,j
          print '(''grdpts(i-1:i+1,j-1:j+1)'')'
          print '(3I2)', grdpts_id(i-1:i+1,j+1)
          print '(3I2)', grdpts_id(i-1:i+1,j)
          print '(3I2)', grdpts_id(i-1:i+1,j-1)
          print '(''****************************************'')'
          print '()'
          ierror = .not.BF_SUCCESS

        end subroutine error_bc_interior_pt_procedure

      end module bf_layer_bc_procedure_module
