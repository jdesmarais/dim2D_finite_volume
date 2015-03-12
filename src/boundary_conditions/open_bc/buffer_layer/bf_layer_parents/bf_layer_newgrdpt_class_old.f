      !> @file
      !> module encapsulating the bf_layer_time object.
      !> bf_layer_dyn enhanced with time integration functions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_dyn enhanced with time integration functions
      !
      !> @date
      ! 23_02_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_newgrdpt_class

        use bf_layer_time_class, only :
     $       bf_layer_time

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_bf_layer_data

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        private
        public :: bf_layer_newgrdpt


        !> @class bf_layer_newgrdpt_class
        !> bf_layer_time enhanced with new grdpt computation functions
        !
        !> @param get_data_for_newgrdpt
        !> extract the data needed to compute the new grid point
        !
        !> @param compute_newgrdpt
        !> computation of the new grid point
        !
        !> @param get_grdpts_id_part
        !> extract the grdpts_id
        !
        !> @param set_grdpts_id_part
        !> set the grpdts_id
        !-------------------------------------------------------------
        type, extends(bf_layer_time) :: bf_layer_newgrdpt

           contains           

           !procedures for updating the gridpts_id
           procedure,   pass :: get_grdpts_id_part
           procedure,   pass :: set_grdpts_id_part

           !procedures for the computation of new grid points
           procedure,   pass :: get_data_for_newgrdpt
           procedure,   pass :: compute_newgrdpt

        end type bf_layer_newgrdpt

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id matching the gen_coords
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param gen_coords
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        subroutine get_grdpts_id_part(
     $     this,
     $     tmp_grdptsid,
     $     gen_coords,
     $     previous_step)

          implicit none

          class(bf_layer_newgrdpt)            , intent(in)    :: this
          integer             , dimension(:,:), intent(inout) :: tmp_grdptsid
          integer(ikind)      , dimension(2,2), intent(in)    :: gen_coords
          logical             , optional      , intent(in)    :: previous_step

          
          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j

          logical :: previous_step_op

          if(present(previous_step)) then
             previous_step_op = previous_step
          else
             previous_step_op = .false.
          end if


          if(previous_step_op) then

            !get the grdpts_id at the previous time step
             call this%bf_compute_used%get_grdpts_id_part(
     $            tmp_grdptsid,
     $            gen_coords)

          else

             !get the synchronization indices
             call get_indices_to_extract_bf_layer_data(
     $            this%alignment,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)
             
             
             !fill the grid points asked
             do j=1, size_y
                do i=1, size_x
                   
                   tmp_grdptsid(i_recv+i-1,j_recv+j-1) =
     $                  this%grdpts_id(i_send+i-1,j_send+j-1)
             
                end do
             end do

          end if

        end subroutine get_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the grdpts_id matching the gen_coords
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param gen_coords
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        subroutine set_grdpts_id_part(
     $     this,
     $     tmp_grdptsid,
     $     gen_coords)

          implicit none

          class(bf_layer_newgrdpt)            , intent(inout) :: this
          integer             , dimension(:,:), intent(in)    :: tmp_grdptsid
          integer(ikind)      , dimension(2,2), intent(in)    :: gen_coords

          
          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j


          !get the synchronization indices
          call get_indices_to_extract_bf_layer_data(
     $         this%alignment,
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)


          !fill the grid points asked
          do j=1, size_y
             do i=1, size_x
                
                this%grdpts_id(i_send+i-1,j_send+j-1) = 
     $               tmp_grdptsid(i_recv+i-1,j_recv+j-1)

             end do
          end do

        end subroutine set_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id, the coordinate maps and the
        !> nodes at t-dt and t corresponding to the general
        !> coordinates gen_coords
        !    ___________________
        !   |                  _|_________
        !   |    buffer layer |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
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
        subroutine get_data_for_newgrdpt(
     $     this,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_coords)

          implicit none

          class(bf_layer_newgrdpt)        , intent(in)    :: this
          integer       , dimension(:,:)  , intent(inout) :: tmp_grdpts_id0
          real(rkind)   , dimension(:,:,:), intent(inout) :: tmp_nodes0
          real(rkind)   , dimension(:,:,:), intent(inout) :: tmp_nodes1
          integer(ikind), dimension(2,2)  , intent(in)    :: gen_coords


          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j
          integer        :: k


          if(allocated(this%nodes)) then

             !get the synchronization indices
             !for nodes at t
             call get_indices_to_extract_bf_layer_data(
     $            this%alignment,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)        


             !extract nodes at t
             do k=1,ne
                do j=1, size_y
                   do i=1, size_x

                      tmp_nodes1(i_recv+i-1,j_recv+j-1,k) =
     $                     this%nodes(i_send+i-1,j_send+j-1,k)

                   end do
                end do
             end do


             !extract nodes and grdpts_id at t-dt
             call this%bf_compute_used%get_data_for_newgrdpt(
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            gen_coords)

          else

             print '(''bf_layer_time_class'')'
             print '(''get_data_for_newgrdpt'')'
             print '(''nodes is not allocated'')'
             stop ''

          end if

        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i1
        !> index identifying the x-coordinate of the new grid point
        !
        !>@param j1
        !> index identifying the y-coordinate of the new grid point
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        subroutine compute_newgrdpt(this,p_model,t,dt,i1,j1,ierror)

          implicit none

          class(bf_layer_newgrdpt), intent(inout) :: this
          type(pmodel_eq)         , intent(in)    :: p_model
          real(rkind)             , intent(in)    :: t
          real(rkind)             , intent(in)    :: dt
          integer(ikind)          , intent(in)    :: i1
          integer(ikind)          , intent(in)    :: j1
          logical                 , intent(out)   :: ierror

          this%nodes(i1,j1,:) = this%bf_compute_used%compute_newgrdpt(
     $         p_model, t, dt,
     $         this%alignment,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         i1,j1,ierror)

        end subroutine compute_newgrdpt

      end module bf_layer_newgrdpt_class
