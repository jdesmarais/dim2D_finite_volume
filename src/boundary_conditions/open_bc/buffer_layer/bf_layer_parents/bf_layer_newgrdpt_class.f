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

        use parameters_bf_layer, only :
     $       NEWGRDPT_NO_PREVIOUS_DATA_ERROR

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

c$$$           !procedures for updating the gridpts_id
c$$$           procedure,   pass :: get_grdpts_id_part
c$$$           procedure,   pass :: set_grdpts_id_part
c$$$
           !procedures for the computation of new grid points
c$$$           procedure,   pass :: get_data_for_newgrdpt
           procedure,   pass :: compute_newgrdpt

        end type bf_layer_newgrdpt

        contains


c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> get the grdpts_id matching the gen_coords
c$$$        !
c$$$        !> @date
c$$$        !> 27_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer object encapsulating the main
c$$$        !> tables extending the interior domain
c$$$        !
c$$$        !>@param gen_coords
c$$$        !> logical identifying whether the buffer layer should be
c$$$        !> removed or not
c$$$        !--------------------------------------------------------------
c$$$        subroutine get_grdpts_id_part(
c$$$     $     this,
c$$$     $     tmp_grdptsid,
c$$$     $     gen_coords,
c$$$     $     previous_step)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_layer_newgrdpt)            , intent(in)    :: this
c$$$          integer             , dimension(:,:), intent(inout) :: tmp_grdptsid
c$$$          integer(ikind)      , dimension(2,2), intent(in)    :: gen_coords
c$$$          logical             , optional      , intent(in)    :: previous_step
c$$$
c$$$          
c$$$          integer(ikind) :: size_x,size_y
c$$$          integer(ikind) :: i_recv,i_send,j_recv,j_send
c$$$          integer(ikind) :: i,j
c$$$
c$$$          logical :: previous_step_op
c$$$
c$$$          if(present(previous_step)) then
c$$$             previous_step_op = previous_step
c$$$          else
c$$$             previous_step_op = .false.
c$$$          end if
c$$$
c$$$
c$$$          if(previous_step_op) then
c$$$
c$$$            !get the grdpts_id at the previous time step
c$$$             call this%bf_compute_used%get_grdpts_id_part(
c$$$     $            tmp_grdptsid,
c$$$     $            gen_coords)
c$$$
c$$$          else
c$$$
c$$$             !get the synchronization indices
c$$$             call get_indices_to_extract_bf_layer_data(
c$$$     $            this%alignment,
c$$$     $            gen_coords,
c$$$     $            size_x, size_y,
c$$$     $            i_recv, j_recv,
c$$$     $            i_send, j_send)
c$$$             
c$$$             
c$$$             !fill the grid points asked
c$$$             do j=1, size_y
c$$$                do i=1, size_x
c$$$                   
c$$$                   tmp_grdptsid(i_recv+i-1,j_recv+j-1) =
c$$$     $                  this%grdpts_id(i_send+i-1,j_send+j-1)
c$$$             
c$$$                end do
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end subroutine get_grdpts_id_part
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> set the grdpts_id matching the gen_coords
c$$$        !
c$$$        !> @date
c$$$        !> 21_01_2015 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer object encapsulating the main
c$$$        !> tables extending the interior domain
c$$$        !
c$$$        !>@param gen_coords
c$$$        !> logical identifying whether the buffer layer should be
c$$$        !> removed or not
c$$$        !--------------------------------------------------------------
c$$$        subroutine set_grdpts_id_part(
c$$$     $     this,
c$$$     $     tmp_grdptsid,
c$$$     $     gen_coords)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_layer_newgrdpt)            , intent(inout) :: this
c$$$          integer             , dimension(:,:), intent(in)    :: tmp_grdptsid
c$$$          integer(ikind)      , dimension(2,2), intent(in)    :: gen_coords
c$$$
c$$$          
c$$$          integer(ikind) :: size_x,size_y
c$$$          integer(ikind) :: i_recv,i_send,j_recv,j_send
c$$$          integer(ikind) :: i,j
c$$$
c$$$
c$$$          !get the synchronization indices
c$$$          call get_indices_to_extract_bf_layer_data(
c$$$     $         this%alignment,
c$$$     $         gen_coords,
c$$$     $         size_x, size_y,
c$$$     $         i_recv, j_recv,
c$$$     $         i_send, j_send)
c$$$
c$$$
c$$$          !fill the grid points asked
c$$$          do j=1, size_y
c$$$             do i=1, size_x
c$$$                
c$$$                this%grdpts_id(i_send+i-1,j_send+j-1) = 
c$$$     $               tmp_grdptsid(i_recv+i-1,j_recv+j-1)
c$$$
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine set_grdpts_id_part
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> get the grdpts_id, the coordinate maps and the
c$$$        !> nodes at t-dt and t corresponding to the general
c$$$        !> coordinates gen_coords
c$$$        !    ___________________
c$$$        !   |                  _|_________
c$$$        !   |    buffer layer |/|         |
c$$$        !   |                 |/|  tmp    |
c$$$        !   !                 !/!         !
c$$$        !                   overlapping which is copied
c$$$        !                     from buffer layer to tmp
c$$$        !
c$$$        !> @date
c$$$        !> 18_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer object encapsulating the main
c$$$        !> tables extending the interior domain
c$$$        !
c$$$        !>@param tmp_grdpts_id0
c$$$        !> array with the grdpts_id data
c$$$        !
c$$$        !>@param tmp_nodes0
c$$$        !> array with the grid points data at t-dt
c$$$        !
c$$$        !>@param tmp_nodes1
c$$$        !> array with the grid points data at t
c$$$        !
c$$$        !>@param gen_coords
c$$$        !> coordinates of the SW corner and the NE corners of the
c$$$        !> tmp arrays computed
c$$$        !--------------------------------------------------------------
c$$$        subroutine get_data_for_newgrdpt(
c$$$     $     this,
c$$$     $     tmp_grdpts_id0,
c$$$     $     tmp_nodes0,
c$$$     $     tmp_nodes1,
c$$$     $     gen_coords)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_layer_newgrdpt)        , intent(in)    :: this
c$$$          integer       , dimension(:,:)  , intent(inout) :: tmp_grdpts_id0
c$$$          real(rkind)   , dimension(:,:,:), intent(inout) :: tmp_nodes0
c$$$          real(rkind)   , dimension(:,:,:), intent(inout) :: tmp_nodes1
c$$$          integer(ikind), dimension(2,2)  , intent(in)    :: gen_coords
c$$$
c$$$
c$$$          integer(ikind) :: size_x,size_y
c$$$          integer(ikind) :: i_recv,i_send,j_recv,j_send
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$
c$$$          if(allocated(this%nodes)) then
c$$$
c$$$             !get the synchronization indices
c$$$             !for nodes at t
c$$$             call get_indices_to_extract_bf_layer_data(
c$$$     $            this%alignment,
c$$$     $            gen_coords,
c$$$     $            size_x, size_y,
c$$$     $            i_recv, j_recv,
c$$$     $            i_send, j_send)        
c$$$
c$$$
c$$$             !extract nodes at t
c$$$             do k=1,ne
c$$$                do j=1, size_y
c$$$                   do i=1, size_x
c$$$
c$$$                      tmp_nodes1(i_recv+i-1,j_recv+j-1,k) =
c$$$     $                     this%nodes(i_send+i-1,j_send+j-1,k)
c$$$
c$$$                   end do
c$$$                end do
c$$$             end do
c$$$
c$$$
c$$$             !extract nodes and grdpts_id at t-dt
c$$$             call this%bf_compute_used%get_data_for_newgrdpt(
c$$$     $            tmp_grdpts_id0,
c$$$     $            tmp_nodes0,
c$$$     $            gen_coords)
c$$$
c$$$          else
c$$$
c$$$             print '(''bf_layer_time_class'')'
c$$$             print '(''get_data_for_newgrdpt'')'
c$$$             print '(''nodes is not allocated'')'
c$$$             stop ''
c$$$
c$$$          end if
c$$$
c$$$        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point if it is possible using the
        !> data available in the buffer layer or the data in the
        !> interior domain
        !
        !> @date
        !> 12_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param ngrdpt_coords
        !> indices identifying the coordinate of the new grid point
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        function compute_newgrdpt(
     $       this,
     $       p_model,
     $       t,
     $       dt,
     $       bf_newgrdpt_coords1,
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type,
     $       grdpts_available,
     $       data_needed_bounds0)
     $       result(ierror)

          implicit none

          class(bf_layer_newgrdpt)           , intent(inout) :: this
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          integer(ikind), dimension(2)       , intent(in)    :: bf_newgrdpt_coords1
          integer                            , intent(out)   :: nb_procedures
          integer       , dimension(4)       , intent(out)   :: procedure_type
          integer       , dimension(4)       , intent(out)   :: gradient_type
          logical       , dimension(4)       , intent(out)   :: grdpts_available
          integer(ikind), dimension(2,2)     , intent(out)   :: data_needed_bounds0
          integer                                            :: ierror


          logical :: bf_can_exchange_with_neighbor1
          logical :: bf_can_exchange_with_neighbor2


          if(this%does_previous_timestep_exist()) then

             bf_can_exchange_with_neighbor1 = this%can_exchange_with_neighbor1()
             bf_can_exchange_with_neighbor2 = this%can_exchange_with_neighbor2()

             ierror = this%bf_compute_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            this%localization,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            this%alignment,
     $            this%x_map,
     $            this%y_map,
     $            this%nodes,
     $            bf_newgrdpt_coords1,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type,
     $            grdpts_available,
     $            data_needed_bounds0)
             
          else

             ierror = NEWGRDPT_NO_PREVIOUS_DATA_ERROR

          end if

        end function compute_newgrdpt

      end module bf_layer_newgrdpt_class
