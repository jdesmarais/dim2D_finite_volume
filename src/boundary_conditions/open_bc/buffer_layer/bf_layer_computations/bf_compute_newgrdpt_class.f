      !> @file
      !> extension of bf_compute to compute the new grid points
      !> belonging to the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> extension of bf_compute to compute the new grid points
      !> belonging to the buffer layer
      !
      !> @date
      !> 03_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_compute_newgrdpt_class

        use bf_compute_class, only :
     $       bf_compute

        use bf_layer_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_bf_layer_data

        use bf_newgrdpt_class, only : 
     $       bf_newgrdpt

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_compute_newgrdpt


        !> @class bf_compute_newgrdpt
        !> class encapsulating the intrinsic procedures to compute the 
        !> new grid points from the data at the previous steps
        !
        !> @param get_data_for_newgrdpt
        !> extract the data needed for the computation of the new
        !> gridpoint (from the previous time step and the grdpts_id)
        !
        !> @param compute_newgrdpt
        !> attempt to compute the new grid point
        !---------------------------------------------------------------
        type, extends(bf_compute) :: bf_compute_newgrdpt

          contains

          procedure, pass :: get_data_for_newgrdpt
          procedure, pass :: compute_newgrdpt

        end type bf_compute_newgrdpt

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param dt
        !> time
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !
        !>@param integration_step_nopt
        !> function determining the procedure for the time integration
        !--------------------------------------------------------------
        subroutine get_data_for_newgrdpt(
     $       this,
     $       tmp_grdpts_id0,
     $       tmp_nodes0,
     $       gen_coords)

          implicit none

          class(bf_compute_newgrdpt)       , intent(in)    :: this
          integer        , dimension(:,:)  , intent(inout) :: tmp_grdpts_id0
          real(rkind)    , dimension(:,:,:), intent(inout) :: tmp_nodes0
          integer(ikind) , dimension(2,2)  , intent(in)    :: gen_coords

          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j
          integer        :: k


          if(allocated(this%alignment_tmp)) then

             !get the synchronization indices
             call get_indices_to_extract_bf_layer_data(
     $            this%alignment_tmp,
     $            gen_coords,
     $            size_x,
     $            size_y,
     $            i_recv,
     $            j_recv,
     $            i_send,
     $            j_send)
             
             !synchronize the grdpts_id
             do j=1, size_y
                do i=1, size_x
             
                   tmp_grdpts_id0(i_recv+i-1,j_recv+j-1) =
     $                  this%grdpts_id_tmp(i_send+i-1,j_send+j-1)
                   
                end do
             end do             
             
             !synchronize the nodes
             do k=1,ne
                do j=1, size_y
                   do i=1, size_x
             
                      tmp_nodes0(i_recv+i-1,j_recv+j-1,k) =
     $                     this%nodes_tmp(i_send+i-1,j_send+j-1,k)
             
                   end do
                end do
             end do

          end if

        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid points
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param dt
        !> time
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !
        !>@param integration_step_nopt
        !> function determining the procedure for the time integration
        !--------------------------------------------------------------
        function compute_newgrdpt(
     $     this,
     $     p_model, t, dt,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1,ierror)
     $     result(new_grdpt)

          implicit none

          class(bf_compute_newgrdpt)      , intent(in)  :: this
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind)                     , intent(in)  :: dt
          integer(ikind), dimension(2,2)  , intent(in)  :: bf_align1
          real(rkind)   , dimension(:)    , intent(in)  :: bf_x_map1
          real(rkind)   , dimension(:)    , intent(in)  :: bf_y_map1
          real(rkind)   , dimension(:,:,:), intent(in)  :: bf_nodes1
          integer(ikind)                  , intent(in)  :: i1
          integer(ikind)                  , intent(in)  :: j1
          logical                         , intent(out) :: ierror
          real(rkind)   , dimension(ne)                 :: new_grdpt

          integer(ikind)        :: i0,j0
          type(bf_newgrdpt)     :: bf_newgrdpt_used
          integer               :: nb_procedures
          integer, dimension(4) :: procedure_type
          integer, dimension(4) :: gradient_type

          
          print '(''bf_compute_newgrdpt_class'')'
          print '(''compute_newgrdpt'')'
          stop 'TO BE VALIDATED'

          ! get the indices of the new gridpoint for the
          ! grdpts_id_tmp
          i0 = bf_align1(1,1) - this%alignment_tmp(1,1) + i1
          j0 = bf_align1(2,1) - this%alignment_tmp(2,1) + j1


          ! compute the procedure_type and the gradient_type
          ! identifying the procedure that needs to be applied
          ! to compute the new grid point
          ierror = get_newgrdpt_procedure(
     $         i0,j0,
     $         this%grdpts_id_tmp,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type)

          if(ierror.neqv.BF_SUCCESS) then
             print '(''bf_compute_class'')'
             print '(''compute_newgrdpt'')'
             print '(''in get_newgrdpt_procedure()'')'
             print '(''configuration around newgrdpt not recognized'')'
             stop ''
          end if


          ! verify that there are enough grid points around the
          ! new grid point for its computation
          ierror = verify_data_for_newgrdpt(
     $         i0,j0,this%grdpts_id_tmp,
     $         procedure_type,
     $         gradient_type)

          if(ierror.eqv.BF_SUCCESS) then

             ! compute the new grid point
             new_grdpt = bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $         
     $            this%alignment_tmp,
     $            this%x_map_tmp,
     $            this%y_map_tmp,
     $            this%nodes_tmp,
     $            
     $            bf_align1,
     $            bf_x_map1,
     $            bf_y_map1,
     $            bf_nodes1,
     $         
     $            i1,j1,
     $            
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type)

          end if

        end function compute_newgrdpt

      end module bf_compute_newgrdpt_class
