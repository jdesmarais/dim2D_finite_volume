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
     $       get_grdpts_id_from_bf_layer,
     $       get_nodes_from_bf_layer

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

           !procedures for the computation of new grid points
           procedure, pass :: compute_newgrdpt

           !procedure for the extraction of the grdpts_id or nodes
           procedure, pass :: extract_grdpts_id
           procedure, pass :: extract_nodes           

        end type bf_layer_newgrdpt

        contains


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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the grdpts_id at t-dt or t
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data extracted
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !
        !>@param previous_step
        !> logical identifying whether the grdpts_id at the previous step (t-dt)
        !> or the current grdpts_id (t) are extracted
        !--------------------------------------------------------------
        subroutine extract_grdpts_id(
     $     this,
     $     tmp_grdpts_id,
     $     gen_coords,
     $     extract_param_in,
     $     extract_param_out,
     $     previous_step)

          implicit none
          
          class(bf_layer_newgrdpt)                , intent(in)    :: this
          integer       , dimension(:,:)          , intent(inout) :: tmp_grdpts_id
          integer(ikind), dimension(2,2)          , intent(in)    :: gen_coords
          integer(ikind), dimension(6)  , optional, intent(in)    :: extract_param_in
          integer(ikind), dimension(6)  , optional, intent(out)   :: extract_param_out
          logical                       , optional, intent(in)    :: previous_step


          logical :: previous_step_op


          print '(''bf_layer_newgrdpt_class'')'
          print '(''extract_grdpts_id'')'
          stop 'NOT VALIDATED'


          if(present(previous_step)) then
             previous_step_op = previous_step
          else
             previous_step_op = .false.
          end if


          if(previous_step_op) then
             if(this%does_previous_timestep_exist()) then

                if(present(extract_param_in)) then
                   call this%bf_compute_used%extract_grdpts_id(
     $                  tmp_grdpts_id,
     $                  gen_coords,
     $                  extract_param_in=extract_param_in)

                else

                   if(present(extract_param_out)) then
                      call this%bf_compute_used%extract_grdpts_id(
     $                  tmp_grdpts_id,
     $                  gen_coords,
     $                  extract_param_out=extract_param_out)

                   else
                      call this%bf_compute_used%extract_grdpts_id(
     $                  tmp_grdpts_id,
     $                  gen_coords)

                   end if

                end if

             end if

          else

             if(present(extract_param_in)) then
                
                call get_grdpts_id_from_bf_layer(
     $               tmp_grdpts_id,
     $               gen_coords,
     $               this%alignment,
     $               this%grdpts_id,
     $               extract_param_in=extract_param_in)

             else

                if(present(extract_param_out)) then

                   call get_grdpts_id_from_bf_layer(
     $                  tmp_grdpts_id,
     $                  gen_coords,
     $                  this%alignment,
     $                  this%grdpts_id,
     $                  extract_param_out=extract_param_out)

                else

                   call get_grdpts_id_from_bf_layer(
     $                  tmp_grdpts_id,
     $                  gen_coords,
     $                  this%alignment,
     $                  this%grdpts_id)

                end if

             end if

          end if

        end subroutine extract_grdpts_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the nodes at t-dt or t
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_nodes
        !> array with the nodes data extracted
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !
        !>@param previous_step
        !> logical identifying whether the nodes at the previous step (t-dt)
        !> or the current nodes (t) are extracted
        !--------------------------------------------------------------
        subroutine extract_nodes(
     $     this,
     $     tmp_nodes,
     $     gen_coords,
     $     extract_param_in,
     $     extract_param_out,
     $     previous_step)

          implicit none
          
          class(bf_layer_newgrdpt)                , intent(in)    :: this
          real(rkind)   , dimension(:,:,:)        , intent(inout) :: tmp_nodes
          integer(ikind), dimension(2,2)          , intent(in)    :: gen_coords
          integer(ikind), dimension(6)  , optional, intent(in)    :: extract_param_in
          integer(ikind), dimension(6)  , optional, intent(out)   :: extract_param_out
          logical                       , optional, intent(in)    :: previous_step


          logical :: previous_step_op


          print '(''bf_layer_newgrdpt_class'')'
          print '(''extract_nodes'')'
          stop 'NOT VALIDATED'


          if(present(previous_step)) then
             previous_step_op = previous_step
          else
             previous_step_op = .false.
          end if


          if(previous_step_op) then
             if(this%does_previous_timestep_exist()) then

                if(present(extract_param_in)) then
                   call this%bf_compute_used%extract_nodes(
     $                  tmp_nodes,
     $                  gen_coords,
     $                  extract_param_in=extract_param_in)

                else

                   if(present(extract_param_out)) then
                      call this%bf_compute_used%extract_nodes(
     $                  tmp_nodes,
     $                  gen_coords,
     $                  extract_param_out=extract_param_out)

                   else
                      call this%bf_compute_used%extract_nodes(
     $                  tmp_nodes,
     $                  gen_coords)

                   end if

                end if

             end if

          else

             if(present(extract_param_in)) then
                
                call get_nodes_from_bf_layer(
     $               tmp_nodes,
     $               gen_coords,
     $               this%alignment,
     $               this%nodes,
     $               extract_param_in=extract_param_in)

             else

                if(present(extract_param_out)) then

                   call get_nodes_from_bf_layer(
     $                  tmp_nodes,
     $                  gen_coords,
     $                  this%alignment,
     $                  this%nodes,
     $                  extract_param_out=extract_param_out)

                else

                   call get_nodes_from_bf_layer(
     $                  tmp_nodes,
     $                  gen_coords,
     $                  this%alignment,
     $                  this%nodes)

                end if

             end if

          end if

        end subroutine extract_nodes      

      end module bf_layer_newgrdpt_class
