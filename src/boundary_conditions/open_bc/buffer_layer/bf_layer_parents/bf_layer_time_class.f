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
      module bf_layer_time_class

        use bf_layer_bc_sections_class, only :
     $       bf_layer_bc_sections

        use bf_layer_dyn_class, only :
     $       bf_layer_dyn

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_bf_layer_data

        use bc_operators_class, only :
     $       bc_operators

        use bf_compute_newgrdpt_class, only :
     $       bf_compute_newgrdpt

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_bf_layer, only :
     $       BF_SUCCESS,
     $       bc_interior_pt

        use parameters_input, only :
     $       bc_size,
     $       debug,
     $       nx,ny,ne
        
        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators
        

        private
        public :: bf_layer_time


        !> @class bf_layer_time
        !> bf_layer_dyn enhanced with time integration functions
        !
        !> @param bf_compute_used
        !> object encapsulating the time integration functions
        !
        !> @param x_borders
        !> x_borders for the time integration
        !
        !> @param y_borders
        !> y_borders for the time integration
        !
        !> @param bc_sections
        !> identification of the boundary sections in the buffer
        !> layer
        !
        !> @param get_data_for_newgrdpt
        !> extract the data needed to compute the new grid point
        !
        !> @param compute_newgrdpt
        !> computation of the new grid point
        !
        !> @param does_previous_timestep_exist
        !> check whether the previous time step is stored in bf_compute
        !
        !> @param get_grdpts_id_part
        !> extract the grdpts_id
        !
        !> @param set_grdpts_id_part
        !> set the grpdts_id
        !
        !> @param apply_initial_conditions
        !> apply the initial conditions in the buffer layer
        !
        !> @param allocate_before_timeInt
        !> allocate the temporary arrays to compute the time
        !> integration step
        !
        !> @param deallocate_after_timeInt
        !> deallocate the temporary arrays to compute the time
        !> integration step
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
        !
        !> @param get_time_dev
        !> get the time derivatives
        !
        !> @param set_x_borders
        !> set the index borders for time integration in the
        !> x-direction
        !
        !> @param set_y_borders
        !> set the index borders for time integration in the
        !> y-direction
        !
        !> @param get_x_borders
        !> get the index borders for time integration in the
        !> x-direction
        !
        !> @param get_y_borders
        !> get the index borders for time integration in the
        !> y-direction
        !
        !> @param set_N_bc_sections
        !> set the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param set_S_bc_sections
        !> set the borders of the south boundary layer computed
        !> by time integration
        !
        !> @param remove_N_bc_sections
        !> remove the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param remove_S_bc_sections
        !> remove the borders of the south boundary layer computed
        !> by time integration
        !
        !> @param get_N_bc_sections
        !> get the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param get_S_bc_sections
        !> get the borders of the south boundary layer computed
        !> by time integration
        !-------------------------------------------------------------
        type, extends(bf_layer_dyn) :: bf_layer_time

           type(bf_compute_newgrdpt)                   :: bf_compute_used
           integer(ikind), dimension(2)                :: x_borders
           integer(ikind), dimension(2)                :: y_borders
           integer       , dimension(:,:), allocatable :: bc_sections

           contains           

           !procedures for updating the gridpts_id
           procedure,   pass :: get_grdpts_id_part
           procedure,   pass :: set_grdpts_id_part

           !for time integration: interior + boundaries
           procedure,   pass :: apply_initial_conditions
           procedure,   pass :: update_bc_sections
           procedure,   pass :: allocate_before_timeInt
           procedure,   pass :: deallocate_after_timeInt
           procedure,   pass :: compute_time_dev
           procedure,   pass :: compute_integration_step

           !procedures for the computation of new grid points
           procedure,   pass :: does_previous_timestep_exist
           procedure,   pass :: get_data_for_newgrdpt
           procedure,   pass :: compute_newgrdpt

           !procedures for setting the integration borders
           !and the boundary sections
           procedure,   pass :: set_x_borders
           procedure,   pass :: set_y_borders
           procedure,   pass :: get_x_borders
           procedure,   pass :: get_y_borders           

           !for tests
           procedure,   pass :: get_time_dev !only for tests

        end type bf_layer_time

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

          class(bf_layer_time)                , intent(in)    :: this
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

          class(bf_layer_time)                , intent(inout) :: this
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
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine apply_initial_conditions(this,p_model)

          implicit none

          class(bf_layer_time), intent(inout) :: this
          type(pmodel_eq)     , intent(in)    :: p_model

          call p_model%apply_ic(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map)

        end subroutine apply_initial_conditions


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> identify the position of the boundary sections in
        !> the buffer layer
        !
        !> @date
        !> 05_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine update_bc_sections(this)

          implicit none

          class(bf_layer_time), intent(inout) :: this
          

          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used
          integer(ikind)             :: i,j
          logical                    :: ierror
          

          if(allocated(this%grdpts_id)) then
          
             !initialize the object gathering information
             !about the bc_sections identified in the
             !buffer layer
             call bf_layer_bc_sections_used%ini()

             !identify the boundary sections
             do j=2,size(this%grdpts_id,2)-1
                do i=2, size(this%grdpts_id,1)-1

                   if(this%grdpts_id(i,j).eq.bc_interior_pt) then

                      call bf_layer_bc_sections_used%analyse_grdpt(
     $                     i,j,
     $                     this%grdpts_id,
     $                     ierror)

                      if(ierror.neqv.BF_SUCCESS) then
                         print '(''bf_layer_time_class'')'
                         print '(''update_bc_sections'')'
                         print '(''failed identifying the bc_section'')'
                         print '(''at (i,j)=('',2I4,'')'')', i,j
                         stop ''
                      end if

                   end if
                end do
             end do

             !finalize the identification of the bc_sections
             if(allocated(this%bc_sections)) then
                deallocate(this%bc_sections)
             end if
             call bf_layer_bc_sections_used%finalize_bc_sections(
     $            this%x_borders,
     $            this%y_borders,
     $            this%bc_sections)             

          else
             print '(''bf_layer_time_class'')'
             print '(''update_bc_sections'')' 
             print '(''enable to update the bc_sections'')'
             print '(''as grdpts_id is not allocated'')'
             stop ''
          end if

        end subroutine update_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine allocate_before_timeInt(this)

          implicit none

          class(bf_layer_time), intent(inout) :: this

          if(allocated(this%nodes)) then

             call this%bf_compute_used%allocate_tables(
     $            size(this%nodes,1),
     $            size(this%nodes,2),
     $            this%alignment,
     $            this%x_map,
     $            this%y_map,
     $            this%grdpts_id)

          else

             print '(''bf_layer_time_class'')'
             print '(''allocate_before_timeInt'')'
             print '(''nodes is not allocated'')'
             stop ''

          end if

        end subroutine allocate_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine deallocate_after_timeInt(this)

          implicit none

          class(bf_layer_time), intent(inout) :: this

          call this%bf_compute_used%deallocate_tables()

        end subroutine deallocate_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used,
     $     interior_nodes)

          implicit none

          class(bf_layer_time)            , intent(inout) :: this
          type(td_operators)              , intent(in)    :: td_operators_used
          real(rkind)                     , intent(in)    :: t
          type(sd_operators)              , intent(in)    :: s
          type(pmodel_eq)                 , intent(in)    :: p_model
          type(bc_operators)              , intent(in)    :: bc_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes

          if(allocated(this%nodes)) then

             if(
     $            (this%x_borders(1).ge.1).and.
     $            (this%x_borders(2).le.size(this%nodes,1)).and.
     $            (this%y_borders(1).ge.1).and.
     $            (this%y_borders(2).le.size(this%nodes,2))) then
                
                call this%bf_compute_used%compute_time_dev(
     $               td_operators_used,
     $               t, this%nodes, this%x_map, this%y_map,
     $               s,p_model,bc_used,
     $               this%alignment,
     $               this%grdpts_id,
     $               interior_nodes,
     $               this%bc_sections,
     $               this%x_borders, this%y_borders)

             else
                
                print '(''bf_layer_time_class'')'
                print '(''compute_time_dev'')'
                print '(''the time integration borders are not set'')'
                stop ''

             end if

          else

             print '(''bf_layer_time_class'')'
             print '(''compute_time_dev'')'
             print '(''nodes is not allocated'')'
             stop ''

          end if

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param dt
        !> integration time step
        !
        !>@param integration_step_nopt
        !> procedure performing the time integration
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this,
     $     dt,
     $     integration_step_nopt)

          implicit none

          class(bf_layer_time)        , intent(inout) :: this
          real(rkind)                 , intent(in)    :: dt
          procedure(timeInt_step_nopt)                :: integration_step_nopt

          call this%bf_compute_used%compute_integration_step(
     $         this%grdpts_id,
     $         this%nodes,
     $         dt,
     $         this%x_borders,
     $         this%y_borders,
     $         integration_step_nopt)

        end subroutine compute_integration_step
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the previous time step is stored in the
        !> buffer layer
        !
        !> @date
        !> 20_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param exist
        !> says whether the previous time step is stored in the
        !> buffer layer
        !--------------------------------------------------------------
        function does_previous_timestep_exist(this) result(exist)

          implicit none

          class(bf_layer_time), intent(in) :: this
          logical                          :: exist

          exist = this%bf_compute_used%does_previous_timestep_exist()

        end function does_previous_timestep_exist


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

          class(bf_layer_time)            , intent(in)    :: this
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

          class(bf_layer_time), intent(inout) :: this
          type(pmodel_eq)     , intent(in)    :: p_model
          real(rkind)         , intent(in)    :: t
          real(rkind)         , intent(in)    :: dt
          integer(ikind)      , intent(in)    :: i1
          integer(ikind)      , intent(in)    :: j1
          logical             , intent(out)   :: ierror

          this%nodes(i1,j1,:) = this%bf_compute_used%compute_newgrdpt(
     $         p_model, t, dt,
     $         this%alignment,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         i1,j1,ierror)

        end subroutine compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the x-borders for the integration 
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_borders
        !> integration borders along the x-direction
        !--------------------------------------------------------------
        subroutine set_x_borders(this,x_borders)

          implicit none

          class(bf_layer_time)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: x_borders

          this%x_borders = x_borders

        end subroutine set_x_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the y-borders for the integration 
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_borders
        !> integration borders along the y-direction
        !--------------------------------------------------------------
        subroutine set_y_borders(this,y_borders)

          implicit none

          class(bf_layer_time)              , intent(inout) :: this
          integer(ikind)      , dimension(2), intent(in)    :: y_borders

          this%y_borders = y_borders

        end subroutine set_y_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the x-borders for the integration 
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return x_borders
        !> integration borders along the x-direction
        !--------------------------------------------------------------
        function get_x_borders(this) result(x_borders)

          implicit none

          class(bf_layer_time)        , intent(inout) :: this
          integer(ikind), dimension(2)                :: x_borders

          x_borders = this%x_borders

        end function get_x_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the y-borders for the integration 
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return y_borders
        !> integration borders along the y-direction
        !--------------------------------------------------------------
        function get_y_borders(this) result(y_borders)

          implicit none

          class(bf_layer_time)              , intent(inout) :: this
          integer(ikind)      , dimension(2)                :: y_borders

          y_borders = this%y_borders

        end function get_y_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the time derivatives
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_layer_time)                      , intent(in)  :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: time_dev

          call this%bf_compute_used%get_time_dev(time_dev)

        end subroutine get_time_dev


      end module bf_layer_time_class
