      !> @file
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @date
      !> 16_07_2014 - initial version         - J.L. Desmarais
      !> 15_10_2014 - interface modifications - J.L. Desmarais
      !> (to have a unique sd_operators, p_model... shared between the
      !>  field and the buffer layer objects)
      !-----------------------------------------------------------------
      module bf_compute_class

        use bc_operators_class, only :
     $       bc_operators

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_bf_layer_data

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       x_min, x_max,
     $       y_min, y_max,
     $       bc_size,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators


        implicit none

        private
        public :: bf_compute


        !>@class bf_compute
        !> class encapsulating the main temporary tables for the time 
        !> integration of the buffer layer
        !
        !>@param alignment_tmp
        !> temporary alignment saving the alignment of the buffer layer at
        !> the previous time step
        !
        !>@param grdpts_id_tmp
        !> temporary array saving the grdpts_id at the previous time step
        !
        !>@param x_map
        !> temporary array saving the x-coordinates of the buffer layer at
        !> the previous time step
        !
        !>@param y_map
        !> temporary array saving the y-coordinates of the buffer layer at
        !> the previous time step
        !
        !>@param grdpts_id_tmp
        !> temporary array saving the grdpts_id at the previous time step
        !
        !>@param nodes_tmp
        !> temporary array whose size is the same as the array containing
        !> the grid points integrated in time
        !
        !>@param time_dev
        !> temporary array containing the time derivatives and whose size
        !> is the same as the array containing the grid points integrated
        !> in time
        !
        !>@param does_previous_timestep_exist
        !> determine whether the previous time step was stored in this object
        !
        !>@param allocate_tables
        !> allocate the grdpts_id_tmp, nodes_tmp and time_dev attributes
        !
        !>@param deallocate_tables
        !> deallocate the grdpts_id_tmp, nodes_tmp and time_dev attributes
        !
        !>@param compute_time_dev
        !> compute the time_dev attribute
        !
        !>@param compute_integration_step
        !> compute the nodes and nodes_tmp using the integration
        !> procedure
        !---------------------------------------------------------------
        type :: bf_compute

          integer(ikind), dimension(:,:)  , allocatable :: alignment_tmp
          integer       , dimension(:,:)  , allocatable :: grdpts_id_tmp
          real(rkind)   , dimension(:)    , allocatable :: x_map_tmp
          real(rkind)   , dimension(:)    , allocatable :: y_map_tmp
          real(rkind)   , dimension(:,:,:), allocatable :: nodes_tmp
          real(rkind)   , dimension(:,:,:), allocatable :: time_dev

          contains

          procedure,   pass :: does_previous_timestep_exist
                       
          procedure,   pass :: allocate_tables
          procedure,   pass :: deallocate_tables
                       
          procedure,   pass :: compute_time_dev
          procedure,   pass :: compute_integration_step

          procedure,   pass :: get_grdpts_id_part

          procedure,   pass :: set_alignment   !only for tests
          procedure,   pass :: set_grdpts_id   !only for tests
          procedure,   pass :: set_x_map       !only for tests
          procedure,   pass :: set_y_map       !only for tests
          procedure,   pass :: set_nodes       !only for tests
          procedure,   pass :: get_time_dev    !only for tests

        end type bf_compute

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate the nodes_tmp and time_dev tables storing
        !> the intermediate time integration steps for the
        !> buffer layer
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param exist
        !> say whether the previous step is stored in the object
        !--------------------------------------------------------------
        function does_previous_timestep_exist(
     $       this)
     $       result(exist)

          implicit none

          class(bf_compute), intent(in) :: this
          logical                       :: exist

          exist = allocated(this%nodes_tmp)
        
        end function does_previous_timestep_exist


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate the nodes_tmp and time_dev tables storing
        !> the intermediate time integration steps for the
        !> buffer layer
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param size_x
        !> size of the nodes tables for the buffer layer integrated
        !> along the x-direction
        !
        !>@param size_y
        !> size of the nodes tables for the buffer layer integrated
        !> along the y-direction
        !
        !>@param grdpts_id
        !> grdpts_id from the buffer layer at t
        !
        !>@param alignment
        !> alignment of the buffer layer at t
        !--------------------------------------------------------------
        subroutine allocate_tables(
     $       this,
     $       size_x,
     $       size_y,
     $       alignment,
     $       x_map,
     $       y_map,
     $       grdpts_id)

          implicit none

          class(bf_compute)                , intent(inout) :: this
          integer(ikind)                   , intent(in)    :: size_x
          integer(ikind)                   , intent(in)    :: size_y
          integer(ikind)   , dimension(2,2), intent(in)    :: alignment
          real(rkind)      , dimension(:)  , intent(in)    :: x_map
          real(rkind)      , dimension(:)  , intent(in)    :: y_map
          integer          , dimension(:,:), intent(in)    :: grdpts_id

          allocate(this%alignment_tmp(2,2))
          allocate(this%x_map_tmp(size_x))
          allocate(this%y_map_tmp(size_y))
          allocate(this%grdpts_id_tmp(size_x,size_y))
          allocate(this%nodes_tmp(size_x,size_y,ne))
          allocate(this%time_dev(size_x,size_y,ne))

          this%alignment_tmp(:,:) = alignment(:,:)
          this%x_map_tmp(:)       = x_map(:)
          this%y_map_tmp(:)       = y_map(:)
          this%grdpts_id_tmp(:,:) = grdpts_id(:,:)

        end subroutine allocate_tables


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the nodes_tmp and time_dev tables storing
        !> the intermediate time integration steps for the
        !> buffer layer
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine deallocate_tables(this)

          implicit none

          class(bf_compute), intent(inout) :: this

          if(allocated(this%alignment_tmp)) then
             deallocate(this%alignment_tmp)
             deallocate(this%x_map_tmp)
             deallocate(this%y_map_tmp)
             deallocate(this%grdpts_id_tmp)
             deallocate(this%nodes_tmp)
             deallocate(this%time_dev)

          end if

        end subroutine deallocate_tables        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the grid points of the
        !> buffer layer
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param td_operators_used
        !> object encapsulating the functions computing the time
        !> derivatives of the grid points
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param x_map
        !> coordinates along the x-direction of the buffer layer
        !> grid points
        !
        !>@param y_map
        !> coordinates along the y-direction of the buffer layer
        !> grid points
        !
        !>@param s
        !> object encapsulating the space discretisation methods
        !
        !>@param p_model
        !> object encapsulating the governing equations of the
        !> physical model
        !
        !>@param bc_used
        !> object encapsulating the boundary conditions
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t, nodes, x_map, y_map,
     $     s,
     $     p_model,bc_used,
     $     bf_alignment,
     $     grdpts_id,
     $     interior_nodes,
     $     bc_sections,
     $     x_borders, y_borders)

          implicit none

          class(bf_compute)                          , intent(inout) :: this
          type(td_operators)                         , intent(in)    :: td_operators_used
          real(rkind)                                , intent(in)    :: t
          real(rkind), dimension(:,:,:)              , intent(in)    :: nodes
          real(rkind), dimension(:)                  , intent(in)    :: x_map
          real(rkind), dimension(:)                  , intent(in)    :: y_map
          type(sd_operators)                         , intent(in)    :: s
          type(pmodel_eq)                            , intent(in)    :: p_model
          type(bc_operators)                         , intent(in)    :: bc_used
          integer(ikind), dimension(2,2)             , intent(in)    :: bf_alignment
          integer       , dimension(:,:)             , intent(in)    :: grdpts_id
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections
          integer(ikind), dimension(2)               , intent(in)    :: x_borders
          integer(ikind), dimension(2)               , intent(in)    :: y_borders
          
          call td_operators_used%compute_time_dev_nopt(
     $         t,nodes,x_map,y_map,
     $         s,p_model,bc_used,
     $         this%time_dev,
     $         bf_alignment,
     $         grdpts_id,
     $         interior_nodes,
     $         bc_sections,
     $         x_borders=x_borders,
     $         y_borders=y_borders)

        end subroutine compute_time_dev


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
        subroutine compute_integration_step(
     $     this,
     $     grdpts_id, nodes, dt,
     $     x_borders, y_borders,
     $     integration_step_nopt)

          implicit none

          class(bf_compute)            , intent(inout) :: this
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: dt
          integer(ikind), dimension(2) , intent(in)    :: x_borders
          integer(ikind), dimension(2) , intent(in)    :: y_borders
          procedure(timeInt_step_nopt)                 :: integration_step_nopt

          call integration_step_nopt(
     $         nodes,
     $         dt,
     $         this%nodes_tmp,
     $         this%time_dev,
     $         grdpts_id,
     $         x_borders=x_borders,
     $         y_borders=y_borders)

        end subroutine compute_integration_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract partially the grdpts_id
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param tmp_grdpts_id
        !> array of grdpts_id extracted
        !
        !>@param gen_coords
        !> SW and NE corners of the grdpts_id extracted expressed
        !> using general coordinates
        !--------------------------------------------------------------
        subroutine get_grdpts_id_part(
     $     this,
     $     tmp_grdptsid,
     $     gen_coords)

          implicit none

          class(bf_compute)              , intent(in)    :: this
          integer        , dimension(:,:), intent(inout) :: tmp_grdptsid
          integer(ikind) , dimension(2,2), intent(in)    :: gen_coords

          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j


          if(allocated(this%alignment_tmp)) then

             !get the synchronization indices
             call get_indices_to_extract_bf_layer_data(
     $            this%alignment_tmp,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)
          
             !fill the grid points asked
             do j=1, size_y
                do i=1, size_x
                   
                   tmp_grdptsid(i_recv+i-1,j_recv+j-1) =
     $                  this%grdpts_id_tmp(i_send+i-1,j_send+j-1)
                   
                end do
             end do

          end if

        end subroutine get_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the alignment attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param alignment
        !> alignment of the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_alignment(this, alignment)

          implicit none

          class(bf_compute)                   , intent(inout) :: this
          integer, dimension(:,:), allocatable, intent(inout) :: alignment

          call MOVE_ALLOC(alignment,this%alignment_tmp)

        end subroutine set_alignment


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the alignment attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the gridpoints at t-dt
        !--------------------------------------------------------------
        subroutine set_grdpts_id(this, grdpts_id)

          implicit none

          class(bf_compute)                   , intent(inout) :: this
          integer, dimension(:,:), allocatable, intent(inout) :: grdpts_id

          call MOVE_ALLOC(grdpts_id,this%grdpts_id_tmp)

        end subroutine set_grdpts_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the x_map attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_map
        !> x-coordinates for the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_x_map(this, x_map)

          implicit none

          class(bf_compute)                     , intent(inout) :: this
          real(rkind), dimension(:), allocatable, intent(inout) :: x_map

          call MOVE_ALLOC(x_map,this%x_map_tmp)

        end subroutine set_x_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the y_map attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_map
        !> y-coordinates for the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_y_map(this, y_map)

          implicit none

          class(bf_compute)                     , intent(inout) :: this
          real(rkind), dimension(:), allocatable, intent(inout) :: y_map

          call MOVE_ALLOC(y_map,this%y_map_tmp)

        end subroutine set_y_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the alignment attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> grid-points data at t-dt
        !--------------------------------------------------------------
        subroutine set_nodes(this, nodes)

          implicit none

          class(bf_compute)                         , intent(inout) :: this
          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: nodes

          call MOVE_ALLOC(nodes,this%nodes_tmp)

        end subroutine set_nodes

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the time_dev attribute
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param time_dev
        !> array which the time derivatives of the buffer layer
        !--------------------------------------------------------------
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_compute)                         , intent(in)  :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: time_dev


          if(allocated(this%time_dev)) then
             allocate(time_dev(
     $            size(this%time_dev,1),
     $            size(this%time_dev,2),
     $            size(this%time_dev,3)))

             time_dev = this%time_dev

          else
             print '(''bf_compute_class'')'
             print '(''get_time_dev'')'
             stop 'time dev not allocated'

          end if

        end subroutine get_time_dev        

      end module bf_compute_class
