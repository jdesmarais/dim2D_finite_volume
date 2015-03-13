      !> @file
      !> class encapsulating the intermediate tables for the time 
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
      !> 13_03_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_compute_basic_class

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: bf_compute_basic


        !>@class bf_compute_basic
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
        !---------------------------------------------------------------
        type :: bf_compute_basic

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

          procedure,   pass :: set_alignment   !only for tests
          procedure,   pass :: set_grdpts_id   !only for tests
          procedure,   pass :: set_x_map       !only for tests
          procedure,   pass :: set_y_map       !only for tests
          procedure,   pass :: set_nodes       !only for tests
          procedure,   pass :: get_time_dev    !only for tests

        end type bf_compute_basic

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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param exist
        !> say whether the previous step is stored in the object
        !--------------------------------------------------------------
        function does_previous_timestep_exist(
     $       this)
     $       result(exist)

          implicit none

          class(bf_compute_basic), intent(in) :: this
          logical                             :: exist

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
        !> bf_compute_basic object encapsulating the main
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

          class(bf_compute_basic)          , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine deallocate_tables(this)

          implicit none

          class(bf_compute_basic), intent(inout) :: this

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
        !> set the alignment attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param alignment
        !> alignment of the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_alignment(this, alignment)

          implicit none

          class(bf_compute_basic)             , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the gridpoints at t-dt
        !--------------------------------------------------------------
        subroutine set_grdpts_id(this, grdpts_id)

          implicit none

          class(bf_compute_basic)             , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_map
        !> x-coordinates for the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_x_map(this, x_map)

          implicit none

          class(bf_compute_basic)               , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_map
        !> y-coordinates for the buffer layer at t-dt
        !--------------------------------------------------------------
        subroutine set_y_map(this, y_map)

          implicit none

          class(bf_compute_basic)               , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> grid-points data at t-dt
        !--------------------------------------------------------------
        subroutine set_nodes(this, nodes)

          implicit none

          class(bf_compute_basic)                   , intent(inout) :: this
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
        !> bf_compute_basic object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param time_dev
        !> array which the time derivatives of the buffer layer
        !--------------------------------------------------------------
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_compute_basic)                   , intent(in)  :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: time_dev


          if(allocated(this%time_dev)) then
             allocate(time_dev(
     $            size(this%time_dev,1),
     $            size(this%time_dev,2),
     $            size(this%time_dev,3)))

             time_dev = this%time_dev

          else
             print '(''bf_compute_basic_class'')'
             print '(''get_time_dev'')'
             stop 'time dev not allocated'

          end if

        end subroutine get_time_dev        

      end module bf_compute_basic_class
