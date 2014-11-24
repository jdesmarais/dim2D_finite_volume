      !> @file
      !> module implementing the object encapsulating the position of
      !> the increasing detectors and the subroutines extending the
      !> computational domain if they are triggered
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating the position of
      !> the increasing detectors and the subroutines extending the
      !> computational domain if they are triggered
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_icr_class

        use bf_detector_dcr_list_class, only :
     $     bf_detector_dcr_list

        use bf_detector_dcr_list_N_class, only :
     $       bf_detector_dcr_list_N

        use bf_detector_dcr_list_S_class, only :
     $       bf_detector_dcr_list_S

        use bf_detector_dcr_list_E_class, only :
     $       bf_detector_dcr_list_E

        use bf_detector_dcr_list_W_class, only :
     $       bf_detector_dcr_list_W

        use bf_detector_icr_list_class, only :
     $       bf_detector_icr_list

        use bf_path_icr_class, only :
     $       bf_path_icr

        use bf_nbc_template_module
 
        use bf_sublayer_class, only :
     $       bf_sublayer

        use bf_interface_class, only :
     $       bf_interface

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       dct_icr_distance,
     $       dct_icr_N_default,
     $       dct_icr_S_default,
     $       dct_icr_E_default,
     $       dct_icr_W_default

        use parameters_constant, only :
     $       N,S,E,W,interior

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       dt,search_nb_dt

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_interface_icr


        !>@class bf_interface_icr
        !> class encapsulating the position of the increasing detectors
        !> and the subroutines extending the interior domain
        !
        !>@param N_detectors_list
        !> list containing the position of the north detectors
        !
        !>@param S_detectors_list
        !> list containing the position of the south detectors
        !
        !>@param E_detectors_list
        !> list containing the position of the east detectors
        !
        !>@param W_detectors_list
        !> list containing the position of the west detectors
        !
        !>@param ini
        !> initialize the position of the increasing detectors
        !> and the parent object
        !
        !>@param get_modified_grdpts_list
        !> from a detector position, get a list of bc_interior_pt
        !> activated
        !
        !>@param process_idetector_list
        !> process the list of current detectors, modify the buffer
        !> layers. If they are activated by the detectors and
        !> determine the new list of detectors
        !
        !>@param combine_bf_idetector_lists
        !> connect the bf_detector_icr_list objects and determine the
        !> new detector lists
        !
        !>@param update_bf_layers_with_idetectors
        !> update the size of the computational domain based on the
        !> increasing detector activations
        !
        !>@param update_icr_detectors_after_removal
        !> update the position of the increasing detectors
        !> if a sublayer is removed
        !
        !>@param remove_sublayer
        !> remove a sublayer
        !
        !>@param is_detector_icr_activated
        !> check whether the detector is activated
        !
        !>@param get_central_grdpt
        !> check whether the detector is activated
        !
        !>@param check_neighboring_bc_interior_pts
        !> check the identity of the grid points surrounding a central
        !> point: this function encapsulates the function used if the
        !> central point is at the interface between the interior and
        !> the buffer layers and the function checking the grid points
        !> in case all the grid points are inside a buffer layer
        !
        !>@param check_neighboring_bc_interior_pts_for_interior
        !> check the neighboring bc_interior_pt around the central
        !> point identified by its general coordinates
        !
        !>@param is_inside_border_layer
        !> check the identity of the grid points surrounding a central
        !> point: this function encapsulates the function used if the
        !> central point is at the interface between the interior and
        !> the buffer layers and the function checking the grid points
        !> in case all the grid points are inside a buffer layer
        !
        !>@param create_nbc_interior_pt_template
        !> create the template for the neighboring points around
        !> the central point asked: as the central is located in
        !> a layer at the interface between the interior points
        !> and the boundary layers, it is required to initialize
        !> the template using the coordinates as if there were no
        !> boundary layers, then data are exchanged with the
        !> neighboring buffer layers
        !
        !>@param check_nbc_interior_pt_template
        !> check if the grid points around the center point are
        !> bc_interior_pt
        !
        !>@param check_bc_interior_pt
        !> check whether the grid point tested is a bc_interior_pt
        !> and if so save the general coordinates of the grid point
        !> in mgrdpts
        !
        !>@param update_grdpts_id_for_template
        !> update the position of the increasing detectors
        !> if a sublayer is removed
        !
        !>@param print_idetectors_on
        !> print the increasing detector positions on a matrix
        !
        !>@param print_idetectors_on_binary
        !> print the increasing detector positions on binary
        !> output files
        !---------------------------------------------------------------
        type, extends(bf_interface) :: bf_interface_icr

          integer(ikind), dimension(:,:), allocatable, private :: N_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: S_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: E_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: W_detectors_list

          real(rkind)   , dimension(:,:), allocatable, private :: N_detectors_list_coords
          real(rkind)   , dimension(:,:), allocatable, private :: S_detectors_list_coords
          real(rkind)   , dimension(:,:), allocatable, private :: E_detectors_list_coords
          real(rkind)   , dimension(:,:), allocatable, private :: W_detectors_list_coords

          contains

          procedure,   pass :: ini
          procedure,   pass :: get_modified_grdpts_list
          procedure,   pass :: process_idetector_list
          procedure,   pass :: combine_bf_idetector_lists
          procedure,   pass :: update_bf_layers_with_idetectors

          procedure,   pass :: update_icr_detectors_after_removal
          procedure,   pass :: remove_sublayer

          procedure, nopass, private :: is_detector_icr_activated
          procedure, nopass, private :: get_central_grdpt
          procedure,   pass, private :: check_neighboring_bc_interior_pts
          procedure,   pass, private :: check_neighboring_bc_interior_pts_for_interior
          procedure, nopass, private :: is_inside_border_layer
          procedure,   pass          :: create_nbc_interior_pt_template
          procedure, nopass          :: check_nbc_interior_pt_template
          procedure, nopass, private :: check_bc_interior_pt
          procedure,   pass, private :: update_grdpts_id_for_template

          procedure,   pass          :: print_idetectors_on
          procedure,   pass          :: print_idetectors_on_binary

        end type bf_interface_icr


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the position of the increasing detectors
        !> and the parent object
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(bf_interface_icr), intent(inout) :: this

          integer(ikind) :: i


          !initialize the parent attributes
          call this%bf_interface%ini()


          !intialize the attributes specific to bf_interface_icr
          !list with the coordinates of the detectors as (i,j)
          allocate(this%N_detectors_list(2,nx-2*(bc_size+dct_icr_distance)+2))
          allocate(this%S_detectors_list(2,nx-2*(bc_size+dct_icr_distance)+2))
          allocate(this%E_detectors_list(2,ny-2*(bc_size+dct_icr_distance)))
          allocate(this%W_detectors_list(2,nx-2*(bc_size+dct_icr_distance)))

          !list with the coordinates of the detectors as (x,y)
          do i=bc_size+dct_icr_distance, nx-(bc_size+dct_icr_distance)+1
             this%S_detectors_list(1,i-(bc_size+dct_icr_distance)+1) = i
             this%S_detectors_list(2,i-(bc_size+dct_icr_distance)+1) = dct_icr_S_default
          end do

          do i=bc_size+dct_icr_distance, nx-(bc_size+dct_icr_distance)+1
             this%N_detectors_list(1,i-(bc_size+dct_icr_distance)+1) = i
             this%N_detectors_list(2,i-(bc_size+dct_icr_distance)+1) = dct_icr_N_default
          end do

          do i=bc_size+dct_icr_distance+1, ny-(bc_size+dct_icr_distance)
             this%W_detectors_list(1,i-(bc_size+dct_icr_distance)) = dct_icr_W_default
             this%W_detectors_list(2,i-(bc_size+dct_icr_distance)) = i
          end do

          do i=bc_size+dct_icr_distance+1, ny-(bc_size+dct_icr_distance)
             this%E_detectors_list(1,i-(bc_size+dct_icr_distance)) = dct_icr_E_default
             this%E_detectors_list(2,i-(bc_size+dct_icr_distance)) = i
          end do          

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the size of the computational domain based on the
        !> increasing detector activations
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time
        !
        !> @param dt
        !> time step
        !
        !> @param interior_x_map
        !> table encapsulating the coordinates along the x-axis
        !
        !> @param interior_y_map
        !> table encapsulating the coordinates along the y-axis
        !
        !> @param interior_nodes_0
        !> table encapsulating the data of the grid points of the
        !> interior domain at the previous time step (t-dt)
        !
        !> @param interior_nodes_1
        !> table encapsulating the data of the grid points of the
        !> interior domain at the current time step (t)
        !--------------------------------------------------------------
        subroutine update_bf_layers_with_idetectors(
     $     this,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(bf_interface_icr)         , intent(inout) :: this
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1


          type(bf_path_icr)              :: path_update_idetectors
          integer(ikind), dimension(2)   :: cpt_coords_p
          type(bf_detector_icr_list)     :: N_ndt_list
          type(bf_detector_icr_list)     :: S_ndt_list
          type(bf_detector_icr_list)     :: E_ndt_list
          type(bf_detector_icr_list)     :: W_ndt_list


          !initialization of the path that will gather information
          !on the buffer layers to be updated: the path collects
          !the updating information belonging to a specific buffer
          !layer and when another buffer layer is activated by
          !the detectors the path applies all the changes to the
          !previous buffer layer and starts collecting information
          !on the next buffer layer
          call path_update_idetectors%ini()

          !there are 4 layers of detectors corresponding to the 4
          !cardinal points. They create a closed path of detectors
          !we update in this order:
          !1: S, 2: E, 3: W, 4: N
          !in this way the path created by one buffer layer can
          !be continued by the next list of detectors and prevent
          !double operations on buffer layers
          
          !0) initialization of the first point indicating 
          !   where the previous neighboring bc_interior_pt
          !   were analyzed
          cpt_coords_p = [nx/2, ny/2]

          !1) South detectors
          if(allocated(this%S_detectors_list)) then
             call S_ndt_list%ini(S, size(this%S_detectors_list,2))

             call process_idetector_list(
     $            this,
     $            this%S_detectors_list,
     $            S_ndt_list,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            cpt_coords_p,
     $            path_update_idetectors)

          end if

          !2) East detectors
          if(allocated(this%E_detectors_list)) then
             call E_ndt_list%ini(E, size(this%E_detectors_list,2))

             call process_idetector_list(
     $            this,
     $            this%E_detectors_list,
     $            E_ndt_list,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            cpt_coords_p,
     $            path_update_idetectors)

          end if

          !3) West detectors
          if(allocated(this%W_detectors_list)) then
             call W_ndt_list%ini(W, size(this%W_detectors_list,2))

             call process_idetector_list(
     $            this,
     $            this%W_detectors_list,
     $            W_ndt_list,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            cpt_coords_p,
     $            path_update_idetectors)

          end if

          !4) North detectors
          if(allocated(this%N_detectors_list)) then
             call N_ndt_list%ini(N, size(this%N_detectors_list,2))

             call process_idetector_list(
     $            this,
     $            this%N_detectors_list,
     $            N_ndt_list,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            cpt_coords_p,
     $            path_update_idetectors)

          end if

          !5) process the last path as after the detectors of
          !   the last main layer are analysed, it is possible
          !   that the final part of the path has not been
          !   processed
          if(path_update_idetectors%get_nb_pts().gt.0) then
             call path_update_idetectors%process_path(
     $            this,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1)
          end if


          !6) reconnect the detector paths
          !   the buffer layer have been updated
          !   we now have several detector lists but they do
          !   not make a closed path. They are now reconnected
          call combine_bf_idetector_lists(
     $         this,
     $         N_ndt_list, S_ndt_list,
     $         E_ndt_list, W_ndt_list)

        end subroutine update_bf_layers_with_idetectors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> connect the bf_detector_icr_list objects and determine the
        !> new detector lists
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param N_idetectors_list
        !> temporary object where the position of the new north
        !> increasing detectors is stored
        !
        !>@param S_idetectors_list
        !> temporary object where the position of the new south
        !> increasing detectors is stored
        !
        !>@param E_idetectors_list
        !> temporary object where the position of the new east
        !> increasing detectors is stored
        !
        !>@param W_idetectors_list
        !> temporary object where the position of the new west
        !> increasing detectors is stored
        !--------------------------------------------------------------
        subroutine combine_bf_idetector_lists(
     $     this,
     $     N_idetectors_list, S_idetectors_list,
     $     E_idetectors_list, W_idetectors_list)

          implicit none

          class(bf_interface_icr) , intent(inout) :: this
          type(bf_detector_icr_list), intent(in)    :: N_idetectors_list
          type(bf_detector_icr_list), intent(in)    :: S_idetectors_list
          type(bf_detector_icr_list), intent(in)    :: E_idetectors_list
          type(bf_detector_icr_list), intent(in)    :: W_idetectors_list


          integer(ikind), dimension(:,:), allocatable :: N_idetectors_list_n
          integer(ikind), dimension(:,:), allocatable :: S_idetectors_list_n
          integer(ikind), dimension(:,:), allocatable :: E_idetectors_list_n
          integer(ikind), dimension(:,:), allocatable :: W_idetectors_list_n


          integer(ikind), dimension(2) :: n1_coords, n2_coords, inter_coords
          integer                      :: n1_inter_nb, n2_inter_nb
          real(rkind)                  :: n1_x_change, n1_y_change
          real(rkind)                  :: n2_x_change, n2_y_change
          integer                      :: k


          !S detectors recombination
          n1_coords = W_idetectors_list%get_head()
          call S_idetectors_list%get_inter_detector_param(
     $         n1_coords, S_idetectors_list%get_head(),
     $         n1_x_change, n1_y_change, n1_inter_nb)

          n2_coords = S_idetectors_list%get_tail()
          call S_idetectors_list%get_inter_detector_param(
     $         n2_coords, E_idetectors_list%get_head(),
     $         n2_x_change, n2_y_change, n2_inter_nb)

          allocate(S_idetectors_list_n(
     $         2, S_idetectors_list%get_nb_detectors()+n1_inter_nb+n2_inter_nb))

          do k=1, n1_inter_nb
             inter_coords = S_idetectors_list%get_inter_detector_coords(
     $            n1_coords,
     $            n1_x_change, n1_y_change, k)
             S_idetectors_list_n(:,k) = inter_coords
          end do

          call S_idetectors_list%fill_new_detector_table(
     $         n1_inter_nb+1, S_idetectors_list_n)

          do k=1, n2_inter_nb
             inter_coords = S_idetectors_list%get_inter_detector_coords(
     $            n2_coords,
     $            n2_x_change, n2_y_change, k)
             S_idetectors_list_n(
     $            :,n1_inter_nb+S_idetectors_list%get_nb_detectors()+k) =
     $            inter_coords
          end do


          !N detectors recombination
          n1_coords = W_idetectors_list%get_tail()
          call N_idetectors_list%get_inter_detector_param(
     $         n1_coords, N_idetectors_list%get_head(),
     $         n1_x_change, n1_y_change, n1_inter_nb)

          n2_coords = N_idetectors_list%get_tail()
          call N_idetectors_list%get_inter_detector_param(
     $         n2_coords, E_idetectors_list%get_tail(),
     $         n2_x_change, n2_y_change, n2_inter_nb)

          allocate(N_idetectors_list_n(
     $         2, N_idetectors_list%get_nb_detectors()+n1_inter_nb+n2_inter_nb))

          do k=1, n1_inter_nb
             inter_coords = N_idetectors_list%get_inter_detector_coords(
     $            n1_coords,
     $            n1_x_change, n1_y_change, k)
             N_idetectors_list_n(:,k) = inter_coords
          end do

          call N_idetectors_list%fill_new_detector_table(
     $         n1_inter_nb+1, N_idetectors_list_n)

          do k=1, n2_inter_nb
             inter_coords = N_idetectors_list%get_inter_detector_coords(
     $            n2_coords,
     $            n2_x_change, n2_y_change, k)
             N_idetectors_list_n(
     $            :,n1_inter_nb+N_idetectors_list%get_nb_detectors()+k) =
     $            inter_coords
          end do

          
          !E detectors recombination
          allocate(E_idetectors_list_n(2, E_idetectors_list%get_nb_detectors()))
          call E_idetectors_list%fill_new_detector_table(
     $         1, E_idetectors_list_n)
          

          !W detectors recombination
          allocate(W_idetectors_list_n(2, W_idetectors_list%get_nb_detectors()))
          call W_idetectors_list%fill_new_detector_table(
     $         1, W_idetectors_list_n)

          !move the allocations of the previous detector tables
          !to the new detector tables
          call MOVE_ALLOC(N_idetectors_list_n, this%N_detectors_list)
          call MOVE_ALLOC(S_idetectors_list_n, this%S_detectors_list)
          call MOVE_ALLOC(E_idetectors_list_n, this%E_detectors_list)
          call MOVE_ALLOC(W_idetectors_list_n, this%W_detectors_list)
             
        end subroutine combine_bf_idetector_lists      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> process the list of current detectors, modify the buffer layers
        !> if they are activated by the detectors and determine the new list
        !> of detectors
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param dt_list
        !> detector list processed
        !
        !>@param ndt_list
        !> new detector list resulting from the position update
        !> of the detector list
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param dx
        !> grid size along the x-direction
        !
        !>@param dy
        !> grid size along the y-direction
        !
        !>@param cpt_coords_p
        !> general coordinates of the central point triggered by the
        !> last detector
        !
        !>@param path
        !> bf_path_icr object gathering the data when grid points are
        !> activated by the detectors
        !--------------------------------------------------------------
        subroutine process_idetector_list(
     $     this,
     $     dt_list,
     $     ndt_list,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     cpt_coords_p,
     $     path)
        
          implicit none

          class(bf_interface_icr)                      , intent(inout) :: this
          integer(ikind)          , dimension(:,:)     , intent(in)    :: dt_list
          type(bf_detector_icr_list)                   , intent(inout) :: ndt_list
          type(pmodel_eq)                              , intent(in)    :: p_model
          real(rkind)                                  , intent(in)    :: t
          real(rkind)                                  , intent(in)    :: dt
          real(rkind)             , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)             , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)             , dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind)             , dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          integer(ikind)          , dimension(2)       , intent(inout) :: cpt_coords_p
          type(bf_path_icr)                            , intent(inout) :: path

          integer(ikind), dimension(2)   :: cpt_coords
          integer                        :: k,l
          integer                        :: nb_mgrdpts
          integer       , dimension(2,9) :: mgrdpts

          !loop over the detectors
          do k=1, size(dt_list,2)

             !extract the list of bc_interior_pt that should
             !be turned into interior_pt due to the activation
             !of the detector k
             call get_modified_grdpts_list(
     $            this, dt_list(:,k),
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes1,
     $            p_model,
     $            cpt_coords_p, cpt_coords,
     $            nb_mgrdpts, mgrdpts, ndt_list)

             !the point used as center point to determine the
             !neighboring points in the analysis of bc_interior_pt
             !is updated
             cpt_coords_p = cpt_coords

             !loop over the list of grid points to be turned from
             !bc_interior_pt to interior_pt
             do l=1, nb_mgrdpts
                
                !gathering the changes to a buffer layer
                !(allocate/reallocate/merge)
                call path%analyze_pt(mgrdpts(:,l), this)

                !if the grid point that is currently analyzed
                !leads to changes on a buffer layer which is
                !independent of the buffer layer investigated
                !up to now, all the changes are implemented
                !on the current buffer layer
                if(path%is_ended()) then

                   !the buffer layers are updated
                   call path%process_path(
     $                  this,
     $                  p_model,
     $                  t,dt,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes0,
     $                  interior_nodes1)

                   !the grid point that led to the path end
                   !is used to reinitialize the current path
                   call path%analyze_pt(mgrdpts(:,l), this)

                end if

             end do
          end do

        end subroutine process_idetector_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> from a detector position, get a list of bc_interior_pt activated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param d_coords
        !> detector general coordinates
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param cpt_coords_p
        !> general coordinates of the central point triggered by the
        !> last detector
        !
        !>@param nb_mgrdpts
        !> number of grid points to be modified
        !
        !>@param mgrdpts
        !> list of the bc_interior_pt grid points to be modified
        !
        !>@param ndt_list
        !> temporary new detector list after their position update
        !--------------------------------------------------------------
        subroutine get_modified_grdpts_list(
     $     this,
     $     d_coords,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     cpt_coords_p,
     $     cpt_coords,
     $     nb_mgrdpts,
     $     mgrdpts,
     $     ndt_list)

          implicit none

          class(bf_interface_icr)         , intent(inout) :: this
          integer(ikind), dimension(2)    , intent(in)    :: d_coords
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer(ikind), dimension(2)    , intent(in)    :: cpt_coords_p
          integer(ikind), dimension(2)    , intent(out)   :: cpt_coords
          integer                         , intent(out)   :: nb_mgrdpts
          integer(ikind), dimension(2,9)  , intent(out)   :: mgrdpts
          type(bf_detector_icr_list)        , intent(inout) :: ndt_list


          real(rkind), dimension(ne)   :: node_var
          real(rkind)   , dimension(2) :: velocity
          integer(ikind), dimension(2) :: d_coords_n


          !initialization of the number of modified grid points
          nb_mgrdpts = 0


          !extract the nodes at the coordinates of the detector
          node_var = this%get_nodes(d_coords, interior_nodes)


          !if the detector is activated, then we check
          !whether grid points need to be modified
          if(is_detector_icr_activated(node_var, p_model)) then

             !extract the velocity at the coordinates of the detector
             velocity = p_model%get_velocity(node_var)
             

             !get the first point from which we should look for a
             !bc_interior_pt to be activated and the new coordinates
             !from the detector
             cpt_coords = get_central_grdpt(
     $            d_coords,
     $            velocity,
     $            interior_x_map,
     $            interior_y_map,
     $            d_coords_n)
             
             !add the new coordinates of the detector of the ndt_list
             call ndt_list%add_new_detector(d_coords_n)

             !look for a bc_interior_pt around the point previously
             !computed whose coordinates are: cpt_coords
             !we make use of the previously checked neighboring points
             !whose center was cpt_coords_p to reduce the number of
             !grid points checked
             call check_neighboring_bc_interior_pts(
     $            this,
     $            cpt_coords_p,
     $            cpt_coords,
     $            nb_mgrdpts,
     $            mgrdpts)

          !otherwise, the coordinates of the new detector are simply
          !the previous ones, and are saved in the ndt_list
          else
             call ndt_list%add_new_detector(d_coords)
          end if

        end subroutine get_modified_grdpts_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the detector is activated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_var
        !> governing variables at the grid point
        !
        !>@return activated
        !> logical stating whether the detector is activated
        !--------------------------------------------------------------
        function is_detector_icr_activated(nodes, p_model)
     $     result(activated)
        
          implicit none
          
          real(rkind), dimension(ne), intent(in) :: nodes
          type(pmodel_eq)           , intent(in) :: p_model
          logical                                :: activated
          
          activated = p_model%are_openbc_undermined(nodes)

        end function is_detector_icr_activated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the detector is activated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param d_coords
        !> detector general coordinates
        !
        !>@param velocity
        !> velocity vector at the detector position
        !
        !>@param dx
        !> grid size along the x-direction
        !
        !>@param dy
        !> grid size along the y-direction
        !
        !>@param d_coords_n
        !> new detector general coordinates
        !
        !>@return cpt_coords
        !> general coordinates of the grid point activated by the
        !> detector
        !--------------------------------------------------------------
        function get_central_grdpt(
     $     d_coords,
     $     velocity,
     $     interior_x_map,
     $     interior_y_map,
     $     d_coords_n)
     $     result(cpt_coords)

          implicit none

          integer(ikind), dimension(2) , intent(in)  :: d_coords
          real(rkind)   , dimension(2) , intent(in)  :: velocity
          real(rkind)   , dimension(nx), intent(in)  :: interior_x_map
          real(rkind)   , dimension(ny), intent(in)  :: interior_y_map
          integer(ikind), dimension(2) , intent(out) :: d_coords_n
          integer(ikind), dimension(2)               :: cpt_coords

          real(rkind) :: dir_x, dir_y
          real(rkind) :: dx,dy

          dx = interior_x_map(2)-interior_x_map(1)
          dy = interior_y_map(2)-interior_y_map(1)

          !1) get the direction to look for a bc_interior_pt
          !dir_x  = velocity(1)*search_nb_dt*dt/dx
          !dir_y  = velocity(2)*search_nb_dt*dt/dy

          dir_x  = velocity(1)*search_nb_dt*dt/dx
          dir_y  = velocity(2)*search_nb_dt*dt/dy

          
          if(rkind.eq.4) then

             !2) get the point indices in the direction given
             !   by the velocity vector
             cpt_coords(1) = d_coords(1) + nint(dir_x)
             cpt_coords(2) = d_coords(2) + nint(dir_y)
             
             !3) compute the new detector position
             d_coords_n(1) = d_coords(1) + nint(max(min(dir_x,1.0d0),-1.0d0))
             d_coords_n(2) = d_coords(2) + nint(max(min(dir_y,1.0d0),-1.0d0))

          else

             !2) get the point indices in the direction given
             !   by the velocity vector
             cpt_coords(1) = d_coords(1) + idnint(dir_x)
             cpt_coords(2) = d_coords(2) + idnint(dir_y)
             
             !3) compute the new detector position
             d_coords_n(1) = d_coords(1) + idnint(max(min(dir_x,1.0d0),-1.0d0))
             d_coords_n(2) = d_coords(2) + idnint(max(min(dir_y,1.0d0),-1.0d0))

          end if

        end function get_central_grdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the identity of the grid points surrounding a central
        !> point: this function encapsulates the function used if the
        !> central point is at the interface between the interior and
        !> the buffer layers and the function checking the grid points
        !> in case all the grid points are inside a buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param cpt_coords_p
        !> previous central point whose neighboring grid points
        !> have been tested
        !
        !>@param cpt_coords
        !> central point whose neighboring grid points should be 
        !> tested
        !
        !>@param nb_mgrdpts
        !> number of grid points to be modified
        !
        !>@param mgrdpts
        !> list of the bc_interior_pt grid points to be modified
        !--------------------------------------------------------------
        subroutine check_neighboring_bc_interior_pts(
     $     this,
     $     cpt_coords_p, cpt_coords,
     $     nb_mgrdpts, mgrdpts)

          implicit none

          class(bf_interface_icr)     , intent(in)  :: this
          integer(ikind), dimension(2), intent(in)  :: cpt_coords_p
          integer(ikind), dimension(2), intent(in)  :: cpt_coords
          integer                     , intent(out) :: nb_mgrdpts
          integer, dimension(:,:)     , intent(out) :: mgrdpts


          type(bf_sublayer), pointer :: sublayer
          integer(ikind), dimension(2) :: l_coords

          !1) analyze the coordinates of the point
          !   check if the point is located inside the interior domain
          !   or if it is located outside the interior and eventually
          !   in a buffer layer

          !if it is inside the border layer, the procedure creating
          !a template neighboring grid point is used
          if(is_inside_border_layer(cpt_coords)) then
             call check_neighboring_bc_interior_pts_for_interior(
     $            this,
     $            cpt_coords_p,
     $            cpt_coords,
     $            nb_mgrdpts,
     $            mgrdpts)
             
          !otherwise, the neighboring grid points to be tested are
          !outside the interior and only a buffer layer can be used
          !to test it
          else
             sublayer => this%get_sublayer(cpt_coords, l_coords)
             if(associated(sublayer)) then
                call sublayer%check_neighboring_bc_interior_pts(
     $               cpt_coords_p(1), cpt_coords_p(2),
     $               cpt_coords(1), cpt_coords(2),
     $               nb_mgrdpts,
     $               mgrdpts)
             end if

          end if
          
        end subroutine check_neighboring_bc_interior_pts    


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the identity of the grid points surrounding a central
        !> point: this function encapsulates the function used if the
        !> central point is at the interface between the interior and
        !> the buffer layers and the function checking the grid points
        !> in case all the grid points are inside a buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param cpt_coords
        !> central point whose neighboring grid points should be 
        !> tested
        !
        !>@return is_border
        !> check whether the central point to be tested is at the
        !> interface between the interiro domain and the buffer layers
        !--------------------------------------------------------------
        function is_inside_border_layer(cpt_coords) result(is_border)

          implicit none

          integer(ikind), dimension(2), intent(in) :: cpt_coords
          logical                                  :: is_border


          logical :: is_interior_and_border
          logical :: is_interior


          is_interior_and_border =
     $         ((cpt_coords(1)).ge.1).and.
     $         ((cpt_coords(1)).le.nx).and.
     $         ((cpt_coords(2)).ge.1).and.
     $         ((cpt_coords(2)).le.ny)

          is_interior =
     $         (cpt_coords(1).ge.(bc_size+2)).and.
     $         (cpt_coords(1).le.(nx-bc_size-1)).and.
     $         (cpt_coords(2).ge.(bc_size+2)).and.
     $         (cpt_coords(2).le.(ny-bc_size-1))
          
          is_border = is_interior_and_border.and.(.not.is_interior)

        end function is_inside_border_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the neighboring bc_interior_pt around the central
        !> point identified by its general coordinates
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param cpt_coords_p
        !> previous central point whose neighboring grid points
        !> have been tested
        !
        !>@param cpt_coords
        !> central point whose neighboring grid points should be 
        !> tested
        !
        !>@param nb_mgrdpts
        !> number of grid points to be modified
        !
        !>@param mgrdpts
        !> list of the bc_interior_pt grid points to be modified        
        !--------------------------------------------------------------
        subroutine check_neighboring_bc_interior_pts_for_interior(
     $     this,
     $     cpt_coords_p,
     $     cpt_coords,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          class(bf_interface_icr)       , intent(in)    :: this
          integer(ikind), dimension(2)  , intent(in)    :: cpt_coords_p
          integer(ikind), dimension(2)  , intent(in)    :: cpt_coords
          integer                       , intent(inout) :: nb_mgrdpts
          integer(ikind), dimension(:,:), intent(out)   :: mgrdpts


          integer, dimension(3,3) :: nbc_template


          !1) initialize an array containing the grid points
          !   identity around the central point and update this
          !   array using the potential buffer layers overlapping
          !   this array
          nbc_template = create_nbc_interior_pt_template(this, cpt_coords)


          !2) identify which grid points are bc_interior_pt
          !   and save them in the mgrdpts array with their general
          !   coordinates
          call check_nbc_interior_pt_template(
     $         nbc_template,
     $         cpt_coords_p(1), cpt_coords_p(2),
     $         cpt_coords(1), cpt_coords(2),
     $         nb_mgrdpts,
     $         mgrdpts)

        end subroutine check_neighboring_bc_interior_pts_for_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the grid points around the center point are
        !> bc_interior_pt
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param nbc_template
        !> temporary array created to combine information on the
        !> role of the grid points at the interface between the interior
        !> domain and the buffer layers
        !
        !>@param i_prev
        !> x-index of the previous central point whose neighboring grid
        !> points have been tested
        !
        !>@param j_prev
        !> y-index of the previous central point whose neighboring grid
        !> points have been tested
        !
        !>@param i_center
        !> x-index of the central point whose neighboring grid
        !> points should be tested
        !
        !>@param j_center
        !> y-index of the central point whose neighboring grid
        !> points should be tested
        !
        !>@param nb_mgrdpts
        !> number of grid points to be modified
        !
        !>@param mgrdpts
        !> list of the bc_interior_pt grid points to be modified        
        !--------------------------------------------------------------
        subroutine check_nbc_interior_pt_template(
     $     nbc_template,
     $     i_prev, j_prev,
     $     i_center, j_center,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          integer, dimension(3,3)        , intent(in)    :: nbc_template
          integer(ikind)                 , intent(in)    :: i_prev
          integer(ikind)                 , intent(in)    :: j_prev
          integer(ikind)                 , intent(in)    :: i_center
          integer(ikind)                 , intent(in)    :: j_center
          integer                        , intent(inout) :: nb_mgrdpts
          integer(ikind) , dimension(:,:), intent(out)   :: mgrdpts


          !radius for the search of bc_interior_pt around the
          !central point identified by (i_center, j_center)
          integer, parameter :: search_r = 1

          integer(ikind), dimension(2) :: match_table
          integer(ikind) :: min_j, max_j
          integer(ikind) :: size_x, size_y
          integer(ikind) :: i,j


          !get the match table converting the general coords
          !into local coords
          match_table = [i_center-2, j_center-2]

          !get the borders of the loops
          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          size_x = size(nbc_template,1)
          size_y = size(nbc_template,2)


          do j=max(1,-search_r+j_center-match_table(2)),
     $         min(size_y, j_prev-search_r-1-match_table(2))

             do i=max(1,-search_r+i_center-match_table(1)),
     $            min(size_x, i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               nbc_template,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do


          do j=max(1,-search_r+j_center+min_j-match_table(2)),
     $         min(size_y, j_center+search_r-max_j-match_table(2))

             do i=max(1,-search_r+i_center-match_table(1)),
     $            min(size_x,i_prev-search_r-1-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               nbc_template,
     $               nb_mgrdpts,
     $               mgrdpts)

             end do
          end do


          do j=max(1,j_center-search_r-min_j-match_table(2)),
     $         min(size_y,j_center+search_r-max_j-match_table(2))

             do i=max(1,i_prev+search_r+1-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               nbc_template,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do


          do j=max(1,j_prev+search_r+1-match_table(2)),
     $         min(size_y,j_center+search_r-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               nbc_template,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do

        end subroutine check_nbc_interior_pt_template


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid point tested is a bc_interior_pt
        !> and if so save the general coordinates of the grid point
        !> in mgrdpts
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param i
        !> x-index of the grid point tested
        !
        !>@param j
        !> y-index of the grid point tested
        !
        !>@param match_table
        !> table converting the general coordinates into local coordinates
        !
        !>@param nbc_template
        !> temporary array created to combine information on the
        !> role of the grid points at the interface between the interior
        !> domain and the buffer layers
        !
        !>@param nb_mgrdpts
        !> number of grid points to be modified
        !
        !>@param mgrdpts
        !> list of the bc_interior_pt grid points to be modified        
        !--------------------------------------------------------------
        subroutine check_bc_interior_pt(
     $     i,j,
     $     match_table,
     $     nbc_template,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          integer(ikind)                , intent(in)    :: i,j
          integer(ikind), dimension(2)  , intent(in)    :: match_table
          integer(ikind), dimension(3,3), intent(in)    :: nbc_template
          integer                       , intent(inout) :: nb_mgrdpts
          integer(ikind), dimension(:,:), intent(out)   :: mgrdpts

          if(nbc_template(i,j).eq.bc_interior_pt) then

             nb_mgrdpts = nb_mgrdpts+1
             mgrdpts(1,nb_mgrdpts) = i+match_table(1)
             mgrdpts(2,nb_mgrdpts) = j+match_table(2)
             
          end if

        end subroutine check_bc_interior_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create the template for the neighboring points around
        !> the central point asked: as the central is located in
        !> a layer at the interface between the interior points
        !> and the boundary layers, it is required to initialize
        !> the template using the coordinates as if there were no
        !> boundary layers, then data are exchanged with the
        !> neighboring buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param cpt_coords
        !> central point whose neighboring grid points should
        !> be tested
        !
        !>@return nbc_template
        !> temporary array created to combine information on the
        !> role of the grid points at the interface between the interior
        !> domain and the buffer layers
        !--------------------------------------------------------------
        function create_nbc_interior_pt_template(this, cpt_coords)
     $     result(nbc_template)

          implicit none

          class(bf_interface_icr)       , intent(in) :: this
          integer(ikind), dimension(2)  , intent(in) :: cpt_coords
          integer       , dimension(3,3)             :: nbc_template

          integer, parameter :: i_lim0  = 1
          integer, parameter :: i_lim1  = bc_size
          integer, parameter :: i_lim2  = bc_size+1
          integer, parameter :: i_lim31 = bc_size+3
          integer, parameter :: i_lim32 = nx-bc_size-1
          integer, parameter :: i_lim4  = nx-bc_size
          integer, parameter :: i_lim5  = nx-bc_size+1
          integer, parameter :: i_lim6  = nx

          integer, parameter :: j_lim0  = 1
          integer, parameter :: j_lim1  = bc_size
          integer, parameter :: j_lim2  = bc_size+1
          integer, parameter :: j_lim31 = bc_size+3
          integer, parameter :: j_lim32 = ny-bc_size-1
          integer, parameter :: j_lim4  = ny-bc_size
          integer, parameter :: j_lim5  = ny-bc_size+1
          integer, parameter :: j_lim6  = ny

          
          select case(cpt_coords(1))

            case(i_lim0)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_00()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)
                 case(j_lim1)
                    nbc_template = make_nbc_template_01()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2, j_lim31)
                    nbc_template = make_nbc_template_02()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32, j_lim4)
                    nbc_template = make_nbc_template_02()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_05()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_06()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_02()
               end select

               call update_grdpts_id_for_template(
     $              this, W, cpt_coords, nbc_template)

            case(i_lim1)

               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_10()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)
                 case(j_lim1)
                    nbc_template = make_nbc_template_11()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_12()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_13()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_13()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_14()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_15()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_16()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_13()
               end select

               call update_grdpts_id_for_template(
     $              this, W, cpt_coords, nbc_template)

            case(i_lim2)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_20()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_21()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_22()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_23()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_23()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_24()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_25()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_26()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_23()
               end select

               call update_grdpts_id_for_template(
     $              this, W, cpt_coords, nbc_template)

            case(i_lim31)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_20()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_31()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_32()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_34()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_35()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_26()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)
               end select

               call update_grdpts_id_for_template(
     $              this, W, cpt_coords, nbc_template)

            case(i_lim32)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_20()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_31()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_32()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_34()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_35()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_26()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)
               end select

               call update_grdpts_id_for_template(
     $              this, E, cpt_coords, nbc_template)

            case(i_lim4)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_20()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_41()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_42()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_43()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_43()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_44()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_45()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_26()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_43()
               end select

               call update_grdpts_id_for_template(
     $              this, E, cpt_coords, nbc_template)

            case(i_lim5)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_50()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_51()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_52()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_53()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_53()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_54()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_55()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_56()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_53()
               end select

               call update_grdpts_id_for_template(
     $              this, E, cpt_coords, nbc_template)

            case(i_lim6)
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_60()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_61()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_62()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_62()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_62()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_62()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_65()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_66()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case default
                    nbc_template = make_nbc_template_62()
               end select

               call update_grdpts_id_for_template(
     $              this, E, cpt_coords, nbc_template)

            case default
               select case(cpt_coords(2))
                 case(j_lim0)
                    nbc_template = make_nbc_template_20()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim1)
                    nbc_template = make_nbc_template_31()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim2)
                    nbc_template = make_nbc_template_32()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim31)
                    nbc_template = make_nbc_template_32()
                    call update_grdpts_id_for_template(
     $                   this, S, cpt_coords, nbc_template)

                 case(j_lim32)
                    nbc_template = make_nbc_template_32()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim4)
                    nbc_template = make_nbc_template_34()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim5)
                    nbc_template = make_nbc_template_35()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)

                 case(j_lim6)
                    nbc_template = make_nbc_template_26()
                    call update_grdpts_id_for_template(
     $                   this, N, cpt_coords, nbc_template)
               end select

          end select

        end function create_nbc_interior_pt_template


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the content of the grdpts_id table from neighboring
        !> buffer layers to the nbc_template
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param mainlayer_id
        !> cardinal coordinate locating where the neighboring
        !> buffer layers are located
        !
        !>@param cpt_coords
        !> central point whose neighboring grid points should
        !> be tested
        !
        !>@param nbc_template
        !> temporary array created to combine information on the
        !> role of the grid points at the interface between the interior
        !> domain and the buffer layers
        !--------------------------------------------------------------
        subroutine update_grdpts_id_for_template(
     $     this, mainlayer_id, cpt_coords, nbc_template)

          implicit none

          class(bf_interface_icr)       , intent(in) :: this
          integer                       , intent(in) :: mainlayer_id
          integer(ikind), dimension(2)  , intent(in) :: cpt_coords
          integer       , dimension(3,3), intent(out):: nbc_template


          type(bf_sublayer), pointer   :: sublayer
          integer(ikind), dimension(2) :: local_coords


          sublayer => this%get_sublayer(
     $         cpt_coords,
     $         local_coords,
     $         tolerance_i=1,
     $         mainlayer_id_i=mainlayer_id)

          if(associated(sublayer)) then
             call sublayer%copy_grdpts_id_to_temp(
     $            cpt_coords, nbc_template)
          end if

        end subroutine update_grdpts_id_for_template


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove a bf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param sublayer_ptr
        !> reference to the bf_sublayer which is removed
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate locating the buffer layer removed
        !--------------------------------------------------------------
        subroutine remove_sublayer(this, sublayer_ptr, bf_mainlayer_id)

          implicit none
          
          class(bf_interface_icr)   , intent(inout) :: this
          type(bf_sublayer), pointer, intent(inout) :: sublayer_ptr
          integer         , optional, intent(in)    :: bf_mainlayer_id
          
          integer :: mainlayer_id


          !cardinal coordinate of the sublayer removed
          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = sublayer_ptr%get_localization()
          end if


          !reoganize the increasing detectors belonging
          !to the sublayer removed
          call update_icr_detectors_after_removal(
     $         this,
     $         sublayer_ptr%get_alignment_tab(),
     $         mainlayer_id)
          
          !remove the sublayer
          call this%bf_interface%remove_sublayer(
     $         sublayer_ptr,
     $         mainlayer_id)

        end subroutine remove_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the position of the increasing detectors
        !> if a sublayer is removed
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param bf_align
        !> alignment of the bf_sublayer which is removed
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate locating the buffer layer removed
        !--------------------------------------------------------------
        subroutine update_icr_detectors_after_removal(
     $     this, bf_align, bf_mainlayer_id)

          class(bf_interface_icr)       , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: bf_align
          integer                       , intent(in)    :: bf_mainlayer_id

          type(bf_detector_dcr_list_N) :: dcr_param_N
          type(bf_detector_dcr_list_S) :: dcr_param_S
          type(bf_detector_dcr_list_E) :: dcr_param_E
          type(bf_detector_dcr_list_W) :: dcr_param_W


          select case(bf_mainlayer_id)
            case(N)
               call update_icr_detectors(this, bf_align, dcr_param_N)
            case(S)
               call update_icr_detectors(this, bf_align, dcr_param_S)
            case(E)
               call update_icr_detectors(this, bf_align, dcr_param_E)
            case(W)
               call update_icr_detectors(this, bf_align, dcr_param_W)
            case default
               call error_mainlayer_id(
     $              'bf_interface_icr_class.f',
     $              'update_icr_detectors_after_removal',
     $              bf_mainlayer_id)
          end select

        end subroutine update_icr_detectors_after_removal

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the detectors by using a temporary object
        !> where the parameters for the changes are stored
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param bf_align
        !> alignment of the bf_sublayer which is removed
        !
        !>@param dcr_param
        !> temporary object used to identify the detectors that should
        !> be removed from the list
        !--------------------------------------------------------------
        subroutine update_icr_detectors(this, bf_align, dcr_param)

          implicit none

          class(bf_interface_icr)       , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: bf_align
          class(bf_detector_dcr_list)   , intent(in)    :: dcr_param

          class(bf_detector_dcr_list), allocatable :: N_dcr_param
          class(bf_detector_dcr_list), allocatable :: S_dcr_param
          class(bf_detector_dcr_list), allocatable :: E_dcr_param
          class(bf_detector_dcr_list), allocatable :: W_dcr_param

          integer(ikind), dimension(2) :: first_pt_linked
          integer(ikind), dimension(2) :: last_pt_linked
          
          !create the temporary objects saving the parameters
          !for the update of the detector lists due to the removal
          !of a sublayer
          allocate(N_dcr_param, source=dcr_param)
          allocate(S_dcr_param, source=dcr_param)
          allocate(E_dcr_param, source=dcr_param)
          allocate(W_dcr_param, source=dcr_param)


          !initialize the objects saving the parameters
          !when constructing the new detector lists
          call N_dcr_param%ini()
          call S_dcr_param%ini()
          call E_dcr_param%ini()
          call W_dcr_param%ini()

          
          !compute the parameters for the construction
          !of the new detector lists
          call N_dcr_param%compute_new_list_param(bf_align, this%N_detectors_list)
          call S_dcr_param%compute_new_list_param(bf_align, this%S_detectors_list)
          call E_dcr_param%compute_new_list_param(bf_align, this%E_detectors_list)
          call W_dcr_param%compute_new_list_param(bf_align, this%W_detectors_list)
          

          !check the overlap of detectors between the lists
          call check_detector_overlap(
     $         N_dcr_param, S_dcr_param, E_dcr_param, W_dcr_param)


          !compute the new detector lists and link the lists
          first_pt_linked = W_dcr_param%get_last_detector()
          last_pt_linked  = E_dcr_param%get_last_detector()
          call N_dcr_param%compute_new_list(
     $         this%N_detectors_list,
     $         first_pt_linked,
     $         last_pt_linked)

          first_pt_linked = W_dcr_param%get_first_detector()
          last_pt_linked  = E_dcr_param%get_first_detector()
          call S_dcr_param%compute_new_list(
     $         this%S_detectors_list,
     $         first_pt_linked,
     $         last_pt_linked)

          first_pt_linked = S_dcr_param%get_last_detector()
          last_pt_linked  = N_dcr_param%get_last_detector()
          call E_dcr_param%compute_new_list(
     $         this%E_detectors_list,
     $         first_pt_linked,
     $         last_pt_linked)

          first_pt_linked = S_dcr_param%get_first_detector()
          last_pt_linked  = N_dcr_param%get_first_detector()
          call W_dcr_param%compute_new_list(
     $         this%W_detectors_list,
     $         first_pt_linked,
     $         last_pt_linked)


          !remove the temporary objects
          deallocate(N_dcr_param)
          deallocate(S_dcr_param)
          deallocate(E_dcr_param)
          deallocate(W_dcr_param)

        end subroutine update_icr_detectors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< check if detectors overlap and modify the
        !> the borders for the detectors
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param N_dcr_param
        !> temporary object with the parameters for updating the
        !> increasing north detectors after the removal of a buffer
        !> layer
        !
        !>@param S_dcr_param
        !> temporary object with the parameters for updating the
        !> increasing south detectors after the removal of a buffer
        !> layer
        !
        !>@param E_dcr_param
        !> temporary object with the parameters for updating the
        !> increasing east detectors after the removal of a buffer
        !> layer
        !
        !>@param W_dcr_param
        !> temporary object with the parameters for updating the
        !> increasing west detectors after the removal of a buffer
        !> layer
        !--------------------------------------------------------------
        subroutine check_detector_overlap(
     $     N_dcr_param, S_dcr_param, E_dcr_param, W_dcr_param)

          implicit none

          class(bf_detector_dcr_list), intent(inout) :: N_dcr_param
          class(bf_detector_dcr_list), intent(inout) :: S_dcr_param
          class(bf_detector_dcr_list), intent(inout) :: E_dcr_param
          class(bf_detector_dcr_list), intent(inout) :: W_dcr_param        
        
          integer(ikind), dimension(2) :: checked_pt1
          integer(ikind), dimension(2) :: checked_pt2

          !south - west overlap
          checked_pt1 = S_dcr_param%get_first_detector()
          checked_pt2 = W_dcr_param%get_first_detector()
          if(overlap(checked_pt1,checked_pt2)) then
             checked_pt2(2) = checked_pt2(2)+1
             call W_dcr_param%set_first_detector(checked_pt2)
          end if

          !south - east overlap
          checked_pt1 = S_dcr_param%get_last_detector()
          checked_pt2 = E_dcr_param%get_first_detector()
          if(overlap(checked_pt1,checked_pt2)) then
             checked_pt2(2) = checked_pt2(2)+1
             call E_dcr_param%set_first_detector(checked_pt2)
          end if

          !north - west overlap
          checked_pt1 = N_dcr_param%get_first_detector()
          checked_pt2 = W_dcr_param%get_last_detector()
          if(overlap(checked_pt1,checked_pt2)) then
             checked_pt2(2) = checked_pt2(2)-1
             call W_dcr_param%set_last_detector(checked_pt2)
          end if

          !north - east overlap
          checked_pt1 = N_dcr_param%get_first_detector()
          checked_pt2 = E_dcr_param%get_last_detector()
          if(overlap(checked_pt1,checked_pt2)) then
             checked_pt2(2) = checked_pt2(2)-1
             call E_dcr_param%set_last_detector(checked_pt2)
          end if

        end subroutine check_detector_overlap


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether two grid points are the same
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param checked_pt1
        !> general coordinates of the first grid point checked
        !
        !>@param checked_pt2
        !> general coordinates of the second grid point checked
        !
        !>@return overlap
        !> logical stating whether the two grid points are the same
        !--------------------------------------------------------------
        function overlap(checked_pt1, checked_pt2)

          implicit none

          integer(ikind), dimension(2), intent(in) :: checked_pt1
          integer(ikind), dimension(2), intent(in) :: checked_pt2
          logical                                  :: overlap

          overlap = (checked_pt1(1).eq.checked_pt2(1)).and.
     $              (checked_pt1(2).eq.checked_pt2(2))
          
        end function overlap


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the increasing detector positions on a matrix
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param matrix
        !> array on which the position of the increasing detectors
        !> is indicated
        !--------------------------------------------------------------
        subroutine print_idetectors_on(this, matrix)

          implicit none

          class(bf_interface_icr)    , intent(in)  :: this
          real(rkind), dimension(:,:), intent(out) :: matrix

          real(rkind), parameter :: N_detector_color = 0.9d0
          real(rkind), parameter :: S_detector_color = 0.7d0
          real(rkind), parameter :: E_detector_color = 0.5d0
          real(rkind), parameter :: W_detector_color = 0.2d0
          integer(ikind) :: i,j,k


          if(allocated(this%N_detectors_list)) then
             do k=1, size(this%N_detectors_list,2)
                i = this%N_detectors_list(1,k)
                j = this%N_detectors_list(2,k)
                matrix(i,j) = N_detector_color
             end do
          end if

          if(allocated(this%S_detectors_list)) then
             do k=1, size(this%S_detectors_list,2)
                i = this%S_detectors_list(1,k)
                j = this%S_detectors_list(2,k)
                matrix(i,j) = S_detector_color
             end do
          end if

          if(allocated(this%E_detectors_list)) then
             do k=1, size(this%E_detectors_list,2)
                i = this%E_detectors_list(1,k)
                j = this%E_detectors_list(2,k)
                matrix(i,j) = E_detector_color
             end do
          end if
          
          if(allocated(this%W_detectors_list)) then
             do k=1, size(this%W_detectors_list,2)
                i = this%W_detectors_list(1,k)
                j = this%W_detectors_list(2,k)
                matrix(i,j) = W_detector_color
             end do
          end if

        end subroutine print_idetectors_on


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the increasing detector positions on binary
        !> output files
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_icr object encapsulating the position of
        !> the increasing detectors and the subroutine controlling
        !> the extension of the computational domain
        !
        !>@param index
        !> integer identifying the file
        !--------------------------------------------------------------
        subroutine print_idetectors_on_binary(this, index)

          implicit none

          class(bf_interface_icr), intent(in) :: this
          integer                , intent(in) :: index

          character(len=11) :: filename_format
          character(len=20) :: filename
          integer :: ios
          integer :: index_format
          

          if(index.le.9) then
             index_format = 1
          else
             if(index.le.99) then
                index_format = 2
             else
                index_format = 3
             end if
          end if

          write(filename_format, '(''(A11,I'',I1,'',A4)'')') index_format


          !N detectors
          write(filename, filename_format) 'N_detectors', index, '.dat'

          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) this%N_detectors_list
              close(unit=1)
           else
              stop 'file opening pb'
           end if


          !S detectors
          write(filename, filename_format) 'S_detectors', index, '.dat'

          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) this%S_detectors_list
              close(unit=1)
           else
              stop 'file opening pb'
           end if


          !E detectors
          write(filename, filename_format) 'E_detectors', index, '.dat'

          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) this%E_detectors_list
              close(unit=1)
           else
              stop 'file opening pb'
           end if


          !W detectors
          write(filename, filename_format) 'W_detectors', index, '.dat'

          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) this%W_detectors_list
              close(unit=1)
           else
              stop 'file opening pb'
           end if


        end subroutine print_idetectors_on_binary

      end module bf_interface_icr_class
