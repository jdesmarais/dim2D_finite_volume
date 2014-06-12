      module bf_interface_icr_class

        use bf_detector_i_list_class, only : bf_detector_i_list
        use bf_path_icr_class     , only : bf_path_icr
        use bf_nbc_template_module
        use bf_sublayer_class       , only : bf_sublayer
        use bf_interface_class      , only : bf_interface
        use parameters_bf_layer     , only : bc_interior_pt
        use parameters_constant     , only : N,S,E,W,interior
        use parameters_input        , only : nx,ny,ne,bc_size,
     $                                       dt,search_nb_dt
        use parameters_kind         , only : ikind, rkind

        implicit none

        private
        public :: bf_interface_icr


        type, extends(bf_interface) :: bf_interface_icr

          integer(ikind), dimension(:,:), allocatable, private :: N_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: S_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: E_detectors_list
          integer(ikind), dimension(:,:), allocatable, private :: W_detectors_list

          contains

          procedure,   pass :: ini
          procedure,   pass :: get_modified_grdpts_list
          procedure,   pass :: process_idetector_list
          procedure,   pass :: combine_bf_idetector_lists
          procedure,   pass :: update_bf_layers_with_idetectors

          procedure, nopass, private :: is_i_detector_activated
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


        !< initialize the detector lists
        subroutine ini(this)

          implicit none

          class(bf_interface_icr), intent(inout) :: this

          integer(ikind) :: i, dbf_distance


          !initialize the parent attributes
          call this%bf_interface%ini()


          !intialize the attributes specific to bf_interface_icr
          dbf_distance = bc_size

          allocate(this%N_detectors_list(2,nx-2*(bc_size+dbf_distance)+2))
          allocate(this%S_detectors_list(2,nx-2*(bc_size+dbf_distance)+2))
          allocate(this%E_detectors_list(2,ny-2*(bc_size+dbf_distance)))
          allocate(this%W_detectors_list(2,nx-2*(bc_size+dbf_distance)))

          do i=bc_size+dbf_distance, nx-(bc_size+dbf_distance)+1
             this%S_detectors_list(1,i-(bc_size+dbf_distance)+1) = i
             this%S_detectors_list(2,i-(bc_size+dbf_distance)+1) = bc_size+dbf_distance
          end do

          do i=bc_size+dbf_distance, nx-(bc_size+dbf_distance)+1
             this%N_detectors_list(1,i-(bc_size+dbf_distance)+1) = i
             this%N_detectors_list(2,i-(bc_size+dbf_distance)+1) = ny-(bc_size+dbf_distance)+1
          end do

          do i=bc_size+dbf_distance+1, ny-(bc_size+dbf_distance)
             this%W_detectors_list(1,i-(bc_size+dbf_distance)) = bc_size+dbf_distance
             this%W_detectors_list(2,i-(bc_size+dbf_distance)) = i
          end do

          do i=bc_size+dbf_distance+1, ny-(bc_size+dbf_distance)
             this%E_detectors_list(1,i-(bc_size+dbf_distance)) = nx-(bc_size+dbf_distance)+1
             this%E_detectors_list(2,i-(bc_size+dbf_distance)) = i
          end do          

        end subroutine ini


        !< update the size of the buffer layers using the information 
        !> given by the increasing detectors
        subroutine update_bf_layers_with_idetectors(this, 
     $     interior_nodes, dx, dy)

          implicit none

          class(bf_interface_icr)         , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy

          type(bf_path_icr)          :: path_update_idetectors
          integer(ikind), dimension(2) :: cpt_coords_p
          type(bf_detector_i_list)     :: N_ndt_list
          type(bf_detector_i_list)     :: S_ndt_list
          type(bf_detector_i_list)     :: E_ndt_list
          type(bf_detector_i_list)     :: W_ndt_list


          !initialization of the path that will gather information
          !on the buffer layers to be updated: the path collects
          !the updating information belonging to a specific buffer
          !layer and when another buffer layer is activated by
          !the detectors the path applies all the changes to the
          !previous buffer layer and starts collecting information
          !on the next buffer layer
          call path_update_idetectors%ini()

          !there are 4 layers of detectors cooresponding to the 4
          !cardinal points. They create a closed path of detectors
          !we update in this order:
          !1: S, 2: E, 3: W, 4: N
          !in this way the path created by one buffer layer can
          !be continued by the next list of detectors and prevent
          !double operations on buffer layers
          
          !0) initialization of the first point indicating 
          !   where the previous neighboring bc_interior-pt
          !   were analyzed
          cpt_coords_p = [nx/2, ny/2]

          !1) South detectors
          if(allocated(this%S_detectors_list)) then
             call S_ndt_list%ini(S, size(this%S_detectors_list,2))

             call process_idetector_list(
     $            this, this%S_detectors_list, S_ndt_list,
     $            interior_nodes, dx, dy,
     $            cpt_coords_p, path_update_idetectors)

          end if

          !2) East detectors
          if(allocated(this%E_detectors_list)) then
             call E_ndt_list%ini(E, size(this%E_detectors_list,2))

             call process_idetector_list(
     $            this, this%E_detectors_list, E_ndt_list,
     $            interior_nodes, dx, dy,
     $            cpt_coords_p, path_update_idetectors)

          end if

          !3) West detectors
          if(allocated(this%W_detectors_list)) then
             call W_ndt_list%ini(W, size(this%W_detectors_list,2))

             call process_idetector_list(
     $            this, this%W_detectors_list, W_ndt_list,
     $            interior_nodes, dx, dy,
     $            cpt_coords_p, path_update_idetectors)

          end if

          !4) North detectors
          if(allocated(this%N_detectors_list)) then
             call N_ndt_list%ini(N, size(this%N_detectors_list,2))

             call process_idetector_list(
     $            this, this%N_detectors_list, N_ndt_list,
     $            interior_nodes, dx, dy,
     $            cpt_coords_p, path_update_idetectors)

          end if

          !5) process the last path as after the detectors of
          !   the last main layer are analysed, it is possible
          !   that the final part of the path has not been
          !   processed
          if(path_update_idetectors%get_nb_pts().gt.0) then
             call path_update_idetectors%process_path(this, interior_nodes)
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


        !< connect the bf_detector_i_list objects and determine the
        !> new detector lists
        subroutine combine_bf_idetector_lists(
     $     this,
     $     N_idetectors_list, S_idetectors_list,
     $     E_idetectors_list, W_idetectors_list)

          implicit none

          class(bf_interface_icr) , intent(inout) :: this
          type(bf_detector_i_list), intent(in)    :: N_idetectors_list
          type(bf_detector_i_list), intent(in)    :: S_idetectors_list
          type(bf_detector_i_list), intent(in)    :: E_idetectors_list
          type(bf_detector_i_list), intent(in)    :: W_idetectors_list


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


        !< process the list of current detectors, modify the buffer layers
        !> if they are activated by the detectors and determine the new list
        !> of detectors
        !> dt_list  : detector list
        !> ndt_list : new detector list
        subroutine process_idetector_list(
     $     this, dt_list, ndt_list,
     $     interior_nodes, dx, dy,
     $     cpt_coords_p, path)
        
          implicit none

          class(bf_interface_icr)                      , intent(inout) :: this
          integer(ikind)          , dimension(:,:)     , intent(in)    :: dt_list
          type(bf_detector_i_list)                     , intent(inout) :: ndt_list
          real(rkind)             , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                                  , intent(in)    :: dx
          real(rkind)                                  , intent(in)    :: dy
          integer(ikind)          , dimension(2)       , intent(inout) :: cpt_coords_p
          type(bf_path_icr)                          , intent(inout) :: path

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
     $            interior_nodes,
     $            dx, dy,
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
                   call path%process_path(this, interior_nodes)

                   !the grid point that led to the path end
                   !is used to reinitialize the current path
                   call path%analyze_pt(mgrdpts(:,l), this)

                end if

             end do
          end do

        end subroutine process_idetector_list


        !< from a detector position, get the bc_interior_pt activated
        !> d_coords     : detector general coordinates
        !> cpt_coords_p : previous grid point whose neighboring grid points
        !>                were checked to find bc_interior_pt
        !> nb_mgrdpts   : number of grid points to be modified
        !> mgrdpts      : table containing the coordinates of the grid
        !>                points to be modified
        !> ndt_list     : temporary new detector list
        subroutine get_modified_grdpts_list(
     $     this, d_coords,
     $     interior_nodes,
     $     dx, dy,
     $     cpt_coords_p, cpt_coords,
     $     nb_mgrdpts, mgrdpts, ndt_list)

          implicit none

          class(bf_interface_icr)         , intent(inout) :: this
          integer(ikind), dimension(2)    , intent(in)    :: d_coords
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          integer(ikind), dimension(2)    , intent(in)    :: cpt_coords_p
          integer(ikind), dimension(2)    , intent(out)   :: cpt_coords
          integer                         , intent(out)   :: nb_mgrdpts
          integer(ikind), dimension(2,9)  , intent(out)   :: mgrdpts
          type(bf_detector_i_list)        , intent(inout) :: ndt_list


          real(rkind), dimension(ne)   :: node_var
          real(rkind)   , dimension(2) :: velocity
          integer(ikind), dimension(2) :: d_coords_n


          !initialization of the number of modified grid points
          nb_mgrdpts = 0


          !if the detector is activated, then we check
          !whether grid points need to be modified
          if(is_i_detector_activated()) then

             !extract the velocity at the coordinates of the detector
             node_var = this%get_nodes(d_coords, interior_nodes)
             velocity(1) = node_var(2)/node_var(1)
             velocity(2) = node_var(3)/node_var(1)
             

             !get the first point from which we should look for a
             !bc_interior_pt to be activated and teh new coordinates
             !from the detector
             cpt_coords = get_central_grdpt(
     $            d_coords, velocity, dx, dy, d_coords_n)
             
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


        !> check whether the detector is activated
        function is_i_detector_activated()
        
          implicit none
          
          logical :: is_i_detector_activated
          
          is_i_detector_activated = .true.

        end function is_i_detector_activated



        !> get the general coordinates of the point activated by an
        !> increasing detector
        function get_central_grdpt(
     $     d_coords, velocity, dx, dy, d_coords_n)
     $     result(cpt_coords)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: d_coords
          real(rkind)   , dimension(2), intent(in)  :: velocity
          real(rkind)                 , intent(in)  :: dx
          real(rkind)                 , intent(in)  :: dy
          integer(ikind), dimension(2), intent(out) :: d_coords_n
          integer(ikind), dimension(2)              :: cpt_coords

          real(rkind) :: dir_x, dir_y

          !1) get the direction to look for a bc_interior_pt
          dir_x  = velocity(1)*search_nb_dt*dt/dx
          dir_y  = velocity(2)*search_nb_dt*dt/dy
          
          !2) get the point indices in the direction given
          !   by the velocity vector
          cpt_coords(1) = d_coords(1) + nint(dir_x)
          cpt_coords(2) = d_coords(2) + nint(dir_y)
          
          !3) compute the new detector position
          d_coords_n(1) = d_coords(1) + nint(max(min(dir_x,1.0d0),-1.0d0)) !nint(dir_x/SQRT(dir_x**2+dir_y**2))
          d_coords_n(2) = d_coords(2) + nint(max(min(dir_y,1.0d0),-1.0d0)) !nint(dir_y/SQRT(dir_x**2+dir_y**2))

        end function get_central_grdpt


        !> check the identity of the grid points surrounding a central
        !> point: this function encapsulates the function used if the
        !> central point is at the interface between the interior and
        !> the buffer layers and the function checking the grid points
        !> in case all the grid points are inside a buffer layer
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


        !> check if a grid point is at the interface between the buffer
        !> layers and the interior domain
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


        !< check the neighboring bc_interior_pt around the central
        !> point identified by its general coordinates
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


        !< check if the grid points around the center point are
        !> bc_interior_pt
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


        !< check whether the grid point tested is a bc_interior_pt
        !> and if so save the general coordinates of the grid point
        !> in mgrdpts
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


        !< create the template for the neighboring points around
        !> the central point asked: as the central is located in
        !> a layer at the interface between the interior points
        !> and the boundary layers, it is required to initialize
        !> the template using the coordinates as if there were no
        !> boundary layers, then data are exchanged with the
        !> neighboring buffer layers
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


        !< copy the content of the grdpts_id table from neighboring
        !> buffer layers to the nbc_template
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


        !< print the detector positions on a matrix
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


        !< print the position of the increasing detectors
        !> on external binary files
        subroutine print_idetectors_on_binary(this, index)

          implicit none

          class(bf_interface_icr), intent(in) :: this
          integer                , intent(in) :: index

          character(len=11) :: filename_format
          character(len=20) :: filename
          integer :: ios
          

          write(filename_format, '(''(A11,I'',I1,'',A4)'')') floor(index/10.0)+1


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
