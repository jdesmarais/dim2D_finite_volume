      !> @file
      !> object gathering the increase operations to be applied on
      !> the domain extension
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> object gathering the increase operations to be applied on
      !> the domain extension
      !
      !> @date
      !> 23_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_interface_class

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_sublayer_class, only :
     $       bf_sublayer

        use icr_path_chain_class, only :
     $       icr_path_chain

        use icr_path_list_class, only :
     $       icr_path_list

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       N,S,E,W
        
        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: icr_interface


        !> @class icr_interface
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !> @param paths
        !> objects gathering the update operations to be applied on
        !> each buffer layer
        !
        !> @param ini
        !> initialize the paths
        !
        !> @param stage
        !> stage the grid-point for the update of the buffer layers
        !
        !> @param finalize_domain_increase
        !> commit the non-empty paths for domain update
        !---------------------------------------------------------------
        type :: icr_interface

          type(icr_path_list)           :: paths
          type(icr_path_chain), pointer :: current_path

          contains

          procedure, pass :: ini
          procedure, pass :: stage
          procedure, pass :: commit
          procedure, pass :: finalize_domain_increase          

          procedure, pass :: check_path_before_commit
          procedure, pass :: check_paths_before_commit

        end type icr_interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the paths
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(icr_interface), intent(inout) :: this

          call this%paths%ini()
          nullify(this%current_path)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> stage the grid-point for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !>@param gen_coords
        !> general coordinates of the bc_interior_pt analyzed
        !--------------------------------------------------------------
        subroutine stage(
     $     this,
     $     gen_coords,
     $     bf_interface_used)

          implicit none

          class(icr_interface)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: gen_coords
          class(bf_interface_coords)  , intent(in)    :: bf_interface_used

          type(icr_path_chain), pointer :: path_checked
          logical :: create_new_path
          integer :: k


          !1) stage the grid-point activated identified by gen_coords
          !   in the current path or create a new path to stage the
          !   grid-point if there is no path in this%paths
          if(this%paths%get_nb_paths().eq.0) then
             this%current_path => this%paths%add_path()
          end if
          call this%current_path%stage(gen_coords,bf_interface_used)


          !2) if the grid point cannot be staged in the current path,
          !   we check whether it can be stage in another path stored
          !   in this%paths
          if(this%current_path%is_ended()) then

             create_new_path = .true.

             !check whether a new path should be created to store the
             !activated grid-point or find an existing  path where it
             !can be saved
             path_checked => this%paths%get_head_path()
             do k=1, this%paths%get_nb_paths()
                
                call path_checked%stage(gen_coords,bf_interface_used)
                if(.not.path_checked%is_ended()) then
                   create_new_path = .false.
                   this%current_path => path_checked
                   exit
                end if

                path_checked => path_checked%get_next()

             end do

             !create a new path to store the activated grid-point
             if(create_new_path) then
                this%current_path => this%paths%add_path()
                call this%current_path%stage(gen_coords,bf_interface_used)
             end if
             
          end if

        end subroutine stage


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the current path should be merged with other
        !> paths of this%paths
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !>@param path_checked
        !> pointer to the path checked
        !--------------------------------------------------------------
        subroutine check_path_before_commit(this,path_checked)

          implicit none

          class(icr_interface)         , intent(inout) :: this
          type(icr_path_chain), pointer, intent(inout) :: path_checked

          type(bf_sublayer)   , pointer :: merge_sublayer
          type(icr_path_chain), pointer :: current_path
          type(icr_path_chain), pointer :: path_removed


          !check whether the current path lead to the
          !merge of buffer layers
          if(path_checked%lead_to_merge_sublayers()) then


             !get the buffer layer with which the current buffer
             !layer will be merged
             merge_sublayer => path_checked%get_matching_sublayer()
             merge_sublayer => merge_sublayer%get_next()


             !if another path stored in this%paths as a matching
             !sublayer pointing to merge_sublayer, the path
             !should be merged with the current path
             current_path => this%paths%get_head_path()
             do while(associated(current_path))

                if(associated(current_path%get_matching_sublayer(),
     $                        merge_sublayer)) then

                   call path_checked%merge(current_path,check_for_merge=.false.)

                   if(associated(current_path%get_next())) then
                      path_removed => current_path
                      current_path => current_path%get_next()
                      call this%paths%remove_path(path_removed)

                   else
                      path_removed => current_path
                      nullify(current_path)
                      call this%paths%remove_path(path_removed)

                   end if

                else
                   current_path => current_path%get_next()
                end if

             end do             

          end if                     

        end subroutine check_path_before_commit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the paths saved in this%paths should be
        !> merged with other paths of this%paths for a specific
        !> mainlayer_id
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !>@param mainlayer_id
        !> cardinal coordinate identfying the buffer main layer
        !> investigated
        !--------------------------------------------------------------
        subroutine check_paths_before_commit(this,mainlayer_id)

          implicit none

          class(icr_interface), intent(inout) :: this
          integer             , intent(in)    :: mainlayer_id

          type(icr_path_chain), pointer :: path_checked


          path_checked => this%paths%get_head_path()

          do while(associated(path_checked))

             if(path_checked%get_mainlayer_id().eq.mainlayer_id) then

                call this%check_path_before_commit(path_checked)

             end if

             path_checked => path_checked%get_next()

          end do

        end subroutine check_paths_before_commit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> commit the grid-points belonging to a main layer for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the buffer main layer updated
        !
        !>@param bf_interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
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
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !--------------------------------------------------------------
        subroutine commit(
     $     this,
     $     mainlayer_id,
     $     bf_interface_used,
     $     p_model,
     $     t,
     $     dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(icr_interface)            , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          class(bf_interface_coords)      , intent(inout) :: bf_interface_used
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1


          type(icr_path_chain), pointer :: path_checked
          integer :: k

          
          call this%check_paths_before_commit(mainlayer_id)
          

          path_checked => this%paths%get_head_path()
          
          do k=1, this%paths%get_nb_paths()

             if(path_checked%get_mainlayer_id().eq.mainlayer_id) then

                call path_checked%commit(
     $               bf_interface_used,
     $               p_model,
     $               t,
     $               dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1)

             end if

             path_checked => path_checked%get_next()

          end do

          this%current_path => this%paths%get_head_path()

        end subroutine commit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> commmit the remaining paths for the update of the
        !> buffer layers
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !--------------------------------------------------------------
        subroutine finalize_domain_increase(
     $     this,
     $     bf_interface_used,
     $     p_model,
     $     t,
     $     dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(icr_interface)               , intent(inout) :: this
          class(bf_interface_coords)         , intent(inout) :: bf_interface_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes1


          integer, dimension(4) :: mainlayer_id
          integer               :: k


          mainlayer_id = [N,S,E,W]

          do k=1,4

             call this%commit(
     $            mainlayer_id(k),
     $            bf_interface_used,
     $            p_model,
     $            t,
     $            dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1)          
             
          end do

          nullify(this%current_path)

        end subroutine finalize_domain_increase

      end module icr_interface_class
