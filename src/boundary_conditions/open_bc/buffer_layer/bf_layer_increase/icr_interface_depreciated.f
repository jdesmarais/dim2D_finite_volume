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
      !> 20_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_interface_class

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use icr_path_class, only :
     $       icr_path

        use pmodel_eq_class, only :
     $       pmodel_eq

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
        !> @param current_path_is_head_path
        !> logical indicating whether the path to consider to store
        !> update procedures is the head_path or the current_path
        !
        !> @param head_path
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !> @param current_path
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !> @param ini
        !> initialize the paths
        !
        !> @param finalize_domain_increase
        !> commit the non-empty paths for domain update
        !---------------------------------------------------------------
        type :: icr_interface

          logical        :: current_path_is_head_path
          type(icr_path) :: head_path
          type(icr_path) :: current_path

          contains

          procedure, pass :: ini
          procedure, pass :: stage
          procedure, pass :: finalize_domain_increase
          

        end type icr_interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the paths
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(icr_interface), intent(inout) :: this

          this%current_path_is_head_path = .true.
          call this%head_path%ini()
          call this%current_path%ini()

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> stage the grid-point for update
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the increase operations to be applied on
        !> the domain extension
        !
        !>@param gen_coords
        !> general coordinates of the bc_interior_pt analyzed
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
        subroutine stage(
     $     this,
     $     gen_coords,
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
          integer(ikind), dimension(2)       , intent(in)    :: gen_coords
          class(bf_interface_coords)         , intent(inout) :: bf_interface_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes1


          !1) if the path used to gather the activated grid points
          !   is the head_path, then the grid-point is staged in
          !   the head_path
          if(this%current_path_is_head_path) then
             
             call this%head_path%stage(gen_coords,bf_interface_used)


             !1.1) if the grid-point can not belong to the head_path
             !     the head_path is ended and the grid-point is staged
             !     in the current_path
             if(this%head_path%is_ended()) then

                this%current_path_is_head_path = .false.

             end if

          end if


          !2) if the path used to gather the activated grid points
          !   is the current_path, thne the grid-point is staged in
          !   the current_path
          if(.not.(this%current_path_is_head_path)) then
             
             call this%current_path%stage(gen_coords,bf_interface_used)


             !2.1) if the grid-point can not belong to the current_path
             !     the current_path is ended and the grid-point is staged
             !     once the current_path has been commit
             if(this%current_path%is_ended()) then

                
                !2.1.1) before committing the current_path, we need to
                !       check whether both paths can be merged to have
                !       one common update operation on the same buffer
                !       layer
                if(this%current_path%share_update_operations_with(this%head_path)) then
                   
                   call this%current_path%merge(this%head_path)

                   call this%head_path%reinitialize()
                   this%current_path_is_head_path = .true.

                end if

                !2.1.2) the current_path is commit to apply the update
                !       operations on the buffer layer
                call this%current_path%commit(
     $               bf_interface_used,
     $               p_model,
     $               t,
     $               dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1)

                !2.1.3) the gen_coords still need to be stored in either
                !       the head_path or the current_path
                if(this%current_path_is_head_path) then
                   call this%head_path%stage(gen_coords,bf_interface_used)
                else
                   call this%current_path%stage(gen_coords,bf_interface_used)
                end if

             end if

          end if

        end subroutine stage


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> finalize the domain increase by applying the remaining
        !> update operations on the domain extension
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
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

          
          if(this%head_path%share_update_operations_with(this%current_path)) then

             call this%head_path%merge(this%current_path)

             call this%head_path%commit(
     $            bf_interface_used,
     $            p_model,
     $            t,
     $            dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1)

          else

             call this%head_path%commit(
     $            bf_interface_used,
     $            p_model,
     $            t,
     $            dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1)

             call this%current_path%commit(
     $            bf_interface_used,
     $            p_model,
     $            t,
     $            dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1)

          end if

          call this%head_path%remove()
          call this%current_path%remove()

        end subroutine finalize_domain_increase        

      end module icr_interface_class
