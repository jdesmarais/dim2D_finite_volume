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

        use icr_path_class, only :
     $     icr_path

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

          logical        :: current_path_is_head-path
          type(icr_path) :: head_path
          type(icr_path) :: current_path

          contains

          procedure, pass :: ini
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

          
          if(this%head_path%shares_update_operations_with(this%current_path)) then

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
