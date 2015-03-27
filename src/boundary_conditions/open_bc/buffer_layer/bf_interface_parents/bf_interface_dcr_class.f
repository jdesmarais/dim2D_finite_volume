      !> @file
      !> bf_interface_icr object augmented with procedures to remove
      !> the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_icr object augmented with procedures to remove
      !> the buffer layers
      !
      !> @date
      ! 26_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_dcr_class

        use bf_interface_icr_class, only :
     $       bf_interface_icr

        use bf_sublayer_class, only :
     $       bf_sublayer

        use dcr_interface_class, only :
     $       dcr_interface

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_interface_dcr
        

        !>@class bf_interface
        !> bf_interface_icr object augmented with procedures to remove
        !> the buffer layers
        !
        !>@param remove_inactivated_bf_layers
        !> remove the buffer layers that do not have activated grid-points
        !> at the interface with the interior domain
        !
        !>@param adapt_domain_extension
        !> adapt the configuration and extents of the domain extension
        !---------------------------------------------------------------
        type, extends(bf_interface_icr) :: bf_interface_dcr

           contains

           procedure, pass :: remove_inactivated_bf_layers
           procedure, pass :: adapt_domain_extension

        end type bf_interface_dcr

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the extent of the domain extension by updating
        !> the bc_interior_pt at the edge of the computational
        !> domain and removing unactivated buffer layers
        !
        !> @date
        !> 26_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !>@param interior_nodes1
        !> nodes of the interior domain at t
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> timestep
        !
        !>@param interior_bc_sections
        !> boundary sections for the interior domain
        !--------------------------------------------------------------
        subroutine remove_inactivated_bf_layers(
     $       this,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model)

          implicit none

          class(bf_interface_dcr)         , intent(inout) :: this
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model

          
          type(dcr_interface) :: dcr_interface_used

          type(bf_sublayer), pointer :: bf_sublayer_ptr
          type(bf_sublayer), pointer :: bf_sublayer_next

          integer :: k
          integer :: m
          integer :: nb_sublayers
          logical :: can_be_checked
          logical :: can_remain


          !> initialize the dcr_interface by allocating
          !> space for references to the buffer layers
          !> analyzed for removal
          call dcr_interface_used%ini(this)


          !> the buffer layers of the domain extension
          !> are analyzed to see whether they can be removed
          !> or whether their removal is conditioned by its
          !> neighbors
          do k=1,4

             nb_sublayers = this%mainlayer_pointers(k)%get_nb_sublayers()

             bf_sublayer_ptr  => this%mainlayer_pointers(k)%get_head_sublayer()

             do m=1, nb_sublayers

                !> the reference to the next sublayer in the main layer
                !> is initialized at the beginning as if the buffer layer
                !> is removed, the reference is lost
                bf_sublayer_next => bf_sublayer_ptr%get_next()

                
                !> check whether the removal of the buffer layer
                !> should be checked
                can_be_checked = dcr_interface_used%not_in_no_check_list(k,bf_sublayer_ptr)


                !> if the removal of the buffer layer can be checked
                if(can_be_checked) then

                   !> check whether the buffer layer can be removed
                   can_remain = bf_sublayer_ptr%should_remain(
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes,
     $                  p_model)

                   !> if it can not be removed, we prevent its neighbors
                   !> from being removed either
                   if(can_remain) then
                      
                      call dcr_interface_used%prevent_neighbor_removal(
     $                     this,
     $                     k,
     $                     bf_sublayer_ptr)

                   !> if it can be removed, we stage the buffer layer
                   !> for removal: if the buffer layer does not depends
                   !> on its neighbors, it is removed immediately;
                   !> otheriwse, its removal can be unlocked in a
                   !> second check later (finalize_domain_decrease)
                   else

                      call dcr_interface_used%stage(
     $                     this,
     $                     k,
     $                     bf_sublayer_ptr)

                   end if
                      
                end if

                !> get the next buffer layer in the main layer
                bf_sublayer_ptr => bf_sublayer_next

             end do

          end do             


          !> the buffer layers of the domain extension
          !> that could have been removed if they did not
          !> depend on thier neighbors are checked for the
          !> second time and the removal can be unlocked
          !> the dcr_interface object is also finalized:
          !> the array allocated to store the references
          !> to the buffer layer to be removed are
          !> deallocated
          call dcr_interface_used%finalize_domain_decrease(
     $         this)


        end subroutine remove_inactivated_bf_layers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the extent of the domain extension by updating
        !> the bc_interior_pt at the edge of the computational
        !> domain and removing unactivated buffer layers
        !
        !> @date
        !> 26_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !>@param interior_nodes1
        !> nodes of the interior domain at t
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> timestep
        !
        !>@param interior_bc_sections
        !> boundary sections for the interior domain
        !--------------------------------------------------------------
        subroutine adapt_domain_extension(
     $       this,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       p_model,
     $       t,
     $       dt,
     $       interior_bc_sections)

          implicit none

          class(bf_interface_dcr)                    , intent(inout) :: this
          real(rkind), dimension(nx)                 , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)                 , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes1
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: interior_bc_sections


          !> remove the buffer layers whose grid-points at the edge
          !> with the interior domain are desactivated
          call remove_inactivated_bf_layers(
     $         this,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1,
     $         p_model)


          !> increase the buffer layers whose bc_interior_pt are
          !> activated
          call this%bf_interface_icr%adapt_domain_extension(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         p_model,
     $         t,
     $         dt,
     $         interior_bc_sections)

        end subroutine adapt_domain_extension

      end module bf_interface_dcr_class
