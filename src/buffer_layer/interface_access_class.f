      module interface_access_class
      
        !> @class interface_access
        !> class encapsulating the last data access by the
        !> interface
        !
        !> @param nodes
        !> pointer to the last nodes access
        !
        !> @param grdptid
        !> pointer to the last grdptid access
        !
        !> @param general_coord
        !> integer table identifying the last general coordinates
        !> access
        !
        !> @param local_coord
        !> integer table identifying the last local coordinates
        !> access
        !
        !> @param interior
        !> logical identifying whether the last access was to
        !> interior or buffer layers
        !---------------------------------------------------------------
        type :: interface_access

          real(rkind)   , dimension(:,:,:), pointer :: nodes
          real(rkind)   , dimension(:,:)  , pointer :: grdptid
          integer(ikind), dimension(2)              :: general_coord
          integer(ikind), dimension(2)              :: local_coord

          contains

          procedure, pass :: compare_access

        end type interface_access


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function comparing the current access with the previous
        !> indices access
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_access encapsulating the information about
        !> to the previous data access
        !--------------------------------------------------------------
        function compare_access(this, general_coord_asked) result(match)

          implicit none

          class(interface_access)     , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: general_coord_asked

          integer(ikind), dimension(2) :: local_coord

          !compute the new local coordinates due to this change
          local_coord(1) = this%local_coord(1) + general_coord_asked(1) - this%general_coord(1)
          local_coord(2) = this%local_coord(2) + general_coord_asked(2) - this%general_coord(2)
          
          !compare to the size of the pointer previously accessed
          match = (local_coord(1).ge.1).and.(local_coord(1).le.size(this%nodes,1))
     $         .and.(local_coord(2).ge.1).and.(local_coord(2).le.size(this%nodes,2))

        end function compare_access

      end module interface_access_class
