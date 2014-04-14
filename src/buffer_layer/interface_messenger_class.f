      !> @file
      !> module implementing the object encapsulating
      !> functions to access data in the buffer layers
      !> or inside the interior domain using the interface
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating
      !> functions to access data in the buffer layers
      !> or inside the interior domain
      !
      !> @date
      ! 11_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module interface_messenger_class

        use interface_abstract_class, only : abstract_interface
        use interface_access_class  , only : interface_access
        use parameters_input        , only : nx,ny,ne

        implicit none

        private
        public :: interface_messenger

        
        !> @class interface_abstract
        !> class encapsulating the bf_layer/interior interface
        !> object but only its main attributes are implemented
        !
        !> @param b_last_access
        !> pointer to the last before least data access
        !
        !> @param last_access
        !> pointer to the last data access
        !---------------------------------------------------------------
        type :: interface_messenger

          type(interface_access), pointer :: b_last_access
          type(interface_access), pointer :: last_access

          contains

          procedure, pass :: ini

        end type interface_messenger

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the interface messenger
        !> by nullifying all the pointers to the previous
        !> memory access
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_messenger encapsulating the pointers
        !> to the previous data access
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none
          
          class(interface_messenger), intent(inout) :: this

          nullify(b_last_access)
          nullify(last_access)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine exchanging the previous access
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_messenger encapsulating the pointers
        !> to the previous data access
        !--------------------------------------------------------------
        subroutine exchange_access(this)

          implicit none
          
          class(interface_messenger), intent(inout) :: this

          type(interface_access), pointer :: temp_ptr

          if(associated(this%b_last_access)) then

             temp_ptr           => this%b_last_access
             this%b_last_access => this%last_access
             this%last_access   => temp_ptr

          end if

        end subroutine exchange_access


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine delivering the nodes gridpoints asked
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_messenger encapsulating the pointers
        !> to the previous data access
        !--------------------------------------------------------------
        function get_nodes(
     $     this,
     $     interface_used,
     $     interior_nodes,
     $     indices)
     $     result(grdpt)

          implicit none
          
          class(interface_messenger)  , intent(inout) :: this
          class(abstract_interface)   , intent(in)    :: interface_used
          real(rkind)                 , intent(in)    :: interior_nodes
          integer(ikind), dimension(3), intent(in)    :: general_coord
          real(rkind)                                 :: grdpt

          if(interior_grdpt(general_coord)) then
             grdpt = nodes(
     $            general_coord(1),
     $            general_coord(2),
     $            general_coord(3))
          else

             !< compare previous access, if it was previously access,
             !> use the pointers pointing to the previous access
             !> otherwise, look for the data in the buffer
             !> layers
             if(associated(this%last_access)) then
                match = this%last_access%compare_access(general_coord)
                
                if(match) then
                   grdpt = this%last_access%nodes(
     $                  general_coord(1),
     $                  general_coord(2),
     $                  general_coord(3))
                else
                   if(associated(this%b_last_access)) then
                      match = this%b_last_access%compare_access(general_coord)
                
                      if(match) then
                         grdpt = this%b_last_access%nodes(
     $                        general_coord(1),
     $                        general_coord(2),
     $                        general_coord(3))
                      
                         call this%exchange_access()

                      end if

                   !< the data were not access in the previous access saved
                   !> it is required to go through the different buffer
                   !> layers to find the data
                   else
                      
                      sublayer => interface_used%get_sublayer(general_coord)
                      
                      local_coord
                      
                      grdpt = sublayer%element%nodes(
     $                        general_coord(1),
     $                        general_coord(2),
     $                        general_coord(3))

                   end if
                end if
                   

          end if

        end function get_nodes


        function interior_grdpt(general_coord)

          implicit none

          integer(ikind), dimension(3), intent(in) :: general_coord
          logical                                  :: interior_grdpt


          interior_grdpt = (general_coord(1).ge.1).and.(general_coord(1).le.nx)
     $         .and.(general_coord(2).ge.1).and.(general_coord(2).le.ny)

        end function interior_grdpt


        

      end module interface_messenger_class
