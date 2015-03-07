      !> @file
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @date
      ! 07_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_basic_class

        use bf_sublayer_class, only :
     $       bf_sublayer

        use bf_layer_errors_module, only :
     $       error_mainlayer_interface_type,
     $       error_mainlayer_interface_incompatible

        use parameters_bf_layer, only :
     $       NE_interface_type,
     $       NW_interface_type,
     $       SE_interface_type,
     $       SW_interface_type

        use parameters_constant, only :
     $       N,S,E,W

        implicit none

        private
        public :: mainlayer_interface_basic


        !>@class mainlayer_interface_basic
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !>@param NE_interface_N_ptr
        !> pointer to the North main layer at the NE interface
        !
        !>@param NE_interface_E_ptr
        !> pointer to the East main layer at the NE interface
        !
        !>@param NW_interface_N_ptr
        !> pointer to the North main layer at the NW interface
        !
        !>@param NW_interface_W_ptr
        !> pointer to the West main layer at the NW interface
        !
        !>@param SE_interface_S_ptr
        !> pointer to the South main layer at the SE interface
        !
        !>@param SE_interface_E_ptr
        !> pointer to the East main layer at the SE interface
        !
        !>@param SW_interface_S_ptr
        !> pointer to the South main layer at the SW interface
        !
        !>@param SW_interface_W_ptr
        !> pointer to the West main layer at the SW interface
        !--------------------------------------------------------------
        type :: mainlayer_interface_basic

          !NE interface
          type(bf_sublayer), pointer :: NE_interface_N_ptr
          type(bf_sublayer), pointer :: NE_interface_E_ptr

          !NW interface
          type(bf_sublayer), pointer :: NW_interface_N_ptr
          type(bf_sublayer), pointer :: NW_interface_W_ptr

          !SE interface
          type(bf_sublayer), pointer :: SE_interface_S_ptr
          type(bf_sublayer), pointer :: SE_interface_E_ptr

          !SW interface
          type(bf_sublayer), pointer :: SW_interface_S_ptr
          type(bf_sublayer), pointer :: SW_interface_W_ptr

          contains

          procedure, pass :: ini
          procedure, pass :: set_mainlayer_interface_bf_layer
          procedure, pass :: remove_mainlayer_interface_bf_layer

          procedure, pass :: enable_mainlayer_interface
          procedure, pass :: disable_mainlayer_interface

        end type mainlayer_interface_basic


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the links to the buffer layers at the interface
        !> between main layers
        !
        !> @date
        !> 07_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(mainlayer_interface_basic), intent(inout) :: this

          !initialize the links to no buffer layer
          nullify(this%NE_interface_N_ptr)
          nullify(this%NE_interface_E_ptr)

          nullify(this%NW_interface_N_ptr)
          nullify(this%NW_interface_W_ptr)

          nullify(this%SE_interface_S_ptr)
          nullify(this%SE_interface_E_ptr)

          nullify(this%SW_interface_S_ptr)
          nullify(this%SW_interface_W_ptr)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set as buffer layer at the interface the buffer layer
        !> passed as argument
        !
        !> @date
        !> 07_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param mainlayer_interface_type
        !> interger identifying the main layer interface
        !
        !> @param bf_sublayer_ptr
        !> pointer to the sublayer set at the interface
        !------------------------------------------------------------
        subroutine set_mainlayer_interface_bf_layer(
     $     this,
     $     mainlayer_interface_type,
     $     bf_sublayer_ptr)

          implicit none

          class(mainlayer_interface_basic), intent(inout) :: this
          integer                         , intent(in) :: mainlayer_interface_type
          type(bf_sublayer), pointer      , intent(in) :: bf_sublayer_ptr

          integer :: localization


          localization = bf_sublayer_ptr%get_localization()


          select case(mainlayer_interface_type)

            case(NE_interface_type)
               select case(localization)
                 case(N)
                    this%NE_interface_N_ptr => bf_sublayer_ptr
                    if(associated(this%NE_interface_E_ptr)) then
                       call enable_mainlayer_interface(this,NE_interface_type)
                    else
                       call disable_mainlayer_interface(this,NE_interface_type)
                    end if

                 case(E)
                    this%NE_interface_E_ptr => bf_sublayer_ptr
                    if(associated(this%NE_interface_N_ptr)) then
                       call enable_mainlayer_interface(this,NE_interface_type)
                    else
                       call disable_mainlayer_interface(this,NE_interface_type)
                    end if

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)

               end select

            case(NW_interface_type)
               select case(localization)
                 case(N)
                    this%NW_interface_N_ptr => bf_sublayer_ptr
                    if(associated(this%NW_interface_W_ptr)) then
                       call enable_mainlayer_interface(this,NW_interface_type)
                    else
                       call disable_mainlayer_interface(this,NW_interface_type)
                    end if

                 case(W)
                    this%NW_interface_W_ptr => bf_sublayer_ptr
                    if(associated(this%NW_interface_N_ptr)) then
                       call enable_mainlayer_interface(this,NW_interface_type)
                    else
                       call disable_mainlayer_interface(this,NW_interface_type)
                    end if

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)

               end select

            case(SE_interface_type)
               select case(localization)
                 case(S)
                    this%SE_interface_S_ptr => bf_sublayer_ptr
                    if(associated(this%SE_interface_E_ptr)) then
                       call enable_mainlayer_interface(this,SE_interface_type)
                    else
                       call disable_mainlayer_interface(this,SE_interface_type)
                    end if

                 case(E)
                    this%SE_interface_E_ptr => bf_sublayer_ptr
                    if(associated(this%SE_interface_S_ptr)) then
                       call enable_mainlayer_interface(this,SE_interface_type)
                    else
                       call disable_mainlayer_interface(this,SE_interface_type)
                    end if

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)
               end select

            case(SW_interface_type)
               select case(localization)
                 case(S)
                    this%SW_interface_S_ptr => bf_sublayer_ptr
                    if(associated(this%SW_interface_W_ptr)) then
                       call enable_mainlayer_interface(this,SW_interface_type)
                    else
                       call disable_mainlayer_interface(this,SW_interface_type)
                    end if

                 case(W)
                    this%SW_interface_W_ptr => bf_sublayer_ptr
                    if(associated(this%SW_interface_S_ptr)) then
                       call enable_mainlayer_interface(this,SW_interface_type)
                    else
                       call disable_mainlayer_interface(this,SW_interface_type)
                    end if

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)
               end select

            case default
               call error_mainlayer_interface_type(
     $              'mainlayer_interface_basic_class',
     $              'set_interface_bf_layer',
     $              mainlayer_interface_type)

          end select

        end subroutine set_mainlayer_interface_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the buffer layer at the interface
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param mainlayer_interface_type
        !> interger identifying the main layer interface
        !
        !> @param bf_sublayer_ptr
        !> pointer to the sublayer set at the interface
        !------------------------------------------------------------
        subroutine remove_mainlayer_interface_bf_layer(
     $     this,
     $     mainlayer_interface_type,
     $     bf_sublayer_ptr)

          implicit none

          class(mainlayer_interface_basic), intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_interface_type
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr

          integer :: localization


          !disable the exchange of nodes at the interface
          !between main layer interfaces
          call disable_mainlayer_interface(this,mainlayer_interface_type)


          !nullify the pointer corresponding to the buffer layer
          localization = bf_sublayer_ptr%get_localization()

          select case(mainlayer_interface_type)

            case(NE_interface_type)

               select case(localization)
                 case(N)
                    nullify(this%NE_interface_N_ptr)

                 case(E)
                    nullify(this%NE_interface_E_ptr)

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)

               end select

            case(NW_interface_type)

               select case(localization)
                 case(N)
                    nullify(this%NW_interface_N_ptr)

                 case(W)
                    nullify(this%NW_interface_W_ptr)

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)

               end select

            case(SE_interface_type)

               select case(localization)
                 case(S)
                    nullify(this%SE_interface_S_ptr)

                 case(E)
                    nullify(this%SE_interface_E_ptr)

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)
               end select

            case(SW_interface_type)

               select case(localization)
                 case(S)
                    nullify(this%SW_interface_S_ptr)

                 case(W)
                    nullify(this%SW_interface_W_ptr)

                 case default
                    call error_mainlayer_interface_incompatible(
     $                   'mainlayer_interface_basic_class',
     $                   'set_interface_bf_layer',
     $                   mainlayer_interface_type,
     $                   localization)
               end select

            case default
               call error_mainlayer_interface_type(
     $              'mainlayer_interface_basic_class',
     $              'set_interface_bf_layer',
     $              mainlayer_interface_type)

          end select

        end subroutine remove_mainlayer_interface_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> enable the exchange of nodes at a mainlayer interface
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param mainlayer_interface_type
        !> interger identifying the main layer interface
        !------------------------------------------------------------
        subroutine enable_mainlayer_interface(this,mainlayer_interface_type)

          implicit none

          class(mainlayer_interface_basic), intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_interface_type


          select case(mainlayer_interface_type)

            case(NE_interface_type)
               call this%NE_interface_N_ptr%set_neighbor2_share(.true.)
               call this%NE_interface_E_ptr%set_neighbor2_share(.true.)

            case(NW_interface_type)
               call this%NW_interface_N_ptr%set_neighbor1_share(.true.)
               call this%NW_interface_W_ptr%set_neighbor2_share(.true.)

            case(SE_interface_type)
               call this%SE_interface_S_ptr%set_neighbor2_share(.true.)
               call this%SE_interface_E_ptr%set_neighbor1_share(.true.)

            case(SW_interface_type)
               call this%SW_interface_S_ptr%set_neighbor1_share(.true.)
               call this%SW_interface_W_ptr%set_neighbor1_share(.true.)

            case default
               call error_mainlayer_interface_type(
     $              'mainlayer_interface_basic_class',
     $              'set_interface_bf_layer',
     $              mainlayer_interface_type)
          end select

        end subroutine enable_mainlayer_interface


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> disable the exchange of nodes at a mainlayer interface
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param mainlayer_interface_type
        !> interger identifying the main layer interface
        !------------------------------------------------------------
        subroutine disable_mainlayer_interface(this,mainlayer_interface_type)

          implicit none

          class(mainlayer_interface_basic), intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_interface_type


          select case(mainlayer_interface_type)

            case(NE_interface_type)
               if(associated(this%NE_interface_N_ptr)) then
                  call this%NE_interface_N_ptr%set_neighbor2_share(.false.)
               end if
               if(associated(this%NE_interface_E_ptr)) then
                  call this%NE_interface_E_ptr%set_neighbor2_share(.false.)
               end if

            case(NW_interface_type)
               if(associated(this%NW_interface_N_ptr)) then
                  call this%NW_interface_N_ptr%set_neighbor1_share(.false.)
               end if
               if(associated(this%NW_interface_W_ptr)) then
                  call this%NW_interface_W_ptr%set_neighbor2_share(.false.)
               end if

            case(SE_interface_type)
               if(associated(this%SE_interface_S_ptr)) then
                  call this%SE_interface_S_ptr%set_neighbor2_share(.false.)
               end if
               if(associated(this%SE_interface_E_ptr)) then
                  call this%SE_interface_E_ptr%set_neighbor1_share(.false.)
               end if

            case(SW_interface_type)
               if(associated(this%SW_interface_S_ptr)) then
                  call this%SW_interface_S_ptr%set_neighbor1_share(.false.)
               end if
               if(associated(this%SW_interface_W_ptr)) then
                  call this%SW_interface_W_ptr%set_neighbor1_share(.false.)
               end if

            case default
               call error_mainlayer_interface_type(
     $              'mainlayer_interface_basic_class',
     $              'set_interface_bf_layer',
     $              mainlayer_interface_type)
          end select

        end subroutine disable_mainlayer_interface
        

      end module mainlayer_interface_basic_class
