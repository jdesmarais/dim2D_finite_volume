      !> @file
      !> module implementing the object encapsulating data
      !> needed for the exchange of information between the
      !> interior domain and the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the abstract version of the
      !> interface between the interior domain and the buffer
      !> layers: only the attributes are implemented
      !
      !> @date
      ! 09_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module interface_abstract_class

        use bf_layer_class        , only : bf_layer
        !use bf_layer_pointer_class, only : bf_layer_pointer
        use parameters_constant   , only : N,S,E,W,N_E,N_W,S_E,S_W

        implicit none


        private
        public :: interface_abstract


        logical, parameter :: debug = .true.


        type :: bf_mainlayer

          type(bf_layer), dimension(:), allocatable :: mainlayer

        end type bf_mainlayer


        type :: bf_mainlayer_pointer

          type(bf_mainlayer), pointer :: pt

        end type bf_mainlayer_pointer



        !> @class interface_abstract
        !> class encapsulating the bf_layer/interior interface
        !> object but only its main attributes are implemented
        !>
        !> @param N_layers
        !> allocatable table which contains the different
        !> sublayers for the north buffer layer
        !>
        !> @param S_layers
        !> allocatable table which contains the different
        !> sublayers for the south buffer layer
        !>
        !> @param E_layers
        !> allocatable table which contains the different
        !> sublayers for the east buffer layer
        !>
        !> @param W_layers
        !> allocatable table which contains the different
        !> sublayers for the west buffer layer
        !>
        !> @param NW_layer
        !> buffer layer representing the NW corner
        !>
        !> @param NE_layer
        !> buffer layer representing the NE corner
        !>
        !> @param SW_layer
        !> buffer layer representing the SW corner
        !>
        !> @param SE_layer
        !> buffer layer representing the SE corner
        !---------------------------------------------------------------
        type :: interface_abstract

          type(bf_mainlayer_pointer), dimension(8) :: mainlayer_pointers

          type(bf_mainlayer), pointer :: N_bf_layers
          type(bf_mainlayer), pointer :: S_bf_layers
          type(bf_mainlayer), pointer :: E_bf_layers
          type(bf_mainlayer), pointer :: W_bf_layers

          type(bf_mainlayer), pointer :: NW_bf_layers
          type(bf_mainlayer), pointer :: NE_bf_layers
          type(bf_mainlayer), pointer :: SW_bf_layers
          type(bf_mainlayer), pointer :: SE_bf_layers

          contains

          procedure, pass :: ini
          procedure, pass :: update_mainlayers_pointers
          procedure, pass :: allocate_bf_mainlayer
          procedure, pass :: get_bf_layer

        end type interface_abstract

        contains


        subroutine ini(this)

          implicit none

          class(interface_abstract), intent(inout) :: this

          this%mainlayer_pointers(N)%pt => this%N_bf_layers
          this%mainlayer_pointers(S)%pt => this%S_bf_layers
          this%mainlayer_pointers(E)%pt => this%E_bf_layers
          this%mainlayer_pointers(W)%pt => this%W_bf_layers

          this%mainlayer_pointers(N_E)%pt => this%NE_bf_layers
          this%mainlayer_pointers(N_W)%pt => this%NW_bf_layers
          this%mainlayer_pointers(S_E)%pt => this%SE_bf_layers
          this%mainlayer_pointers(S_W)%pt => this%SW_bf_layers

        end subroutine ini


        subroutine update_mainlayers_pointers(this)

          implicit none

          class(interface_abstract), intent(inout) :: this

          call this%ini()

        end subroutine update_mainlayers_pointers


        subroutine allocate_space_mainlayer(this, mainlayer_id)

          class(interface_abstract), intent(inout) :: this
          integer                  , intent(in)    :: mainlayer_id

          select case(mainlayer_id)
            case(N)
               allocate(this%N_bf_layers)

            case(S)
               allocate(this%S_bf_layers)

            case(E)
               allocate(this%E_bf_layers)

            case(W)
               allocate(this%W_bf_layers)

            case(N_E)
               allocate(this%NE_bf_layers)

            case(N_W)
               allocate(this%NW_bf_layers)

            case(S_E)
               allocate(this%SE_bf_layers)

            case(S_W)
               allocate(this%SW_bf_layers)

            case default
               print '(''interface_abstract_class'')'
               print '(''allocate_space_mainlayer'')'
               print '(''mainlayer_id not recognized'')'
               print '(''mainlayer_id: '',I2)', mainlayer_id
               stop 'is there a wrong initialization in mainlayer_id ?'
               
          end select

        end subroutine allocate_space_mainlayer


        subroutine allocate_bf_mainlayer(this,mainlayer_id,nb_bf_layers)

          implicit none

          class(interface_abstract), intent(inout) :: this
          integer                  , intent(in)    :: mainlayer_id
          integer                  , intent(in)    :: nb_bf_layers

          integer :: i

          !check if the space for the mainlayer is already associated
          call allocate_space_mainlayer(this,mainlayer_id)
            
          !associate the pointer from the table of pointers to the
          !space just allocated
          call this%update_mainlayers_pointers()          

          !use the pointer just associated to allocate the table
          !containing the buffer layers
          allocate(this%mainlayer_pointers(mainlayer_id)%pt%mainlayer(nb_bf_layers))

          !initialize the buffer layers inside the table using the
          !localization of the main layer
          do i=1, nb_bf_layers
             call this%mainlayer_pointers(mainlayer_id)%pt%mainlayer(nb_bf_layers)%ini([mainlayer_id,1])
          end do

        end subroutine allocate_bf_mainlayer


        function get_bf_layer(this,mainlayer_id, sublayer_id) result(pointer_to_bf_layer)

          class(interface_abstract), intent(in) :: this
          integer                  , intent(in) :: mainlayer_id
          integer, optional        , intent(in) :: sublayer_id
          type(bf_layer), pointer               :: pointer_to_bf_layer
          
          integer :: sublayer_id_I

          if(.not.present(sublayer_id)) then
             sublayer_id_I = 1
          else
             sublayer_id_I = sublayer_id
          end if

          if(debug) then

             if(.not.associated(this%mainlayer_pointers(mainlayer_id)%pt)) then
                print '(''interface_abstract_class'')'
                print '(''get_bf_layer'')'
                print '(''pointer to mainlayer not associated'')'
                stop 'was abstract interface object initialized ?'
             end if

             if(.not.allocated(this%mainlayer_pointers(mainlayer_id)%pt%mainlayer)) then
                print '(''interface_abstract_class'')'
                print '(''get_bf_layer'')'
                print '(''mainlayer not allocated'')'
                print '(''mainlayer asked'', I2)', mainlayer_id                
                stop 'was mainlayer allocated before ?'
             end if

             if(sublayer_id>size(this%mainlayer_pointers(mainlayer_id)%pt%mainlayer)) then
                print '(''interface_abstract_class'')'
                print '(''get_bf_layer'')'
                print '(''sublayer id does not exist'')'
                print '(''mainlayer asked'', I2)', mainlayer_id
                print '(''sublayer asked'', I2)', sublayer_id
                stop 'were sublayers merged and pointer not updated ?'
                stop 'are you really looking for this sublayer ?'
             end if
          end if

          pointer_to_bf_layer => this%mainlayer_pointers(mainlayer_id)%pt%mainlayer(sublayer_id)

        end function get_bf_layer
        
      end module interface_abstract_class
