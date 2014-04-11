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

        use bf_sublayer_class  , only : bf_sublayer
        use bf_mainlayer_class , only : bf_mainlayer
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : rkind

        implicit none


        private
        public :: interface_abstract


        logical, parameter :: debug = .true.    



        !> @class bf_mainlayer_pointer
        !> class encapsulating a pointer to a buffer main layer
        !>
        !> @param ptr
        !> pointer to a buffer main layer
        !---------------------------------------------------------------
        type :: bf_mainlayer_pointer

          type(bf_mainlayer), pointer :: ptr

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

          type(bf_mainlayer), pointer :: NE_bf_layers
          type(bf_mainlayer), pointer :: NW_bf_layers
          type(bf_mainlayer), pointer :: SE_bf_layers
          type(bf_mainlayer), pointer :: SW_bf_layers

          contains

          procedure, pass :: ini
          procedure, pass :: update_mainlayers_pointers
          procedure, pass :: add_sublayer

          !procedure, pass :: allocate_bf_mainlayer
          !procedure, pass :: get_bf_layer

        end type interface_abstract

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the interface by nullifying
        !> all the pointers to the buffer main layers
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_abstract class encapsulating the pointers
        !> to the buffer main layers
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(interface_abstract), intent(inout) :: this

          nullify(this%N_bf_layers)
          nullify(this%S_bf_layers)
          nullify(this%E_bf_layers)
          nullify(this%W_bf_layers)
          nullify(this%NE_bf_layers)
          nullify(this%NW_bf_layers)
          nullify(this%SE_bf_layers)
          nullify(this%SW_bf_layers)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine updating the interface pointers
        !> to the main layers
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_abstract class encapsulating the pointers
        !> to the buffer main layers
        !--------------------------------------------------------------
        subroutine update_mainlayers_pointers(this)

          implicit none

          class(interface_abstract), intent(inout) :: this

          if(associated(this%N_bf_layers)) this%mainlayer_pointers(N)%ptr => this%N_bf_layers
          if(associated(this%S_bf_layers)) this%mainlayer_pointers(S)%ptr => this%S_bf_layers
          if(associated(this%E_bf_layers)) this%mainlayer_pointers(E)%ptr => this%E_bf_layers
          if(associated(this%W_bf_layers)) this%mainlayer_pointers(W)%ptr => this%W_bf_layers

          if(associated(this%NE_bf_layers)) this%mainlayer_pointers(N_E)%ptr => this%NE_bf_layers
          if(associated(this%NW_bf_layers)) this%mainlayer_pointers(N_W)%ptr => this%NW_bf_layers
          if(associated(this%SE_bf_layers)) this%mainlayer_pointers(S_E)%ptr => this%SE_bf_layers
          if(associated(this%SW_bf_layers)) this%mainlayer_pointers(S_W)%ptr => this%SW_bf_layers

        end subroutine update_mainlayers_pointers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine updating the interface pointers
        !> to the main layers
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> interface_abstract class encapsulating the pointers
        !> to the buffer main layers
        !
        !>@param mainlayer_id
        !> integer identifying the main layer
        !
        !>@param alignment
        !> table if integers identifying the position of the buffer layer
        !> compared to the interior nodes table
        !
        !>@param added_sublayer
        !> pointer to the newly added sublayer
        !--------------------------------------------------------------
        function add_sublayer(
     $     this,
     $     mainlayer_id,
     $     alignment,
     $     nodes,
     $     neighbors)
     $     result(added_sublayer)
        
          class(interface_abstract)       , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          integer, dimension(2,2)         , intent(in)    :: alignment
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          logical, dimension(4)           , intent(in)    :: neighbors
          type(bf_sublayer), pointer                      :: added_sublayer

          type(bf_mainlayer), pointer :: mainlayer_asked

          !first check if the mainlayer corresponding to the cardinal
          !point is indeed allocated
          !if the memory space is not allocated, the space is first
          !allocated, the pointer identifying the mainlayer is initialized
          !and the main layer itself is initialized
          select case(mainlayer_id)
            case(N)
               if(.not.associated(this%N_bf_layers)) then
                  allocate(this%N_bf_layers)
                  this%mainlayer_pointers(N)%ptr => this%N_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(S)
               if(.not.associated(this%S_bf_layers)) then
                  allocate(this%S_bf_layers)
                  this%mainlayer_pointers(S)%ptr => this%S_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(E)
               if(.not.associated(this%E_bf_layers)) then
                  allocate(this%E_bf_layers)
                  this%mainlayer_pointers(E)%ptr => this%E_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(W)
               if(.not.associated(this%W_bf_layers)) then
                  allocate(this%W_bf_layers)
                  this%mainlayer_pointers(W)%ptr => this%W_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(N_W)
               if(.not.associated(this%NW_bf_layers)) then
                  allocate(this%NW_bf_layers)
                  this%mainlayer_pointers(N_W)%ptr => this%NW_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(N_E)
               if(.not.associated(this%NE_bf_layers)) then
                  allocate(this%NE_bf_layers)
                  this%mainlayer_pointers(N_E)%ptr => this%NE_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(S_W)
               if(.not.associated(this%SW_bf_layers)) then
                  allocate(this%SW_bf_layers)
                  this%mainlayer_pointers(S_W)%ptr => this%SW_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case(S_E)
               if(.not.associated(this%SE_bf_layers)) then
                  allocate(this%SE_bf_layers)
                  this%mainlayer_pointers(S_E)%ptr => this%SE_bf_layers
                  call this%mainlayer_pointers(mainlayer_id)%ptr%ini(mainlayer_id)
               end if
            case default
               print '(''interface_abstract_class'')'
               print '(''add_sublayer'')'
               print '(''mainlayer_id not recognized'')'
               print '(''mainlayer asked:'',I2)', mainlayer_id

          end select
                  

          !now that we are sure that space is allocated for the main layer,
          !the sublayer can be integrated to the mainlayer and the buffer
          !layer can be initialized with the localization of the main
          !buffer layer
          mainlayer_asked => this%mainlayer_pointers(mainlayer_id)%ptr
          added_sublayer  => mainlayer_asked%add_sublayer(alignment)
          call added_sublayer%element%ini(mainlayer_id)
          call added_sublayer%element%allocate_bf_layer(
     $         alignment, nodes, neighbors)

        end function add_sublayer


c$$$        function get_bf_layer(this,mainlayer_id, sublayer_id) result(pointer_to_bf_layer)
c$$$
c$$$          class(interface_abstract), intent(in) :: this
c$$$          integer                  , intent(in) :: mainlayer_id
c$$$          integer, optional        , intent(in) :: sublayer_id
c$$$          type(bf_layer), pointer               :: pointer_to_bf_layer
c$$$          
c$$$          integer :: sublayer_id_I
c$$$
c$$$          if(.not.present(sublayer_id)) then
c$$$             sublayer_id_I = 1
c$$$          else
c$$$             sublayer_id_I = sublayer_id
c$$$          end if
c$$$
c$$$          if(debug) then
c$$$
c$$$             if(.not.associated(this%mainlayer_pointers(mainlayer_id)%pt)) then
c$$$                print '(''interface_abstract_class'')'
c$$$                print '(''get_bf_layer'')'
c$$$                print '(''pointer to mainlayer not associated'')'
c$$$                stop 'was abstract interface object initialized ?'
c$$$             end if
c$$$
c$$$             if(.not.allocated(this%mainlayer_pointers(mainlayer_id)%pt%mainlayer)) then
c$$$                print '(''interface_abstract_class'')'
c$$$                print '(''get_bf_layer'')'
c$$$                print '(''mainlayer not allocated'')'
c$$$                print '(''mainlayer asked'', I2)', mainlayer_id                
c$$$                stop 'was mainlayer allocated before ?'
c$$$             end if
c$$$
c$$$             if(sublayer_id>size(this%mainlayer_pointers(mainlayer_id)%pt%mainlayer)) then
c$$$                print '(''interface_abstract_class'')'
c$$$                print '(''get_bf_layer'')'
c$$$                print '(''sublayer id does not exist'')'
c$$$                print '(''mainlayer asked'', I2)', mainlayer_id
c$$$                print '(''sublayer asked'', I2)', sublayer_id
c$$$                stop 'were sublayers merged and pointer not updated ?'
c$$$                stop 'are you really looking for this sublayer ?'
c$$$             end if
c$$$          end if
c$$$
c$$$          pointer_to_bf_layer => this%mainlayer_pointers(mainlayer_id)%pt%mainlayer(sublayer_id)
c$$$
c$$$        end function get_bf_layer
        
      end module interface_abstract_class
