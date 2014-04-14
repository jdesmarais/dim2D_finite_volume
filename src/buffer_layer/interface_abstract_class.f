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
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W,interior
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind

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
        !
        !> @param N_layers
        !> allocatable table which contains the different
        !> sublayers for the north buffer layer
        !
        !> @param S_layers
        !> allocatable table which contains the different
        !> sublayers for the south buffer layer
        !
        !> @param E_layers
        !> allocatable table which contains the different
        !> sublayers for the east buffer layer
        !
        !> @param W_layers
        !> allocatable table which contains the different
        !> sublayers for the west buffer layer
        !
        !> @param NW_layer
        !> buffer layer representing the NW corner
        !
        !> @param NE_layer
        !> buffer layer representing the NE corner
        !
        !> @param SW_layer
        !> buffer layer representing the SW corner
        !
        !> @param SE_layer
        !> buffer layer representing the SE corner
        !
        !> @param ini
        !> initialize the interface by nullifying all the 
        !> pointers
        !
        !> @param update_mainlayers_pointers
        !> update the pointers to the main layers of the
        !> interface
        !
        !> @param add_sublayer
        !> add a new sublayer to the main layer identified
        !> by its cardinal coordinate
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
          procedure, pass :: get_sublayer

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

          integer :: i

          !nullify all the pointers of the mainlayer_pointers table
          do i=1, size(this%mainlayer_pointers,1)
             nullify(this%mainlayer_pointers(i)%ptr)
          end do

          !nullify all the pointers to mainlayers
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
        function get_sublayer(this, general_coord, local_coord) result(sublayer)

          implicit none

          class(interface_abstract)   , intent(in)  :: this
          integer(ikind), dimension(2), intent(in)  :: general_coord
          integer(ikind), dimension(2), intent(out) :: local_coord
          type(bf_sublayer), pointer                :: sublayer


          integer                     :: mainlayer_id
          type(bf_mainlayer), pointer :: mainlayer
          logical                     :: grdpt_in_sublayer

          !identification of the main layer
          mainlayer_id = get_mainlayer_id(general_coord)

          !identification of the sublayer
          if(debug) then
             if(.not.associated(this%mainlayer_pointers(mainlayer_id)%ptr)) then
                print '(''interface_abstract_class'')'
                print '(''get_sublayer'')'
                print '(''mainlayer not associated'')'
                stop 'the coords do not match any existing buffer layer' 
             end if
          end if
          mainlayer => this%mainlayer_pointers(mainlayer_id)%ptr

          if(debug) then
             if(.not.associated(mainlayer%head_sublayer)) then
                print '(''interface_abstract_class'')'
                print '(''get_sublayer'')'
                print '(''mainlayer%head not associated'')'
                stop 'the coords do not match any existing buffer layer'
             end if
          end if
          sublayer => mainlayer%head_sublayer

          select case(mainlayer_id)
            case(N,S)

               !check if the grid point belongs to the current sublayer
               grdpt_in_sublayer =
     $              (general_coord(1).ge.(sublayer%element%alignment(1,1)-bc_size))
     $              .and.(general_coord(1).le.(sublayer%element%alignment(1,2)+bc_size))

               !go through the different sublayers
               do while(.not.grdpt_in_sublayer)
                  
                  if(.not.associated(sublayer%next)) then
                     print '(''interface_abstract_class'')'
                     print '(''get_sublayer'')'
                     print '(''mainlayer%head not associated'')'
                     stop 'no match for existing buffer layer' 
                  end if

                  sublayer => sublayer%next
                  grdpt_in_sublayer = (general_coord(1).ge.(sublayer%element%alignment(1,1)-bc_size))
     $              .and.(general_coord(1).le.(sublayer%element%alignment(1,2)+bc_size))

               end do

            case(E,W)
               
               !check if the grid point belongs to the current sublayer
               grdpt_in_sublayer =
     $              (general_coord(2).ge.(sublayer%element%alignment(2,1)-bc_size))
     $              .and.(general_coord(2).le.(sublayer%element%alignment(2,2)+bc_size))

               !go through the different sublayers
               do while(.not.grdpt_in_sublayer)
                  
                  if(.not.associated(sublayer%next)) then
                     print '(''interface_abstract_class'')'
                     print '(''get_sublayer'')'
                     print '(''mainlayer%head not associated'')'
                     stop 'no match for existing buffer layer' 
                  end if

                  sublayer => sublayer%next
                  grdpt_in_sublayer = (general_coord(2).ge.(sublayer%element%alignment(2,1)-bc_size))
     $              .and.(general_coord(2).le.(sublayer%element%alignment(2,2)+bc_size))

               end do
          end select

          !compute the local coordinates
          local_coord = sublayer%element%get_local_coord(general_coord)

        end function get_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine converting general coordinates into
        !> the main layer ID (N,S,E,W,N_E,N_W,S_E,S_W)
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param general_coord
        !> integer table giving the general coordinates
        !
        !>@param mainlayer_id
        !> main layer cardinal coordinates
        !--------------------------------------------------------------
        function get_mainlayer_id(general_coord) result(mainlayer_id)

          implicit none

          integer(ikind), dimension(2), intent(in) :: general_coord
          integer                                  :: mainlayer_id

          if(general_coord(2).lt.1) then
             if(general_coord(1).lt.1) then
                mainlayer_id = S_W
             else
                if(general_coord(1).le.ny) then
                   mainlayer_id = S
                else
                   mainlayer_id = S_E
                end if
             end if

          else
             if(general_coord(2).le.ny) then
                if(general_coord(1).lt.1) then
                   mainlayer_id = W
                else
                   if(general_coord(1).le.nx) then
                      mainlayer_id = interior
                      print '(''interface_abstract_class'')'
                      print '(''get_mainlayer_id'')'
                      print '(''main_layer_id = interior'')'
                      stop 'the interior should not be access this way'
                   else
                      mainlayer_id = E
                   end if
                end if
                      
             else
                if(general_coord(1).lt.1) then
                   mainlayer_id = N_W
                else
                   if(general_coord(1).le.nx) then
                      mainlayer_id = N
                   else
                      mainlayer_id = N_E
                   end if
                end if

             end if
          end if

        end function get_mainlayer_id

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
