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

        use bf_sublayer_class     , only : bf_sublayer
        use bf_mainlayer_class    , only : bf_mainlayer
        use parameters_constant   , only : N,S,E,W,
     $                                     N_E,N_W,S_E,S_W,
     $                                     interior
        use parameters_input      , only : nx,ny,ne,bc_size
        use parameters_kind       , only : ikind, rkind

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

          procedure, pass   :: ini
          procedure, pass   :: update_mainlayers_pointers
          procedure, pass   :: get_mainlayer
          procedure, pass   :: add_sublayer
          procedure, pass   :: get_sublayer
          !procedure, pass   :: merge_sublayers

          procedure, nopass :: get_mainlayer_id

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
        !>@param nodes
        !> table of the interior domain grid points
        !
        !>@param neighbors
        !> table identifying the neighbors for the allocation of the
        !> added sublayer
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
        !>@param general_coord
        !> table giving the general coordinates of the point analyzed
        !
        !>@param local_coord
        !> table giving the local coordinates of the point analyzed
        !> in the corresponding sublayer
        !
        !>@param tolerance_i
        !> integer indicating how far the gridpoint can be from the
        !> closest sublayer to be considered inside
        !
        !>@param sublayer
        !> pointer to the sublayer matching the general coordinates
        !> of the grid point
        !--------------------------------------------------------------
        function get_sublayer(
     $     this,
     $     general_coord,
     $     local_coord,
     $     tolerance_i,
     $     mainlayer_id_i)
     $     result(sublayer)

          implicit none

          class(interface_abstract)   , intent(in)  :: this
          integer(ikind), dimension(2), intent(in)  :: general_coord
          integer(ikind), dimension(2), intent(out) :: local_coord
          integer       , optional    , intent(in)  :: tolerance_i
          integer       , optional    , intent(in)  :: mainlayer_id_i
          type(bf_sublayer), pointer                :: sublayer


          integer                     :: mainlayer_id
          type(bf_mainlayer), pointer :: mainlayer
          integer                     :: tolerance
          logical                     :: grdpt_in_sublayer

          !< identification of the main layer
          if(present(mainlayer_id_i)) then
             mainlayer_id = mainlayer_id_i
          else
             mainlayer_id = get_mainlayer_id(general_coord)
          end if

          !< if the general coordinates match the interior,
          !> no sublayer matches the general coordinates
          !> and the sublayer pointer is nullified
          if(mainlayer_id.eq.interior) then
             nullify(sublayer)

          !< otherwise, the mainlayers are analyzed
          else

             !< check that the main layer exists
             !< if it does not exist, no sublayer can be saved inside
             !< and the pointer to the sublayer is nullified
             if(.not.associated(this%mainlayer_pointers(mainlayer_id)%ptr)) then
                nullify(sublayer)
                  
             !< if the main layer exists, the sublayers saved inside are
             !< checked to decide whether the grid point asked belongs to
             !< one of them or not
             else
                mainlayer => this%mainlayer_pointers(mainlayer_id)%ptr
             
                !< check if sublayers are saved inside the mainlayer
                !> if no sublayers are saved inside the mainlayer,
                !> no existing sublayer can match the general coord
                !> and so the pointer to sublayer is nullified
                if(.not.associated(mainlayer%head_sublayer)) then
                   nullify(sublayer)
             
                !< otherwise, the sublayer corresponding to the general
                !> coordinates is searched by going through the different
                !> element of the doubled chained list
                else
                   sublayer => mainlayer%head_sublayer
             
             	!< processing the tolerance for matching a sublayer
             	!> if no tolerence is provided, the default option is 0
             	if(.not.present(tolerance_i)) then
             	   tolerance=0
             	else
                   tolerance=tolerance_i
                end if
             	
                   !< if the mainlayer investigated is N,S,E or W, there can
                   !> be sublayers to be investigated, otherwise, there are no
                   !> several sublayers in the main layers for NE,NW,SE,SW and
                   !> local coordinates can be directly investigated
             	select case(mainlayer_id)
             	  case(N,S)
             	
             	     !check if the grid point belongs to the current sublayer
             	     grdpt_in_sublayer =
     $                    (general_coord(1).ge.(sublayer%element%alignment(1,1)-bc_size-tolerance))
     $                    .and.(general_coord(1).le.(sublayer%element%alignment(1,2)+bc_size+tolerance))
             	
             	     !go through the different sublayers
             	     do while(.not.grdpt_in_sublayer)
             	        
             	        !if no matching sublayer can be found
             	        !nullify the corresponding pointer
             	        if(.not.associated(sublayer%next)) then
             	           nullify(sublayer)
             	           exit
             	        end if
             	
             	        sublayer => sublayer%next
             	        grdpt_in_sublayer = (general_coord(1).ge.(sublayer%element%alignment(1,1)-bc_size-tolerance))
     $                       .and.(general_coord(1).le.(sublayer%element%alignment(1,2)+bc_size+tolerance))
             	
             	     end do
             	
             	  case(E,W)
             	     
             	     !check if the grid point belongs to the current sublayer
             	     grdpt_in_sublayer =
     $                    (general_coord(2).ge.(sublayer%element%alignment(2,1)-bc_size-tolerance))
     $                    .and.(general_coord(2).le.(sublayer%element%alignment(2,2)+bc_size+tolerance))
             	
             	     !go through the different sublayers
             	     do while(.not.grdpt_in_sublayer)
             	        
             	        !if no matching sublayer can be found
             	        !nullify the corresponding pointer
             	        if(.not.associated(sublayer%next)) then
             	           nullify(sublayer)
             	           exit
             	        end if
             	
             	        sublayer => sublayer%next
             	        grdpt_in_sublayer = (general_coord(2).ge.(sublayer%element%alignment(2,1)-bc_size-tolerance))
     $                       .and.(general_coord(2).le.(sublayer%element%alignment(2,2)+bc_size+tolerance))
             	
             	     end do
             	end select
             
                end if
             end if
          end if

          !< if a sublayer matching the general coordinates was found
          !> compute the local coordinates in this sublayer
          if(associated(sublayer)) then
             local_coord = sublayer%element%get_local_coord(general_coord)
          end if           

        end function get_sublayer


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
        !>@param id
        !> cardinal coordinate of the main layer
        !--------------------------------------------------------------
        function get_mainlayer(this, id)

          implicit none

          class(interface_abstract), intent(in) :: this
          integer                  , intent(in) :: id
          type(bf_mainlayer), pointer           :: get_mainlayer

          
          if(associated(this%mainlayer_pointers(id)%ptr)) then
             get_mainlayer => this%mainlayer_pointers(id)%ptr
          else
             nullify(get_mainlayer)
          end if

        end function get_mainlayer



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

          if(general_coord(2).le.bc_size) then
             if(general_coord(1).le.bc_size) then
                mainlayer_id = S_W
             else
                if(general_coord(1).le.(nx-bc_size)) then
                   mainlayer_id = S
                else
                   mainlayer_id = S_E
                end if
             end if

          else
             if(general_coord(2).le.(ny-bc_size)) then
                if(general_coord(1).le.bc_size) then
                   mainlayer_id = W
                else
                   if(general_coord(1).le.(nx-bc_size)) then
                      mainlayer_id = interior
                   else
                      mainlayer_id = E
                   end if
                end if
                      
             else
                if(general_coord(1).le.bc_size) then
                   mainlayer_id = N_W
                else
                   if(general_coord(1).le.(nx-bc_size)) then
                      mainlayer_id = N
                   else
                      mainlayer_id = N_E
                   end if
                end if

             end if
          end if

        end function get_mainlayer_id

      end module interface_abstract_class
