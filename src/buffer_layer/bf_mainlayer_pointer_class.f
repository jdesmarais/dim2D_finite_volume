      module bf_mainlayer_pointer_class

        use bf_sublayer_class , only : bf_sublayer
        use bf_mainlayer_class, only : bf_mainlayer
        use parameters_input  , only : nx,ny,ne, debug
        use parameters_kind   , only : ikind, rkind


        implicit none


        private
        public :: bf_mainlayer_pointer


        !> @class bf_mainlayer_pointer
        !> pointer to a buffer main layer
        !>
        !> @param ptr
        !> pointer to a buffer main layer
        !---------------------------------------------------------------
        type :: bf_mainlayer_pointer

          type(bf_mainlayer), pointer, private :: ptr

          contains

          procedure, pass :: ini
          procedure, pass :: get_ptr
          procedure, pass :: set_ptr
          procedure, pass :: nullify_ptr
          procedure, pass :: allocate_ptr
          procedure, pass :: deallocate_ptr
          procedure, pass :: associated_ptr

          procedure, pass :: ini_mainlayer
          procedure, pass :: get_mainlayer_id
          procedure, pass :: get_nb_sublayers
          procedure, pass :: get_head_sublayer
          procedure, pass :: get_tail_sublayer
          procedure, pass :: add_sublayer
          procedure, pass :: merge_sublayers
          procedure, pass :: print_binary
          
        end type bf_mainlayer_pointer


        contains


        !< initialize the ptr attribute by nullifying it
        subroutine ini(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          call this%nullify_ptr()

        end subroutine ini        


        !< get the buffer main layer to which the
        !> ptr attribute is pointing at 
        function get_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          type(bf_mainlayer), pointer             :: get_ptr

          if(associated(this%ptr)) then
             get_ptr => this%ptr
          else
             nullify(get_ptr)
          end if

        end function get_ptr


        !< set the buffer main layer to which the
        !> ptr attribute is pointing at 
        subroutine set_ptr(this, bf_mainlayer_ptr)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this
          type(bf_mainlayer), pointer, intent(in)    :: bf_mainlayer_ptr

          this%ptr => bf_mainlayer_ptr

        end subroutine set_ptr


        !< nullify the ptr attribute
        subroutine nullify_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          nullify(this%ptr)

        end subroutine nullify_ptr


        !< allocate space for the ptr attribute
        subroutine allocate_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          allocate(this%ptr)

        end subroutine allocate_ptr


        !< deallocate space for the ptr attribute
        subroutine deallocate_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          deallocate(this%ptr)
          nullify(this%ptr)

        end subroutine deallocate_ptr


        !< check if the ptr attribute is associated
        function associated_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          logical                                 :: associated_ptr

          associated_ptr = associated(this%ptr)

        end function associated_ptr


        !< initialize the mainlayer corresponding to the pointer
        subroutine ini_mainlayer(this, mainlayer_id)

          class(bf_mainlayer_pointer), intent(inout) :: this
          integer                    , intent(in)    :: mainlayer_id

          !debug: check mainlayer_id
          if(debug) then
             if((mainlayer_id.lt.1).or.(mainlayer_id.gt.8)) then
                print '(''bf_interface_class'')'
                print '(''add_sublayer'')'
                print '(''mainlyer_id not recognized'')'
                print '(''mainlayer_id: '',I2)', mainlayer_id
                stop 'change mainlayer_id'
             end if
          end if

          !check ptr attribute association
          if(this%associated_ptr()) then
             print '(''bf_mainlayer_pointer_class'')'
             print '(''ini'')'
             print '(''ptr attribute is already associated'')'
             print '(''- either the ptr attribute is really'')'
             print '(''  associated to space in memory and'')'
             print '(''  the initialization will lead to a)'')'
             print '(''  memory leak'')'
             print '(''- or the object has not been initialized'')'
             print '(''  before begin used which is bad practice'')'
             stop 'check the call order to this object'
          end if

          !allocate space for the mainlayer and initialize it
          call this%allocate_ptr()
          call this%ptr%ini(mainlayer_id)          

        end subroutine ini_mainlayer

      
        !< get the cardinal coordinate of the buffer main layer
        !> encapsulated
        function get_mainlayer_id(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          integer                                 :: get_mainlayer_id

          if(this%associated_ptr()) then
             get_mainlayer_id = this%ptr%get_mainlayer_id()
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''get_mainlayer_id'')'
             stop 'ptr attribute not associated'
          end if             

        end function get_mainlayer_id


        !< get the number of sublayers contained in the buffer
        !> main layer encapsulated
        function get_nb_sublayers(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          integer                                 :: get_nb_sublayers

          if(this%associated_ptr()) then
             get_nb_sublayers = this%ptr%get_nb_sublayers()
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''get_nb_sublayers'')'
             stop 'ptr attribute not associated'
          end if

        end function get_nb_sublayers


        !< get the first sublayer in the chained list of the 
        !> main layer encapsulated
        function get_head_sublayer(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          type(bf_sublayer)          , pointer    :: get_head_sublayer

          if(this%associated_ptr()) then
             if(associated(this%ptr%get_head_sublayer())) then
                get_head_sublayer => this%ptr%get_head_sublayer()
             else
                nullify(get_head_sublayer)
             end if
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''get_head_sublayer'')'
             stop 'ptr attribute not associated'
          end if

        end function get_head_sublayer


        !< get the last sublayer in the chained list of the 
        !> main layer encapsulated
        function get_tail_sublayer(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          type(bf_sublayer)          , pointer    :: get_tail_sublayer

          if(this%associated_ptr()) then
             if(associated(this%ptr%get_tail_sublayer())) then
                get_tail_sublayer => this%ptr%get_tail_sublayer()
             else
                nullify(get_tail_sublayer)
             end if
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''get_tail_sublayer'')'
             stop 'ptr attribute not associated'
          end if

        end function get_tail_sublayer


        !< add sublayer to the chained list of the mainlayer
        !> encapsulated
        function add_sublayer(this, nodes, alignment)
     $     result(added_sublayer_ptr)

          implicit none

          class(bf_mainlayer_pointer)        , intent(inout) :: this
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: alignment

          type(bf_sublayer), pointer                         :: added_sublayer_ptr

          if(this%associated_ptr()) then
             added_sublayer_ptr => this%ptr%add_sublayer(
     $            nodes, alignment)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''add_sublayer'')'
             stop 'ptr attribute not associated'
          end if
          
        end function add_sublayer


        !< merge two sublayers of the mainlayer encapsulated
        function merge_sublayers(
     $     this,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     interior_nodes,
     $     alignment)
     $     result(merged_sublayer)
        
          implicit none        
        
          class(bf_mainlayer_pointer)             , intent(inout) :: this
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer2
          real(rkind)   , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2), optional, intent(in)    :: alignment
          type(bf_sublayer), pointer                              :: merged_sublayer
        

          if(this%associated_ptr()) then
             merged_sublayer => this%ptr%merge_sublayers(
     $            bf_sublayer1, bf_sublayer2,
     $            interior_nodes, alignment)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''merge_sublayers'')'
             stop 'ptr attribute not associated'
          end if

        end function merge_sublayers

      
        !< print the content of the mainlayer encapsulated
        !> on an output binary file
        subroutine print_binary(
     $     this, suffix_nodes, suffix_grdid, suffix_sizes)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          character(*)               , intent(in) :: suffix_nodes
          character(*)               , intent(in) :: suffix_grdid
          character(*)               , intent(in) :: suffix_sizes


          if(this%associated_ptr()) then
             call this%ptr%print_binary(
     $            suffix_nodes, suffix_grdid, suffix_sizes)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''merge_sublayers'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine print_binary

      end module bf_mainlayer_pointer_class
