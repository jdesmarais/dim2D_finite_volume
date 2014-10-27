      !> @file
      !> module implementing the pointer to a bf_mainlayer object
      !> such that a table of bf_mainlayer pointers can be created
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the pointer to a bf_mainlayer object
      !> such that a table of bf_mainlayer pointers can be created
      !
      !> @date
      ! 11_04_2014 - initial version      - J.L. Desmarais
      ! 26_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_pointer_class

        use bc_operators_class, only :
     $     bc_operators

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_sublayer_class, only :
     $       bf_sublayer

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_input, only :
     $       nx,ny,ne, debug

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators

        implicit none


        private
        public :: bf_mainlayer_pointer


        !>@class bf_mainlayer_pointer
        !> pointer to a bf_mainlayer object
        !
        !>@param ptr
        !> pointer to a bf_mainlayer object
        !
        !>@param ini
        !> initialize the bf_mainlayer_ptr object by nullify its
        !> ptr attribute
        !
        !>@param get_ptr
        !> get the ptr attribute
        !
        !>@param set_ptr
        !> set the ptr attribute
        !
        !>@param nullify_ptr
        !> nullify the ptr attribute
        !
        !>@param allocate_ptr
        !> allocate space for the ptr attribute
        !
        !>@param deallocate_ptr
        !> deallocate space for the ptr attribute
        !
        !>@param associated_ptr
        !> check if the ptr attribute is associated
        !
        !>@param ini_mainlayer
        !> initialize the mainlayer corresponding to the
        !> ptr attribute
        !
        !>@param get_mainlayer_id
        !> get the mainlayer_id attribute
        !
        !>@param get_nb_sublayers
        !> number of sublayers stored in the main layer
        !
        !>@param get_head_sublayer
        !> get the head_sublayer attribute
        !
        !>@param get_tail_sublayer
        !> get the tail_sublayer attribute
        !
        !>@param add_sublayer
        !> allocate space for a new buffer sublayer in the double
        !> chained list and organize the bf_mainlayer using the
        !> alignment of the buffer layers. A pointer to the newly
        !> added buffer sublayer is returned
        !
        !>@param merge_sublayers
        !> combine two sublayers of the main layer
        !
        !>@param remove_sublayer
        !> remove a sublayer from the doubled chained list
        !
        !>@param print_binary
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !>@param print_netcdf
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate netcdf output files
        !
        !> @param allocate_before_timeInt
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param deallocate_after_timeInt
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
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
          procedure, pass :: remove_sublayer

          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

          procedure, pass :: allocate_before_timeInt
          procedure, pass :: deallocate_after_timeInt
          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step
          
        end type bf_mainlayer_pointer


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the bf_mainlayer_ptr object by nullify its
        !> ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          call this%nullify_ptr()

        end subroutine ini        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@return get_ptr
        !> pointer to the ptr attribute
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@param bf_mainlayer_ptr
        !> pointer to the bf_mainlayer object
        !--------------------------------------------------------------
        subroutine set_ptr(this, bf_mainlayer_ptr)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this
          type(bf_mainlayer), pointer, intent(in)    :: bf_mainlayer_ptr

          this%ptr => bf_mainlayer_ptr

        end subroutine set_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> nullify the ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !--------------------------------------------------------------
        subroutine nullify_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          nullify(this%ptr)

        end subroutine nullify_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate space for the ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !--------------------------------------------------------------
        subroutine allocate_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          allocate(this%ptr)

        end subroutine allocate_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate space for the ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !--------------------------------------------------------------
        subroutine deallocate_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          deallocate(this%ptr)
          nullify(this%ptr)

        end subroutine deallocate_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the ptr attribute is associated
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !--------------------------------------------------------------        
        function associated_ptr(this)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          logical                                 :: associated_ptr

          associated_ptr = associated(this%ptr)

        end function associated_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the mainlayer corresponding to the
        !> ptr attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the position of the buffer
        !> sublayers stored in the mainlayer
        !--------------------------------------------------------------        
        subroutine ini_mainlayer(this, mainlayer_id)

          class(bf_mainlayer_pointer), intent(inout) :: this
          integer                    , intent(in)    :: mainlayer_id

          !debug: check mainlayer_id
          if(debug) then
             if((mainlayer_id.lt.1).or.(mainlayer_id.gt.4)) then
                call error_mainlayer_id(
     $               'bf_mainlayer_pointer_class.f',
     $               'ini_mainlayer',
     $               mainlayer_id)
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

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the mainlayer_id attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@return get_mainlayer_id
        !> cardinal coordinate identifying the position of the buffer
        !> sublayers stored in the mainlayer
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_sublayers attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@return get_nb_sublayers
        !> number of sublayers in the mainlayer chained list
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the head_sublayer attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@return get_head_sublayer
        !> pointer to the first sublayer in the chained list
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the tail_sublayer attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@return get_head_sublayer
        !> pointer to the last sublayer in the chained list
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate space for a new buffer sublayer in the double
        !> chained list and organize the bf_mainlayer using the
        !> alignment of the buffer layers. A pointer to the newly
        !> added buffer sublayer is returned
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer_ptr object encapsulating a pointer to
        !> a bf_mainlayer object
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table(2,2) of integers identifying the position of the buffer
        !> layer compared to the interior nodes
        !
        !> @return added_sublayer_ptr
        !> pointer to the newly added bf_sublayer
        !--------------------------------------------------------------
        function add_sublayer(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(added_sublayer_ptr)

          implicit none

          class(bf_mainlayer_pointer)        , intent(inout) :: this
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: alignment

          type(bf_sublayer), pointer                         :: added_sublayer_ptr

          if(this%associated_ptr()) then
             added_sublayer_ptr => this%ptr%add_sublayer(
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            alignment)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''add_sublayer'')'
             stop 'ptr attribute not associated'
          end if
          
        end function add_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> combine two sublayers of the main layer
        !
        !> @date
        !> 09_05_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param bf_sublayer1
        !> pointer to the first sublayer to merge
        !
        !>@param bf_sublayer2
        !> pointer to the second sublayer to merge
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table identifying the final position of the merged sublayer
        !> compared to the interior domain
        !
        !>@return merged_sublayer
        !> pointer to the bf_sublayer resulting from the merge of the
        !> buffer sublayers
        !--------------------------------------------------------------
        function merge_sublayers(
     $     this,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(merged_sublayer)
        
          implicit none        
        
          class(bf_mainlayer_pointer)             , intent(inout) :: this
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer2
          real(rkind)   , dimension(nx)           , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)           , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2), optional, intent(in)    :: alignment
          type(bf_sublayer), pointer                              :: merged_sublayer
        

          if(this%associated_ptr()) then
             merged_sublayer => this%ptr%merge_sublayers(
     $            bf_sublayer1,
     $            bf_sublayer2,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            alignment)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''merge_sublayers'')'
             stop 'ptr attribute not associated'
          end if

        end function merge_sublayers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove a sublayer from the doubled chained list
        !
        !> @date
        !> 09_05_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param bf_sublayer_ptr
        !> pointer to the sublayer to be removed
        !--------------------------------------------------------------
        subroutine remove_sublayer(
     $     this,
     $     sublayer_ptr)
        
          implicit none        
        
          class(bf_mainlayer_pointer), intent(inout) :: this
          type(bf_sublayer), pointer , intent(inout) :: sublayer_ptr
        

          if(this%associated_ptr()) then
             call this%ptr%remove_sublayer(sublayer_ptr)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''remove_sublayer'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine remove_sublayer

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !> @date
        !> 09_05_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param suffix_nodes
        !> suffix for the name of the output files storing the nodes
        !> of the bf_sublayers
        !
        !>@param suffix_grdid
        !> suffix for the name of the output files storing the grdpts_id
        !> of the bf_sublayers
        !
        !>@param suffix_sizes
        !> suffix for the name of the output files storing the sizes
        !> of the bf_sublayers        
        !--------------------------------------------------------------
        subroutine print_binary(
     $     this,
     $     suffix_x_map,
     $     suffix_y_map,
     $     suffix_nodes,
     $     suffix_grdid,
     $     suffix_sizes)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          character(*)               , intent(in) :: suffix_x_map
          character(*)               , intent(in) :: suffix_y_map
          character(*)               , intent(in) :: suffix_nodes
          character(*)               , intent(in) :: suffix_grdid
          character(*)               , intent(in) :: suffix_sizes


          if(this%associated_ptr()) then
             call this%ptr%print_binary(
     $            suffix_x_map,
     $            suffix_y_map,
     $            suffix_nodes,
     $            suffix_grdid,
     $            suffix_sizes)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''merge_sublayers'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine print_binary


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on a netcdf output file
        !
        !> @date
        !> 11_07_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param timestep_written
        !> integer identifying the timestep written
        !
        !>@param name_var
        !> table with the short name for the governing variables saved
        !> in the netcdf file
        !
        !>@param longname_var
        !> table with the long name for the governing variables saved
        !> in the netcdf file
        !
        !>@param unit_var
        !> table with the units of the governing variables saved
        !> in the netcdf file
        !
        !>@param x_min_interior
        !> x-coordinate of interior grid point next to the left
        !> boundary layer
        !
        !>@param y_min_interior
        !> y-coordinate of interior grid point next to the lower
        !> boundary layer
        !
        !>@param dx
        !> grid size along the x-coordinate
        !
        !>@param dy
        !> grid size along the y-coordinate
        !
        !>@param time
        !> time corresponding to the data for the grdpts and the nodes
        !--------------------------------------------------------------
        subroutine print_netcdf(
     $     this,
     $     timestep_written,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     time)

          implicit none

          class(bf_mainlayer_pointer), intent(in) :: this
          integer                    , intent(in) :: timestep_written
          character(*), dimension(ne), intent(in) :: name_var
          character(*), dimension(ne), intent(in) :: longname_var
          character(*), dimension(ne), intent(in) :: unit_var
          real(rkind)                , intent(in) :: time


          if(this%associated_ptr()) then
             call this%ptr%print_netcdf(
     $            timestep_written,
     $            name_var,
     $            longname_var,
     $            unit_var,
     $            time)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''print_netcdf'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine print_netcdf


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine allocate_before_timeInt(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          if(this%associated_ptr()) then
             call this%ptr%allocate_before_timeInt()
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''allocate_before_timeInt'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine allocate_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine deallocate_after_timeInt(this)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this

          if(this%associated_ptr()) then
             call this%ptr%deallocate_after_timeInt()
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''deallocate_after_timeInt'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine deallocate_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the sublayers
        !> contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this
          type(td_operators)         , intent(in)    :: td_operators_used
          real(rkind)                , intent(in)    :: t
          type(sd_operators)         , intent(in)    :: s
          type(pmodel_eq)            , intent(in)    :: p_model
          type(bc_operators)         , intent(in)    :: bc_used

          if(this%associated_ptr()) then
             call this%ptr%compute_time_dev(
     $            td_operators_used,
     $            t,s,p_model,bc_used)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''compute_time_dev'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param dt
        !> integration time step
        !
        !>@param integration_step_nopt
        !> procedure performing the time integration
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt)

          implicit none

          class(bf_mainlayer_pointer), intent(inout) :: this
          real(rkind)                , intent(in)    :: dt
          procedure(timeInt_step_nopt) :: integration_step_nopt

          if(this%associated_ptr()) then
             call this%ptr%compute_integration_step(
     $            dt, integration_step_nopt)
          else
             print '(''bf_mainlayer_pointer_class'')'
             print '(''compute_integration_step'')'
             stop 'ptr attribute not associated'
          end if

        end subroutine compute_integration_step        

      end module bf_mainlayer_pointer_class
