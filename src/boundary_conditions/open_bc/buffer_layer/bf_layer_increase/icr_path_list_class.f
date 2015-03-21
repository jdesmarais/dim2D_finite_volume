      !> @file
      !> double chained list of icr_path_chain elements
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> double chained list of icr_path_chain elements
      !
      !> @date
      ! 21_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_path_list_class

        use icr_path_chain_class, only :
     $     icr_path_chain

        implicit none

        private
        public :: icr_path_list


        !> @class icr_path_list
        !> double chained list of icr_path_chain elements
        !
        !> @param nb_paths
        !> number of paths stored in the double chain list
        !
        !> @param head_path
        !> pointer to the first path of the list
        !
        !> @param tail_path
        !> pointer to the last path of the list
        !
        !> @param ini
        !> initialize the path list by initializing the
        !> number of paths to 0 and nullifying the head
        !> and tail attributes
        !
        !> @param get_nb_paths
        !> number of paths stored in the main layer
        !
        !> @param get_head_path
        !> get the head_path attribute
        !
        !> @param get_tail_path
        !> get the tail_path attribute
        !
        !> @param add_path
        !> allocate space for a new path in the double
        !> chained list. A pointer to the newly added
        !> path is returned
        !
        !> @param remove_path
        !> remove a path from the doubled chained list
        !---------------------------------------------------------------
        type :: icr_path_list

          integer :: nb_paths

          type(icr_path_chain), pointer :: head_path
          type(icr_path_chain), pointer :: tail_path

          contains

          !initialization
          procedure, pass :: ini

          !get attributes
          procedure, pass :: get_nb_paths
          procedure, pass :: get_head_path
          procedure, pass :: get_tail_path

          !add/remove a path
          procedure, pass :: add_path
          procedure, pass :: remove_path

        end type icr_path_list

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the path list by initializing the
        !> number of paths to 0 and nullifying the head
        !> and tail attributes
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> double chained list of icr_path_chain elements
        !
        !> @param mainlayer_id
        !> cardinal coordinate identifying the position of the buffer
        !> sublayers stored in the mainlayer
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(icr_path_list), intent(inout) :: this

          this%nb_paths = 0
          nullify(this%head_path)
          nullify(this%tail_path)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_paths attribute
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of paths, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@return get_nb_paths
        !> number of paths in the mainlayer chained list
        !--------------------------------------------------------------
        function get_nb_paths(this)

          implicit none

          class(icr_path_list), intent(in) :: this
          integer                          :: get_nb_paths

          get_nb_paths = this%nb_paths

        end function get_nb_paths


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the head_path attribute
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> double chained list of icr_path_chain elements
        !
        !>@return get_head_path
        !> pointer to the first path in the chained list
        !--------------------------------------------------------------
        function get_head_path(this)

          implicit none

          class(icr_path_list), intent(in) :: this
          type(icr_path_chain), pointer    :: get_head_path

          if(associated(this%head_path)) then
             get_head_path => this%head_path
          else
             nullify(get_head_path)
          end if

        end function get_head_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the tail_path attribute
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> double chained list of icr_path_chain elements
        !
        !>@return get_head_path
        !> pointer to the last path in the chained list
        !--------------------------------------------------------------
        function get_tail_path(this)

          implicit none

          class(icr_path_list), intent(in) :: this
          type(icr_path_chain), pointer    :: get_tail_path

          if(associated(this%tail_path)) then
             get_tail_path => this%tail_path
          else
             nullify(get_tail_path)
          end if

        end function get_tail_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate space for a new buffer path in the double
        !> chained list. A pointer to the newly
        !> added buffer path is returned
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> double chained list of icr_path_chain elements
        !
        !> @return added_path_ptr
        !> pointer to the newly added icr_path_chain
        !--------------------------------------------------------------
        function add_path(
     $     this)
     $     result(added_path_ptr)

          implicit none
          
          class(icr_path_list), intent(inout) :: this
          type(icr_path_chain), pointer       :: added_path_ptr


          !allocate space for the path and initialize it
          allocate(added_path_ptr)
          call added_path_ptr%ini()

          !insert the path at the end of the double chained list
          if(associated(this%tail_path)) then
             call this%tail_path%set_next(added_path_ptr)
             call added_path_ptr%set_prev(this%tail_path)
             this%tail_path => added_path_ptr

          else
             this%head_path => added_path_ptr
             this%tail_path => added_path_ptr

          end if

          !update the number of paths in the mainlayer
          this%nb_paths = this%nb_paths+1

        end function add_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove a path from the doubled chained list
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> double chained list of icr_path_chain elements
        !
        !> @param path_ptr
        !> pointer to the path to be removed
        !--------------------------------------------------------------
        subroutine remove_path(this, path_ptr)

          implicit none

          class(icr_path_list)         , intent(inout) :: this
          type(icr_path_chain), pointer, intent(inout) :: path_ptr


          !check if the head of the mainlayer should be changed
          if(associated(this%head_path,path_ptr)) then
             if(associated(path_ptr%get_next())) then
                this%head_path => path_ptr%get_next()
             else
                nullify(this%head_path)
             end if
          end if


          !check if the tail of the mainlayer should be changed
          if(associated(this%tail_path, path_ptr)) then
             if(associated(path_ptr%get_prev())) then
                this%tail_path => path_ptr%get_prev()
             else
                nullify(this%tail_path)
             end if
          end if


          !remove the buffer layer
          call path_ptr%remove()
          deallocate(path_ptr)
          nullify(path_ptr)


          !update the number of paths in the mainlayer
          this%nb_paths = this%nb_paths-1

        end subroutine remove_path

      end module icr_path_list_class
