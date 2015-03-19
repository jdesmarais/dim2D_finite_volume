      !> @file
      !> bf_mainlayer_time augmented with procedures to
      !> extract sublayer using coordinates
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_mainlayer_time augmented with procedures to
      !> extract sublayer using coordinates
      !
      !> @date
      ! 19_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_class

        use bf_increase_coords_module, only :
     $       get_mainlayer_coord

        use bf_mainlayer_time_class, only :
     $       bf_mainlayer_time

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: bf_mainlayer
        
        
        !> @class bf_mainlayer
        !> bf_mainlayer_time augmented with procedures to
        !> extract sublayer using coordinates
        !
        !> @param get_sublayer_from_coords
        !---------------------------------------------------------------
        type, extends(bf_mainlayer_time) :: bf_mainlayer

          contains

          procedure, pass :: get_sublayer_from_gen_coords

        end type bf_mainlayer


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the sublayer whose coordinates match the
        !> buffer layer inside the mainlayer
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param gen_coords
        !> general coordinates identifying the position of the grid
        !> point in the general reference frame
        !
        !> @param tolerance
        !> integer identifying the distance allowed between the grid
        !> point and the buffer layer such that it is considered inside
        !> (along the x-direction for N and S buffer layers and along the
        !> y-direction for the E and W buffer layers), 0 by default
        !
        !> @param no_check_ID
        !> logical determining whether the mainlayer_id of the grid-point
        !> is computed to see if it matches the mainlayer_interface
        !> coordinate, .false. by default
        !--------------------------------------------------------------
        function get_sublayer_from_gen_coords(
     $       this,
     $       gen_coords,
     $       tolerance,
     $       no_check_ID)
     $       result(bf_sublayer_ptr)

          implicit none
          
          class(bf_mainlayer)         , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: gen_coords
          integer(ikind), optional    , intent(in) :: tolerance
          logical       , optional    , intent(in) :: no_check_ID
          type(bf_sublayer), pointer               :: bf_sublayer_ptr

          integer                    :: tolerance_op
          logical                    :: no_check_ID_op
          integer                    :: mainlayer_id
          integer                    :: dir
          integer(ikind)             :: min_i
          integer(ikind)             :: max_i
          integer                    :: k
          type(bf_sublayer), pointer :: current_sublayer
          

          if(present(tolerance)) then
             tolerance_op = tolerance
          else
             tolerance_op = 0
          end if

          if(present(no_check_ID)) then
             no_check_ID_op = no_check_ID
          else
             no_check_ID_op = .false.
          end if


          if(.not.no_check_ID_op) then
             mainlayer_id = get_mainlayer_coord(gen_coords)
             if(mainlayer_id.ne.this%mainlayer_id) then
                print '(''bf_mainlayer_class'')'
                print '(''get_sublayer_from_coords'')'
                print '(''the mainlayer ID do not match:'')'
                print '(I2,''->'',I2)', mainlayer_id, this%mainlayer_id
                stop ''
             end if
          end if


          select case(this%mainlayer_id)
            case(N,S)
               dir = x_direction
            case(E,W)
               dir = y_direction
            case default
               call error_mainlayer_id(
     $              'bf_mainlayer_class',
     $              'get_sublayer_from_coords',
     $              this%mainlayer_id)
          end select


          nullify(bf_sublayer_ptr)


          current_sublayer => this%head_sublayer

          do k=1, this%nb_sublayers

             min_i = current_sublayer%get_alignment(dir,1) - bc_size
             max_i = current_sublayer%get_alignment(dir,2) + bc_size

             !if there is some overlap between
             ![gen_coords(dir)-tolerance_op,gen_coords(dir)+tolerance_op]
             !and [min_i,max_i], the grid point belongs to the buffer
             !layer
             if((
     $            min(gen_coords(dir)+tolerance_op,max_i)-
     $            max(gen_coords(dir)-tolerance_op,min_i)+1).gt.0) then

                bf_sublayer_ptr => current_sublayer
                exit

             end if

             current_sublayer => current_sublayer%get_next()

          end do

        end function get_sublayer_from_gen_coords

      end module bf_mainlayer_class
