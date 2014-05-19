      module bf_interface_class

        use bf_sublayer_class         , only : bf_sublayer
        use bf_mainlayer_class        , only : bf_mainlayer
        use bf_mainlayer_pointer_class, only : bf_mainlayer_pointer

        use parameters_bf_layer       , only : clockwise, counter_clockwise
        use parameters_constant       , only : N,S,E,W,
     $                                         x_direction, y_direction,
     $                                         interior
        use parameters_input          , only : nx,ny,ne,bc_size,debug
        use parameters_kind           , only : ikind, rkind


        implicit none

        private
        public :: bf_interface
        

        !> @class bf_interface
        !> class encapsulating the bf_layer/interior interface
        !> object
        !
        !> @param mainlayer_pointers
        !> table with reference to the buffer main layers
        !
        !> @param ini
        !> initialize the interface by nullifying all the 
        !> references to the buffer main layers
        !
        !> @param get_mainlayer
        !> get the reference to the main layer corresponding
        !> to the cardinal coordinate passed
        !
        !> @param add_sublayer
        !> add a new sublayer to the main layer corresponding
        !> to the cardinal coordinate passed
        !---------------------------------------------------------------
        type :: bf_interface

          type(bf_mainlayer_pointer), dimension(4) :: mainlayer_pointers

          contains

          procedure, pass :: ini
          procedure, pass :: get_mainlayer
          procedure, pass :: add_sublayer

          procedure, nopass :: get_mainlayer_id
          procedure, pass   :: get_sublayer
          procedure, nopass :: get_neighboring_sublayer

          procedure, pass :: print_binary

        end type bf_interface

        contains


        !< nullify all the pointers to the main layers
        subroutine ini(this)

          implicit none

          class(bf_interface), intent(inout) :: this

          integer :: i
          
          do i=1, size(this%mainlayer_pointers,1)
             call this%mainlayer_pointers(i)%ini()             
          end do

        end subroutine ini


       !< get main layer corresponding to the cardinal point
       function get_mainlayer(this, mainlayer_id)
       
         implicit none
   
         class(bf_interface), intent(in) :: this
         integer            , intent(in) :: mainlayer_id
         type(bf_mainlayer) , pointer    :: get_mainlayer

         if(debug) then
            if((mainlayer_id.lt.1).or.(mainlayer_id.gt.8)) then
               print '(''bf_interface_class'')'
               print '(''get_mainlayer'')'
               print '(''mainlayer_id not recognized'')'
               print '(''mainlayer_id: '',I2)', mainlayer_id
               stop 'change mainlayer_id: N,S,E,W'
            end if
         end if


         if(this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
         else
            nullify(get_mainlayer)
         end if


       end function get_mainlayer


       !< add a buffer sublayer to the existing layers in the corresponding
       !> buffer mainlayer
       function add_sublayer(
     $     this,
     $     mainlayer_id,
     $     nodes,
     $     alignment)
     $     result(added_sublayer)
        
          class(bf_interface)             , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer, dimension(2,2)         , intent(in)    :: alignment

          type(bf_sublayer), pointer                      :: added_sublayer


          !debug : check the mainlayer id
          if(debug) then
             if((mainlayer_id.lt.1).or.(mainlayer_id.gt.8)) then
                print '(''bf_interface_class'')'
                print '(''add_sublayer'')'
                print '(''mainlyer_id not recognized'')'
                print '(''mainlayer_id: '',I2)', mainlayer_id
                stop 'change mainlayer_id'
             end if
          end if

          !first check if the mainlayer corresponding to the cardinal
          !point is indeed allocated: if the memory space is not allocated,
          !the space in memory is first allocated, the pointer identifying
          !the mainlayer is initialized and the main layer itself is initialized
          if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
             call this%mainlayer_pointers(mainlayer_id)%ini_mainlayer(mainlayer_id)
          end if            

          !now that we are sure that space is allocated for the main layer,
          !the sublayer can be integrated to the mainlayer and the buffer
          !layer can be initialized using the nodes, alignment and neighbors
          !arguments
          added_sublayer   => this%mainlayer_pointers(mainlayer_id)%add_sublayer(
     $         nodes, alignment)

       end function add_sublayer


       !> get neighbor of the same mainlayer
       function get_neighboring_sublayer(current_bf_sublayer)
     $     result(neighboring_sublayer)

         implicit none

         type(bf_sublayer), pointer, intent(in) :: current_bf_sublayer
         type(bf_sublayer), pointer             :: neighboring_sublayer

         if(associated(current_bf_sublayer%get_next())) then
            neighboring_sublayer => current_bf_sublayer%get_next()
         else
            nullify(neighboring_sublayer)
         end if

       end function get_neighboring_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> subroutine converting general coordinates into
       !> the main layer ID (N,S,E,W)
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
            mainlayer_id = S

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
               mainlayer_id = N

            end if
         end if

       end function get_mainlayer_id


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
     $    this,
     $    general_coord,
     $    local_coord,
     $    tolerance_i,
     $    mainlayer_id_i)
     $    result(sublayer)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         integer(ikind), dimension(2), intent(in)  :: general_coord
         integer(ikind), dimension(2), intent(out) :: local_coord
         integer       , optional    , intent(in)  :: tolerance_i
         integer       , optional    , intent(in)  :: mainlayer_id_i
         type(bf_sublayer), pointer                :: sublayer


         integer                     :: direction_tested
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
            if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
               nullify(sublayer)
                 
            !< if the main layer exists, the sublayers saved inside are
            !< checked to decide whether the grid point asked belongs to
            !< one of them or not
            else
               mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
            
               !< check if sublayers are saved inside the mainlayer
               !> if no sublayers are saved inside the mainlayer,
               !> no existing sublayer can match the general coord
               !> and so the pointer to sublayer is nullified
               if(.not.associated(mainlayer%get_head_sublayer())) then
                  nullify(sublayer)
            
               !< otherwise, the sublayer corresponding to the general
               !> coordinates is searched by going through the different
               !> element of the doubled chained list
               else
                  sublayer => mainlayer%get_head_sublayer()
            
              	  !< processing the tolerance for matching a sublayer
            	  !> if no tolerence is provided, the default option is 0
                  if(.not.present(tolerance_i)) then
                     tolerance=0
                  else
                     tolerance=tolerance_i
                  end if
                  
                  !< if the mainlayer investigated is N,S,E or W, there can
                  !> be sublayers to be investigated
                  select case(mainlayer_id)
                    case(N,S)
                       direction_tested = x_direction
                    case(E,W)
                       direction_tested = y_direction
                  end select 

                  !check if the grid point belongs to the current sublayer
                  grdpt_in_sublayer =
     $                 (  general_coord(direction_tested).ge.
     $                 (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                 .and.(
     $                 general_coord(direction_tested).le.
     $                 (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  !go through the different sublayers
                  do while(.not.grdpt_in_sublayer)
            	        
            	     !if no matching sublayer can be found
            	     !nullify the corresponding pointer
                     if(.not.associated(sublayer%get_next())) then
                        nullify(sublayer)
                        exit
                     end if
                   
                     sublayer => sublayer%get_next()
                     grdpt_in_sublayer =
     $                    (general_coord(direction_tested).ge.
     $                    (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                    .and.(
     $                    general_coord(1).le.
     $                    (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  end do
            
               end if
            end if
         end if

         !< if a sublayer matching the general coordinates was found
         !> compute the local coordinates in this sublayer
         if(associated(sublayer)) then
            local_coord = sublayer%get_local_coord(general_coord)
         end if           
         
       end function get_sublayer


       !< print the content of the interface on external binary files
       subroutine print_binary(
     $     this,
     $     suffix_nodes, suffix_grdid, suffix_sizes,
     $     suffix_nb_sublayers_max)

         implicit none

         class(bf_interface), intent(in) :: this
         character(*)       , intent(in) :: suffix_nodes
         character(*)       , intent(in) :: suffix_grdid
         character(*)       , intent(in) :: suffix_sizes
         character(*)       , intent(in) :: suffix_nb_sublayers_max
         

         integer           :: i
         integer           :: nb_sublayers_max
         character(len=18) :: filename_format
         character(len=28) :: nb_sublayers_filename
                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files and 
         !determine the maximum number of sublayers
         nb_sublayers_max = 0
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_binary(suffix_nodes,
     $                                            suffix_grdid,
     $                                            suffix_sizes)

               nb_sublayers_max = max(
     $              nb_sublayers_max,
     $              this%mainlayer_pointers(i)%get_nb_sublayers())
            end if            
         end do


         !print the maximum number of sublayers in an output
         !binary file
         write(filename_format,
     $        '(''(A12,A'',I2,'')'')')
     $        len(suffix_nb_sublayers_max)

         write(nb_sublayers_filename, filename_format)
     $        'sublayers_nb',
     $        suffix_nb_sublayers_max

         call print_nb_sublayers_max(
     $        nb_sublayers_filename, nb_sublayers_max)

        end subroutine print_binary


        subroutine print_nb_sublayers_max(filename, nb_sublayers)

          implicit none

          character(*), intent(in) :: filename
          integer     , intent(in) :: nb_sublayers

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers_max

      end module bf_interface_class
