      module bf_interface_class

        use bf_sublayer_class         , only : bf_sublayer
        use bf_mainlayer_class        , only : bf_mainlayer
        use bf_mainlayer_pointer_class, only : bf_mainlayer_pointer

        use parameters_constant       , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input          , only : nx,ny,ne, debug
        use parameters_kind           , only : rkind


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

          type(bf_mainlayer_pointer), dimension(8) :: mainlayer_pointers

          contains

          procedure, pass :: ini
          procedure, pass :: get_mainlayer
          procedure, pass :: add_sublayer

          procedure, pass :: get_neighboring_sublayers

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
               stop 'change mainlayer_id: N,S,E,W,N_E,N_W,S_E,S_W'
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
     $     alignment,
     $     neighbors)
     $     result(added_sublayer)
        
          class(bf_interface)             , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer, dimension(2,2)         , intent(in)    :: alignment
          logical, dimension(4)           , intent(in)    :: neighbors
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
     $         nodes, alignment, neighbors)

       end function add_sublayer


       subroutine get_neighboring_sublayers(
     $     this, corner_id,
     $     neighboring_sublayer1, neighboring_sublayer2)

         implicit none

         class(bf_interface), intent(in) :: this
         integer            , intent(in) :: corner_id
         type(bf_sublayer)  , pointer    :: neighboring_sublayer1
         type(bf_sublayer)  , pointer    :: neighboring_sublayer2

         
         select case(corner_id)
           case(N_E)
              if(this%mainlayer_pointers(N)%associated_ptr()) then
                 neighboring_sublayer1 =>
     $                this%mainlayer_pointers(N)%get_tail_sublayer()
              else
                 nullify(neighboring_sublayer1)
              end if

              if(this%mainlayer_pointers(E)%associated_ptr()) then
                 neighboring_sublayer2 =>
     $                this%mainlayer_pointers(E)%get_tail_sublayer()
              else
                 nullify(neighboring_sublayer2)
              end if

           case(N_W)
              if(this%mainlayer_pointers(W)%associated_ptr()) then
                 neighboring_sublayer1 =>
     $                this%mainlayer_pointers(W)%get_tail_sublayer()
              else
                 nullify(neighboring_sublayer1)
              end if

              if(this%mainlayer_pointers(N)%associated_ptr()) then
                 neighboring_sublayer2 =>
     $                this%mainlayer_pointers(N)%get_head_sublayer()
              else
                 nullify(neighboring_sublayer2)
              end if

           case(S_E)
              if(this%mainlayer_pointers(E)%associated_ptr()) then
                 neighboring_sublayer1 =>
     $                this%mainlayer_pointers(E)%get_head_sublayer()
              else
                 nullify(neighboring_sublayer1)
              end if

              if(this%mainlayer_pointers(S)%associated_ptr()) then
                 neighboring_sublayer2 =>
     $                this%mainlayer_pointers(S)%get_tail_sublayer()
              else
                 nullify(neighboring_sublayer2)
              end if

           case(S_W)
              if(this%mainlayer_pointers(S)%associated_ptr()) then
                 neighboring_sublayer1 =>
     $                this%mainlayer_pointers(S)%get_head_sublayer()
              else
                 nullify(neighboring_sublayer1)
              end if

              if(this%mainlayer_pointers(W)%associated_ptr()) then
                 neighboring_sublayer2 =>
     $                this%mainlayer_pointers(W)%get_head_sublayer()
              else
                 nullify(neighboring_sublayer2)
              end if

           case default
              print '(''bf_interface_class'')'
              print '(''get_neighboring_mainlayers'')'
              print '(''corner_id not recognized'')'
              print '(''corner_id: '',I2)', corner_id
              stop 'check corner_id asked'
         end select

       end subroutine get_neighboring_sublayers
       


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
