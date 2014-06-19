      module test_case_class

        use bf_interface_class    , only : bf_interface
        use bf_interface_icr_class, only : bf_interface_icr
        use bf_mainlayer_class    , only : bf_mainlayer
        use bf_sublayer_class     , only : bf_sublayer
        use bubble_class          , only : bubble
        use rising_bubble_class   , only : rising_bubble
        use growing_bubble_class  , only : growing_bubble
        use parameters_bf_layer   , only : no_pt
        use parameters_input      , only : nx, ny, ne, bc_size
        use parameters_kind       , only : ikind, rkind
        use test_bf_layer_module  , only : print_interior_data
        

        implicit none


        private
        public :: test_case


        type :: test_case

          class(bubble), allocatable :: bubble_used
          integer                    :: vf_choice
          real(rkind)                :: angle

          contains
          
          procedure,   pass          :: ini
          procedure,   pass          :: ini_with_set
          procedure,   pass          :: set_bubble
          procedure, nopass, private :: get_c_coords
          procedure,   pass, private :: get_r_coords
          procedure,   pass          :: get_coords
          procedure,   pass, private :: get_mass
          procedure,   pass, private :: get_velocity
          procedure,   pass          :: get_nodes
          procedure,   pass, private :: update_interior_nodes
          procedure,   pass, private :: update_bf_nodes
          procedure,   pass, private :: update
          procedure,   pass          :: update_nodes          
          procedure, nopass          :: print_state

        end type test_case


        contains


        subroutine ini(this, vf_choice, angle)

          implicit none

          class(test_case)          , intent(inout) :: this
          integer         , optional, intent(in)    :: vf_choice
          real(rkind)     , optional, intent(in)    :: angle

          if(present(vf_choice)) then
             this%vf_choice = vf_choice
          else
             this%vf_choice = 1
          end if

          if(present(angle)) then
             this%angle = angle
          else
             this%angle = 0.0d0
          end if

        end subroutine ini


        subroutine ini_with_set(this, set_index, nb_bubbles)

          implicit none

          class(test_case) , intent(inout) :: this
          integer          , intent(in)    :: set_index
          integer, optional, intent(in)    :: nb_bubbles

          
          type(rising_bubble)  :: rising_bubble_used
          type(growing_bubble) :: growing_bubble_used
          integer              :: nb_bubbles_i


          !rising bubble initialization
          call rising_bubble_used%ini()

          if(present(nb_bubbles)) then
             nb_bubbles_i = nb_bubbles
          else
             nb_bubbles_i = 1
          end if

          call rising_bubble_used%set_nb_bubbles(nb_bubbles_i)
          
          if(nb_bubbles_i.eq.2) then
             call rising_bubble_used%set_center([-0.1d0,0.0d0])
          end if


          !growing bubble initialization
          call growing_bubble_used%ini()
          

          !test case initialization
          select case(set_index)
            case(1)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 1
               this%angle = 0.0d0
            case(2)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 2
               this%angle = 0.0d0
            case(3)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 3
               this%angle = 0.0d0
            case(4)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 4
               this%angle = 0.0d0
            case(5)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 1
               this%angle =  ACOS(-1.0d0)/4.0d0
            case(6)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 1
               this%angle = -ACOS(-1.0d0)/4.0d0
            case(7)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 2
               this%angle =  ACOS(-1.0d0)/4.0d0
            case(8)
               call this%set_bubble(rising_bubble_used)
               this%vf_choice = 2
               this%angle = -ACOS(-1.0d0)/4.0d0
            case(9)
               call this%set_bubble(growing_bubble_used)
               this%vf_choice = 5
               this%angle = 0.0d0
            case default
               stop 'test case : ini_with_set : caes not recognized'
          end select
          
        end subroutine ini_with_set


        subroutine set_bubble(this, bubble_used)
          
          implicit none

          class(test_case), intent(inout) :: this
          class(bubble)   , intent(in)    :: bubble_used

          allocate(this%bubble_used, source=bubble_used)

        end subroutine set_bubble


        function get_c_coords(g_coords,dx,dy) result(c_coords)
        
          implicit none

          integer(ikind), dimension(2), intent(in) :: g_coords
          real(rkind)                 , intent(in) :: dx
          real(rkind)                 , intent(in) :: dy
          real(rkind)   , dimension(2)             :: c_coords

          c_coords(1) = (g_coords(1)-nx/2-3/2)*dx
          c_coords(2) = (g_coords(2)-ny/2-3/2)*dy

        end function get_c_coords



        function get_r_coords(this, c_coords, inv) result(r_coords)
        
          implicit none

          class(test_case)            , intent(in) :: this
          real(rkind)   , dimension(2), intent(in) :: c_coords
          logical       , optional    , intent(in) :: inv
          real(rkind)   , dimension(2)             :: r_coords

          real(rkind) :: angle
          
          if(present(inv)) then
             if(inv) then
                angle = -this%angle
             else
                angle = this%angle
             end if
          else
             angle = this%angle
          end if

          r_coords(1) = c_coords(1)*cos(angle) + c_coords(2)*sin(angle)
          r_coords(2) = c_coords(2)*cos(angle) - c_coords(1)*sin(angle)

        end function get_r_coords


        function get_coords(this, g_coords,dx,dy) result(r_coords)

          implicit none

          class(test_case)            , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: g_coords
          real(rkind)                 , intent(in) :: dx
          real(rkind)                 , intent(in) :: dy
          real(rkind)   , dimension(2)             :: r_coords

          real(rkind), dimension(2) :: c_coords

          c_coords = this%get_c_coords(g_coords,dx,dy)
          r_coords = this%get_r_coords(c_coords)

        end function get_coords


        function get_mass(this,coords) result(mass)

          implicit none

          class(test_case)              , intent(in) :: this
          real(rkind)     , dimension(2), intent(in) :: coords
          real(rkind)                                :: mass

          mass = this%bubble_used%get_mass(coords)

        end function get_mass


        function get_velocity(this,coords,r_ref) result(velocity)

          implicit none

          class(test_case), intent(in) :: this
          real(rkind)     , dimension(2), intent(in) :: coords
          logical             , optional, intent(in) :: r_ref
          real(rkind)     , dimension(2)             :: velocity

          select case(this%vf_choice)
            case(1)
               velocity(1) =  1.0d0
               velocity(2) =  0.0d0
            case(2)
               velocity(1) = -1.0d0
               velocity(2) =  0.0d0
            case(3)
               velocity(1) =  0.0d0
               velocity(2) =  1.0d0
            case(4)
               velocity(1) =  0.0d0
               velocity(2) = -1.0d0
            case(5)
               if((coords(1).eq.0).and.(coords(2).eq.0)) then
                  velocity(1) = 0.0d0
                  velocity(2) = 0.0d0
               else
                  velocity(1) = coords(1)/(SQRT(coords(1)**2+coords(2)**2))
                  velocity(2) = coords(2)/(SQRT(coords(1)**2+coords(2)**2))
               end if
            case default
               stop 'test case not recognized'
          end select

          if(present(r_ref)) then
             if(.not.r_ref) then
                velocity = this%get_r_coords(velocity,inv=.true.)
             end if
          else
             velocity = this%get_r_coords(velocity,.true.)
          end if

        end function get_velocity


        function get_nodes(this, g_coords, dx, dy) result(nodes)

          implicit none

          class(test_case)             , intent(in) :: this
          integer(ikind), dimension(2) , intent(in) :: g_coords
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)   , dimension(ne)             :: nodes

          real(rkind), dimension(2) :: coords
          real(rkind)               :: mass
          real(rkind), dimension(2) :: velocity

          coords   = this%get_coords(g_coords,dx,dy)
          mass     = this%get_mass(coords)
          velocity = this%get_velocity(coords)

          nodes(1) = mass
          nodes(2) = mass*velocity(1)
          nodes(3) = mass*velocity(2)          
          
        end function get_nodes

      
        !< update the data in the system
        subroutine update(this, dx, dy)

          implicit none

          class(test_case), intent(inout) :: this
          real(rkind)     , intent(in)    :: dx
          real(rkind)     , intent(in)    :: dy

          real(rkind), dimension(2) :: velocity

          if(allocated(this%bubble_used)) then
             velocity = this%get_velocity([1.0d0,1.0d0],r_ref=.true.)
             call this%bubble_used%update(velocity, dx, dy)
          else
             stop 'test_case : update: bubble not allocated'
          end if
          

        end subroutine update


        !< write data in the interior nodes to simulate a
        !> vapor bubble moving in the interior domain
        subroutine update_interior_nodes(this, dx, dy, interior_nodes)

          implicit none

          class(test_case)                , intent(in)  :: this
          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes


          integer(ikind) :: i,j


          !initialization of the nodes coresponding
          !to the bubble
          do j=1, ny
             do i=1, nx

                interior_nodes(i,j,:) = this%get_nodes([i,j],dx,dy)

             end do
          end do

        end subroutine update_interior_nodes


        !< write data in the buffer layer nodes to simulate a
        !> vapor bubble moving in the buffer layers
        subroutine update_bf_nodes(this, dx, dy, interface_used)

          implicit none

          class(test_case)   , intent(in)    :: this
          real(rkind)        , intent(in)    :: dx
          real(rkind)        , intent(in)    :: dy
          class(bf_interface), intent(inout) :: interface_used


          integer        :: k,l
          integer(ikind) :: i,j

          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr

          integer                                       :: nb_sublayers
          integer(ikind), dimension(2,2)                :: alignment_tab
          integer(ikind), dimension(2)                  :: sizes
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer(ikind), dimension(:,:)  , allocatable :: grdpts_id
          integer(ikind), dimension(2)                  :: g_coords
          

          !loop over the main layers
          do k=1,4

             !get the mainlayer
             mainlayer_ptr => interface_used%get_mainlayer(k)
             if(associated(mainlayer_ptr)) then

                nb_sublayers =  mainlayer_ptr%get_nb_sublayers()
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()

                !get the sublayers of the mainlayer
                do l=1, nb_sublayers

                   !get the alignment
                   alignment_tab = sublayer_ptr%get_alignment_tab()

                   !get the sizes of the nodes table
                   sizes = sublayer_ptr%get_sizes()

                   !allocate the new nodes
                   allocate(new_nodes(sizes(1), sizes(2), ne))

                   !get the grdpts_id
                   call sublayer_ptr%get_grdpts_id(grdpts_id)
                   
                   !fill the nodes
                   do j=1, sizes(2)
                      do i=1, sizes(1)
                         if(grdpts_id(i,j).ne.no_pt) then
                            g_coords(1) = alignment_tab(1,1)-bc_size+(i-1)
                            g_coords(2) = alignment_tab(2,1)-bc_size+(j-1)

                            new_nodes(i,j,:) = this%get_nodes(g_coords,dx,dy)

                         else
                            new_nodes(i,j,1)=0.0d0
                            new_nodes(i,j,2)=0.0d0
                            new_nodes(i,j,3)=0.0d0
                         end if
                      end do
                   end do

                   !replace the content of the nodes
                   call sublayer_ptr%set_nodes(new_nodes)

                   !deallocate the temporary table for the grid points
                   deallocate(grdpts_id)

                   !next sublayer in the main layer
                   sublayer_ptr => sublayer_ptr%get_next()

                end do

             end if

          end do

        end subroutine update_bf_nodes


        !< update the data in the interior nodes and the buffer layers
        !> to simulate a vapor bubble moving in the computational domain
        subroutine update_nodes(
     $       this,
     $       dx, dy,
     $       interior_nodes, interface_used)


          implicit none

          class(test_case)                , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne), intent(out)   :: interior_nodes
          class(bf_interface)             , intent(inout) :: interface_used

          !update the interior nodes
          call this%update_interior_nodes(dx, dy, interior_nodes)

          !update the data in the buffer layers
          call this%update_bf_nodes(dx, dy, interface_used)

          !update the system
          call this%update(dx,dy)

        end subroutine update_nodes


        !> print the state of the system
        subroutine print_state(nodes, grdpts_id, interface_used, index)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface_icr)         , intent(in)    :: interface_used
          integer                         , intent(inout) :: index

          integer :: format_index

          character(len=11) :: i_format_nodes
          character(len=15) :: i_format_grdpt

          character(len=15) :: bf_format_nodes
          character(len=18) :: bf_format_grdpt
          character(len=10) :: bf_format_nbsbl

          character(len=20) :: i_nodes_filename
          character(len=24) :: i_grdpts_id_filename
          character(len=20) :: i_sizes_filename

          character(len=11) :: bf_nodes_filename
          character(len=14) :: bf_grdpts_id_filename
          character(len=11) :: bf_sizes_filename
          character(len=6)  :: bf_nb_sbf_filename


          !determine the number of integer needed to write the
          !file index
          if(index.le.9) then
             format_index = 1
          else
             if((index.ge.10).and.(index.le.99)) then
                format_index = 2
             else
                print '(''test_bf_interface_prog'')'
                print '(''print_output'')'
                stop 'file_index not supported'
             end if
          end if


          !determine the format for the name of the output files
          write(i_format_nodes, '(''(A14,I'',I1,'',A4)'')') format_index
          write(i_format_grdpt, '(''(A18,I'',I1,'',A4)'')') format_index

          write(bf_format_nodes, '(''(A5,I'',I1,'',A4)'')') format_index
          write(bf_format_grdpt, '(''(A8,I'',I1,'',A4)'')') format_index
          write(bf_format_nbsbl, '(''(I'',I1,'',A4)'')'  )  format_index


          !determine the names of the output files
          write(i_nodes_filename, i_format_nodes)
     $         'interior_nodes', index, '.dat'
          write(i_grdpts_id_filename, i_format_grdpt)
     $         'interior_grdpts_id', index, '.dat'
          write(i_sizes_filename, i_format_nodes)
     $         'interior_sizes', index, '.dat'
                    
          write(bf_nodes_filename, bf_format_nodes)
     $         'nodes', index, '.dat'
          write(bf_grdpts_id_filename, bf_format_grdpt)
     $         'grdpt_id', index, '.dat'
          write(bf_sizes_filename, bf_format_nodes)
     $         'sizes', index, '.dat'
          write(bf_nb_sbf_filename, bf_format_nbsbl)
     $         index, '.dat'

          !write interior data          
          call print_interior_data(
     $         nodes,
     $         grdpts_id,
     $         i_nodes_filename,
     $         i_grdpts_id_filename,
     $         i_sizes_filename)
          
          !write data stored in the buffer layers
          call interface_used%print_binary(
     $         bf_nodes_filename,
     $         bf_grdpts_id_filename,
     $         bf_sizes_filename,
     $         bf_nb_sbf_filename)

          !write the position of the detectors on output file
          call interface_used%print_idetectors_on_binary(
     $         index)

        end subroutine print_state      

      end module test_case_class
