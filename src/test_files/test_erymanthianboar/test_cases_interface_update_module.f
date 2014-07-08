      module test_cases_interface_update_module

        use bf_interface_class    , only : bf_interface
        use bf_interface_icr_class, only : bf_interface_icr
        use bf_mainlayer_class    , only : bf_mainlayer
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_bf_layer   , only : no_pt
        use parameters_input      , only : nx,ny,ne,bc_size
        use parameters_kind       , only : ikind, rkind
        use test_bf_layer_module  , only : print_interior_data

        implicit none

        private
        public :: update_nodes,
     $            print_state


        contains

        !< print the content of the interior, the buffer layer and the
        !> increasing detectors on output files
        subroutine print_state(nodes, grdpts_id, interface_used, index)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface_icr)         , intent(in)    :: interface_used
          integer                         , intent(inout) :: index


          character(len=20) :: i_nodes_filename
          character(len=24) :: i_grdpts_id_filename
          character(len=20) :: i_sizes_filename

          character(len=10) :: bf_nodes_filename
          character(len=13) :: bf_grdpts_id_filename
          character(len=10) :: bf_sizes_filename
          character(len=5) :: bf_nb_sbf_filename


          write(i_nodes_filename,
     $         '(''interior_nodes'',I1,''.dat'')') index
          write(i_grdpts_id_filename,
     $         '(''interior_grdpts_id'',I1,''.dat'')') index
          write(i_sizes_filename,
     $         '(''interior_sizes'',I1,''.dat'')') index
                    
          write(bf_nodes_filename,
     $         '(''nodes'',I1,''.dat'')') index
          write(bf_grdpts_id_filename,
     $         '(''grdpt_id'',I1,''.dat'')') index
          write(bf_sizes_filename,
     $         '(''sizes'',I1,''.dat'')') index
          write(bf_nb_sbf_filename,
     $         '(I1,''.dat'')') index


          call print_interior_data(
     $         nodes,
     $         grdpts_id,
     $         i_nodes_filename,
     $         i_grdpts_id_filename,
     $         i_sizes_filename)
          
          call interface_used%print_binary(
     $         bf_nodes_filename,
     $         bf_grdpts_id_filename,
     $         bf_sizes_filename,
     $         bf_nb_sbf_filename)

          call interface_used%print_idetectors_on_binary(
     $         index)

        end subroutine print_state


        !< update the data in the interior nodes and the buffer layers
        !> to simulate a vapor bubble moving in the computational domain
        subroutine update_nodes(
     $       test_case,
     $       timestep, dx, dy, 
     $       interior_nodes, interface_used)


          implicit none

          integer                         , intent(in)    :: test_case
          integer                         , intent(in)    :: timestep
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne), intent(out)   :: interior_nodes
          class(bf_interface)             , intent(inout) :: interface_used


          !update the interior nodes
          call update_interior_nodes(test_case, timestep, dx, dy, interior_nodes)

          !update the data in the buffer layers
          call update_bf_nodes(test_case, timestep, dx, dy, interface_used)

        end subroutine update_nodes


        !< write data in the interior nodes to simulate a
        !> vapor bubble moving in the interior domain
        subroutine update_interior_nodes(test_case, timestep, dx, dy, interior_nodes)

          implicit none

          integer                         , intent(in)  :: test_case
          integer                         , intent(in)  :: timestep
          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes

          real(rkind), dimension(2) :: c_coords
          real(rkind)               :: d_liq, d_vap
          real(rkind)               :: l_interface
          real(rkind)               :: radius, x_center, y_center
          real(rkind)               :: mass
          real(rkind), dimension(2) :: velocity

          integer(ikind) :: i,j


          !parameters constraining the bubble
          call get_bubble_param(
     $         test_case,
     $         timestep, dx, dy,
     $         d_liq, d_vap, l_interface, radius,
     $         x_center, y_center)
          

          !initialization of the nodes coresponding
          !to the bubble
          do j=1, ny
             do i=1, nx

                c_coords = get_coords(i,j,dx,dy)

                mass     = get_mass(c_coords,
     $                              d_liq, d_vap,
     $                              l_interface,
     $                              radius,
     $                              x_center, y_center)

                velocity = get_velocity(test_case,
     $                                  c_coords,
     $                                  x_center, y_center)
                
                interior_nodes(i,j,1) = mass
                interior_nodes(i,j,2) = mass*velocity(1)
                interior_nodes(i,j,3) = mass*velocity(2)

             end do
          end do

        end subroutine update_interior_nodes


        !< write data in the buffer layer nodes to simulate a
        !> vapor bubble moving in the buffer layers
        subroutine update_bf_nodes(test_case, timestep, dx, dy, interface_used)

          implicit none

          integer            , intent(in)    :: test_case
          integer            , intent(in)    :: timestep
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
          real(rkind)   , dimension(2)                  :: c_coords

          real(rkind)               :: d_liq, d_vap
          real(rkind)               :: l_interface
          real(rkind)               :: radius, x_center, y_center
          real(rkind)               :: mass
          real(rkind), dimension(2) :: velocity
          

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
                   
                   !parameters constraining the bubble
                   call get_bubble_param(
     $                  test_case,
     $                  timestep,
     $                  dx, dy,
     $                  d_liq, d_vap, l_interface, radius,
     $                  x_center, y_center)

                   !fill the nodes
                   do j=1, sizes(2)
                      do i=1, sizes(1)
                         if(grdpts_id(i,j).ne.no_pt) then
                            g_coords(1) = alignment_tab(1,1)-bc_size+(i-1)
                            g_coords(2) = alignment_tab(2,1)-bc_size+(j-1)

                            c_coords = get_coords(g_coords(1),
     $                                            g_coords(2),
     $                                            dx,
     $                                            dy)

                            mass     = get_mass(c_coords,
     $                                          d_liq, d_vap,
     $                                          l_interface,
     $                                          radius,
     $                                          x_center, y_center)

                            velocity = get_velocity(test_case,
     $                                              c_coords,
     $                                              x_center, y_center)
                
                            new_nodes(i,j,1) = mass
                            new_nodes(i,j,2) = mass*velocity(1)
                            new_nodes(i,j,3) = mass*velocity(2)

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


        !< get the parameters contraining the bubble
        subroutine get_bubble_param(
     $     test_case,
     $     timestep,
     $     dx, dy,
     $     d_liq, d_vap, l_interface, radius,
     $     x_center, y_center)

          implicit none

          integer    , intent(in)  :: test_case
          integer    , intent(in)  :: timestep
          real(rkind), intent(in)  :: dx
          real(rkind), intent(in)  :: dy
          real(rkind), intent(out) :: d_liq
          real(rkind), intent(out) :: d_vap
          real(rkind), intent(out) :: l_interface
          real(rkind), intent(out) :: radius
          real(rkind), intent(out) :: x_center
          real(rkind), intent(out) :: y_center

          d_liq       = 1.1
          d_vap       = 0.1
          l_interface = 0.1
          radius      = 0.3

          select case(test_case)
            case(1)
               x_center    = 0.5+timestep*dx
               y_center    = 0.5
            case(2)
               x_center    = 0.5-timestep*dx
               y_center    = 0.5
            case(3)
               x_center    = 0.5
               y_center    = 0.5+timestep*dy
            case(4)
               x_center    = 0.5
               y_center    = 0.5-timestep*dy
            case(5)
               x_center    = 0.5+timestep*dx
               y_center    = 0.5+timestep*dy
            case(6)
               x_center    = 0.5+timestep*dx
               y_center    = 0.5-timestep*dy
            case(7)
               x_center    = 0.5-timestep*dx
               y_center    = 0.5+timestep*dy
            case(8)
               x_center    = 0.5-timestep*dx
               y_center    = 0.5-timestep*dy
            case default
               print '(''test_cases_interface_update_module'')'
               print '(''get_bubble_param'')'
               stop 'test case not recognized'
          end select

        end subroutine get_bubble_param



        !< from teh general coordinates, find the
        !> 2-D Cartesian coordinates
        function get_coords(i,j,dx,dy)

          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)   , intent(in) :: dx
          real(rkind)   , intent(in) :: dy
          real(rkind), dimension(2)  :: get_coords

          get_coords(1) = (i-bc_size+1/2)*dx
          get_coords(2) = (j-bc_size+1/2)*dy

        end function get_coords


        !< compute the mass
        function get_mass(
     $     g_coords,
     $     d_liq, d_vap,
     $     l_interface,
     $     radius,
     $     x_center, y_center)

          implicit none

          real(rkind), dimension(2), intent(in) :: g_coords
          real(rkind)              , intent(in) :: d_liq
          real(rkind)              , intent(in) :: d_vap        
          real(rkind)              , intent(in) :: l_interface  
          real(rkind)              , intent(in) :: radius       
          real(rkind)              , intent(in) :: x_center     
          real(rkind)              , intent(in) :: y_center
          real(rkind)                           :: get_mass

          real(rkind) :: r

          r = SQRT((g_coords(1)-x_center)**2+(g_coords(2)-y_center)**2)

          get_mass = (d_liq+d_vap)/2.0d0 +
     $               (d_liq-d_vap)/2.0d0*Tanh(-2*(r-radius)/l_interface)

        end function get_mass


        !< compute the velocity
        function get_velocity(test_case, g_coords, x_center, y_center)

          implicit none

          integer                  , intent(in) :: test_case
          real(rkind), dimension(2), intent(in) :: g_coords
          real(rkind)              , intent(in) :: x_center     
          real(rkind)              , intent(in) :: y_center
          real(rkind), dimension(2)             :: get_velocity

          real(rkind) :: norm
          real(rkind) :: x,y
          
          x = g_coords(1)-x_center
          y = g_coords(2)-y_center

          norm = 1.0d0
          
          !get_velocity(1) =-1.0 !x/SQRT(x**2+y**2)*norm
          !get_velocity(2) = 0.0 !y/SQRT(x**2+y**2)*norm

          select case(test_case)
            case(1)
               get_velocity(1) =  1.0
               get_velocity(2) =  0.0
            case(2)
               get_velocity(1) = -1.0
               get_velocity(2) =  0.0
            case(3)
               get_velocity(1) =  0.0
               get_velocity(2) =  1.0
            case(4)
               get_velocity(1) =  0.0
               get_velocity(2) = -1.0
            case(5)
               get_velocity(1) =  1.0
               get_velocity(2) =  1.0
            case(6)
               get_velocity(1) =  1.0
               get_velocity(2) = -1.0
            case(7)
               get_velocity(1) = -1.0
               get_velocity(2) =  1.0
            case(8)
               get_velocity(1) = -1.0
               get_velocity(2) = -1.0
            case default
               print '(''test_cases_interface_update_module'')'
               print '(''get_velocity'')'
               stop 'test case not recognized'
          end select          

        end function get_velocity


      end module test_cases_interface_update_module
