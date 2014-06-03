      program test_bf_interface_icr_update_prog

        use bf_detector_i_list_class, only : bf_detector_i_list
        use bf_interface_icr_class  , only : bf_interface_icr
c$$$        use bf_sublayer_class       , only : bf_sublayer
c$$$        use parameters_bf_layer     , only : align_N, interior_pt, no_pt
c$$$        use parameters_constant     , only : N,S,E,W
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind
        use test_bf_layer_module    , only : print_interior_data,
     $                                       ini_grdpts_id

        implicit none


        type(bf_interface_icr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
        integer       , dimension(nx,ny)    :: grdpts_id

        integer :: i

        real(rkind) :: dx
        real(rkind) :: dy

        integer :: file_index

        dx = 1.0d0/(nx-2*bc_size)
        dy = 1.0d0/(ny-2*bc_size)


        file_index = 1

        !initialization of th d
        call ini_nodes(dx,dy,interior_nodes,file_index-1)
        call ini_grdpts_id(grdpts_id)

        !initialization of the interface
        call interface_used%ini()

        !printing the initial state
        call print_state(
     $       interior_nodes, grdpts_id,
     $       interface_used,
     $       file_index)

        !update the buffer layers
        do i=1,5
           call ini_nodes(dx,dy,interior_nodes,file_index-1)

           call interface_used%update_bf_layers_with_idetectors(
     $          interior_nodes, dx, dy)

           call print_state(
     $       interior_nodes, grdpts_id,
     $       interface_used,
     $       file_index)
        end do

        

        contains


        subroutine ini_nodes(dx,dy,interior_nodes,index)

          implicit none

          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer                         , intent(in)  :: index

          real(rkind), dimension(2) :: g_coords
          real(rkind)               :: d_liq, d_vap
          real(rkind)               :: l_interface
          real(rkind)               :: radius, x_center, y_center
          real(rkind)               :: mass
          real(rkind), dimension(2) :: velocity

          integer(ikind) :: i,j

          d_liq       = 1.1
          d_vap       = 0.1
          l_interface = 0.1
          radius      = 0.3
          x_center    = 0.5+index*dx
          y_center    = 0.5+index*dy

          do j=1, ny
             do i=1, nx

                g_coords = get_coords(i,j,dx,dy)

                mass     = get_mass(g_coords,
     $                              d_liq, d_vap,
     $                              l_interface,
     $                              radius,
     $                              x_center, y_center)

                velocity = [1.0, 0.0] !get_velocity(g_coords,
!     $                                  x_center, y_center)
                
                interior_nodes(i,j,1) = mass
                interior_nodes(i,j,2) = mass*velocity(1)
                interior_nodes(i,j,3) = mass*velocity(2)

             end do
          end do

        end subroutine ini_nodes


        function get_coords(i,j,dx,dy)

          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)   , intent(in) :: dx
          real(rkind)   , intent(in) :: dy
          real(rkind), dimension(2)  :: get_coords

          get_coords(1) = (i-bc_size+1/2)*dx
          get_coords(2) = (j-bc_size+1/2)*dy

        end function get_coords


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


        function get_velocity(g_coords, x_center, y_center)

          implicit none

          real(rkind), dimension(2), intent(in) :: g_coords
          real(rkind)              , intent(in) :: x_center     
          real(rkind)              , intent(in) :: y_center
          real(rkind), dimension(2)             :: get_velocity

          real(rkind) :: norm
          real(rkind) :: x,y
          
          x = g_coords(1)-x_center
          y = g_coords(2)-y_center

          norm = 1.0d0

          get_velocity(1) = x/SQRT(x**2+y**2)*norm
          get_velocity(2) = y/SQRT(x**2+y**2)*norm

        end function get_velocity


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

          index = index+1

        end subroutine print_state


      end program test_bf_interface_icr_update_prog
