      program test_bf_layer_remove_prog

        use bf_layer_remove_module, only : get_check_line_param
        use bf_interface_class    , only : bf_interface
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_bf_layer   , only : align_N, align_S,
     $                                     align_E, align_W,
     $                                     interior_pt
        use parameters_input      , only : nx,ny,ne,bc_size
        use parameters_constant   , only : N,S,E,W
        use parameters_kind       , only : ikind, rkind
        use test_bf_layer_module  , only : print_interior_data,
     $                                     ini_grdpts_id!,
!     $                                     ini_nodes


        implicit none

        
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        type(bf_interface)                  :: interface_used

        integer                                    :: k
        integer                                    :: mainlayer_id
        integer(ikind), dimension(2,2)             :: alignment
        type(bf_sublayer), pointer                 :: added_sublayer
        integer(ikind), dimension(2)               :: match_table
        integer(ikind), dimension(2,2)             :: bf_coords
        integer(ikind), dimension(2,2)             :: in_coords
        integer(ikind), dimension(2)               :: sizes
        real(rkind), dimension(:,:,:), allocatable :: bf_nodes
        real(rkind)                                :: color
        integer(ikind)                             :: i,j
        real(rkind) :: dx,dy


        !initialize the nodes
        call ini_cst_nodes(nodes)
        !call ini_nodes(nodes)

        !initialization of the grdpts_id
        call ini_grdpts_id(grdpts_id)

        !initialize the interface
        call interface_used%ini()

        !add the buffer layers and test the indices for
        !checking the grid points
        do k=1,8

           !determine the alignment + mainlayer_id
           call get_alignment(k, mainlayer_id, alignment)

           !add a new buffer layer
           added_sublayer => interface_used%allocate_sublayer(
     $          mainlayer_id,
     $          nodes,
     $          alignment,
     $          dx,
     $          dy)

           !match_table
           match_table = added_sublayer%get_general_to_local_coord_tab()
           
           !get buffer layer extent
           sizes = added_sublayer%get_sizes()
           
           !determine the indices for the check of the removal frontier
           call get_check_line_param(
     $          mainlayer_id, alignment, match_table,
     $          sizes(1), sizes(2),
     $          bf_coords, in_coords)
           
           !allocate and initialize nodes for the buffer layer           
           allocate(bf_nodes(sizes(1), sizes(2), ne))
           call ini_cst_nodes(bf_nodes)
           
           !determine the color
           color = 0.1 + (k-1)*0.1
           
           !set the color in the interior nodes
           do j=in_coords(2,1), in_coords(2,2)
              do i=in_coords(1,1), in_coords(1,2)
                 nodes(i,j,1) = color
              end do
           end do
           
           !set the color in the buffer layer nodes
           do j=bf_coords(2,1), bf_coords(2,2)
              do i=bf_coords(1,1), bf_coords(1,2)
                 bf_nodes(i,j,1) = color
                 call added_sublayer%set_grdpts_id_pt(i,j,interior_pt)
              end do
           end do
           
           !set the nodes to the buffer layer
           call added_sublayer%set_nodes(bf_nodes)

        end do

        !print interior data
        call print_interior_data(
     $       nodes,
     $       grdpts_id,
     $       'interior_nodes1.dat',
     $       'interior_grdpts_id1.dat',
     $       'interior_sizes1.dat')
        
        !print the interface
        call interface_used%print_binary(
     $       'nodes1.dat',
     $       'grdpt_id1.dat',
     $       'sizes1.dat',
     $       '1.dat')


        contains

        !< initialize the nodes using a constant variable
        subroutine ini_cst_nodes(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          integer(ikind) :: i,j

          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                nodes(i,j,1)  = 1.0
             end do
          end do

        end subroutine ini_cst_nodes


        !< alignments for the buffer layers initialized
        subroutine get_alignment(index, mainlayer, alignment)

          implicit none
 
          integer                       , intent(in)  :: index
          integer                       , intent(out) :: mainlayer
          integer(ikind), dimension(2,2), intent(out) :: alignment


          integer :: size_large_layer
          integer :: size_small_layer

          size_large_layer = 10
          size_small_layer = 5


          select case(index)
            case(1)
               mainlayer = N
               alignment(1,1) = - size_large_layer
               alignment(1,2) = - size_large_layer+size_small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N + size_small_layer

            case(2)
               mainlayer = N
               alignment(1,1) = 1
               alignment(1,2) = size_large_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N + size_small_layer+1

            case(3)
               mainlayer = N
               alignment(1,1) = nx-size_small_layer
               alignment(1,2) = alignment(1,1) + size_large_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N + size_small_layer+3

            case(4)
               mainlayer = E
               alignment(1,1) = align_E
               alignment(1,2) = align_E+size_small_layer
               alignment(2,1) = bc_size + size_small_layer
               alignment(2,2) = ny-bc_size - size_small_layer

            case(5)
               mainlayer = W
               alignment(1,1) = align_W-size_small_layer
               alignment(1,2) = align_W
               alignment(2,1) = ny-bc_size-size_large_layer
               alignment(2,2) = ny-bc_size

            case(6)
               mainlayer = S
               alignment(1,1) = - size_large_layer
               alignment(1,2) = - size_large_layer+size_small_layer
               alignment(2,1) = align_S - size_small_layer
               alignment(2,2) = align_S

            case(7)
               mainlayer = S
               alignment(1,1) = 1
               alignment(1,2) = size_large_layer
               alignment(2,1) = align_S - (size_small_layer+1)
               alignment(2,2) = align_S

            case(8)
               mainlayer = S
               alignment(1,1) = nx-size_small_layer
               alignment(1,2) = alignment(1,1) + size_large_layer
               alignment(2,1) = align_S - (size_small_layer+3)
               alignment(2,2) = align_S

          end select

        end subroutine get_alignment

      end program test_bf_layer_remove_prog
