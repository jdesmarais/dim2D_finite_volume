      program test_bf_corner_prog

        use bf_corner_module    , only : is_alignment_compatible_with_corner

        use parameters_constant , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input    , only : nx,ny,ne,bc_size
        use parameters_kind     , only : ikind, rkind

        use test_bf_layer_module, only : print_interior_data,
     $                                   ini_nodes,
     $                                   ini_grdpts_id

        implicit none

        real(rkind), dimension(nx,ny,ne)     :: nodes
        integer    , dimension(nx,ny)        :: grdpts_id
        integer, dimension(2,2)              :: alignment
        logical, dimension(4)                :: neighbors
        logical                              :: compatible
        integer, dimension(2,2)              :: new_alignment
        logical, dimension(4)                :: new_neighbors

        integer :: mainlayer_id
        integer :: corner_id
        integer(ikind) :: i,j
        integer :: k


        !initialize the nodes
        do k=1, size(nodes,3)
           do j=1, size(nodes,2)
              do i=1, size(nodes,1)
                 nodes(i,j,k) = 1.0
              end do
           end do
        end do
        call ini_grdpts_id(grdpts_id)      
        

        !test N & N_E
        call test_y_direction(
     $       N, N_E,
     $       ny-1, nodes)

        !test N & N_W
        call test_y_direction(
     $       N, N_W,
     $       ny, nodes)

        !test S & S_E
        call test_y_direction(
     $       S, S_E,
     $       1, nodes)

        !test S & S_W
        call test_y_direction(
     $       S, S_W,
     $       2, nodes)

        !test E & N_E
        call test_x_direction(
     $       E, N_E,
     $       nx-1, nodes)

        !test E & S_E
        call test_x_direction(
     $       E, S_E,
     $       nx, nodes)

        !test W & N_W
        call test_x_direction(
     $       W, N_W,
     $       1, nodes)

        !test S & S_W
        call test_x_direction(
     $       W, S_W,
     $       2, nodes)

        !print the interior data
        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes.dat',
     $                           'interior_grdpts_id.dat',
     $                           'interior_sizes.dat')


        contains

        subroutine test_y_direction(
     $       mainlayer_id, corner_id,
     $       j, nodes)
        
          implicit none

          integer                         , intent(in)    :: mainlayer_id
          integer                         , intent(in)    :: corner_id
          integer(ikind)                  , intent(in)    :: j
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer(ikind)          :: i
          integer, dimension(2,2) :: alignment
          logical, dimension(4)   :: neighbors
          logical                 :: compatible
          
          do i=bc_size+1, nx-bc_size
             
             alignment(1,1) = i
             alignment(1,2) = i
             alignment(2,1) = j
             alignment(2,2) = j
             
             neighbors = [.false.,.false.,.false.,.false.]
             
             compatible = is_alignment_compatible_with_corner(
     $            corner_id, mainlayer_id,
     $            alignment, neighbors,
     $            new_alignment, new_neighbors)
             
             if(compatible) then
                nodes(i,j,1) = 0.5
             else
                nodes(i,j,1) = 0.9
             end if
           
          end do          

        end subroutine test_y_direction


        subroutine test_x_direction(
     $       mainlayer_id, corner_id,
     $       i, nodes)
        
          implicit none

          integer                         , intent(in)    :: mainlayer_id
          integer                         , intent(in)    :: corner_id
          integer(ikind)                  , intent(in)    :: i
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer(ikind)          :: j
          integer, dimension(2,2) :: alignment
          logical, dimension(4)   :: neighbors
          logical                 :: compatible
          
          do j=bc_size+1, ny-bc_size
             
             alignment(1,1) = i
             alignment(1,2) = i
             alignment(2,1) = j
             alignment(2,2) = j
             
             neighbors = [.false.,.false.,.false.,.false.]
             
             compatible = is_alignment_compatible_with_corner(
     $            corner_id, mainlayer_id,
     $            alignment, neighbors,
     $            new_alignment, new_neighbors)
             
             if(compatible) then
                nodes(i,j,1) = 0.5
             else
                nodes(i,j,1) = 0.9
             end if
           
          end do          

        end subroutine test_x_direction

      end program test_bf_corner_prog
