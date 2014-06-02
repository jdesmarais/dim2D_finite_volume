      program test_bf_interface_icr_prog

        use bf_interface_icr_class, only : bf_interface_icr
        use parameters_input      , only : nx,ny,ne
        use parameters_kind       , only : ikind, rkind
        use test_bf_layer_module  , only : print_interior_data,
     $                                     ini_grdpts_id

        implicit none


        type(bf_interface_icr)           :: interface_used
        real(rkind), dimension(nx,ny,ne) :: nodes
        integer    , dimension(nx,ny)    :: grdpts_id


        !initialization of nodes and grdpts_id
        call ini_cst_nodes(nodes)
        call ini_grdpts_id(grdpts_id)

        !test of ini()
        call interface_used%ini()
        call interface_used%print_idetectors_on(nodes(:,:,1))

        !print the detector positions
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes1.dat',
     $                           'interior_grdpts_id1.dat',
     $                           'interior_sizes1.dat')
        
        contains


        subroutine ini_cst_nodes(nodes)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          integer(ikind) :: i,j

          do j=1, ny
             do i=1, nx
                nodes(i,j,1)  = 1.0
             end do
          end do

        end subroutine ini_cst_nodes

      end program test_bf_interface_icr_prog
