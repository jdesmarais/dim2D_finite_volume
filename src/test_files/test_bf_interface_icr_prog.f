      program test_bf_interface_icr_prog

        use bf_interface_icr_class, only : bf_interface_icr
        use parameters_input      , only : nx,ny,ne
        use parameters_kind       , only : ikind, rkind
        use test_bf_layer_module  , only : print_interior_data,
     $                                     ini_grdpts_id

        implicit none


        type(bf_interface_icr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer(ikind), dimension(7,7,2)    :: cpt_coords_tab
        integer(ikind), dimension(2)        :: cpt_coords
        integer       , dimension(3,3)      :: nbc_template
        integer :: i,j


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

        !test of create_nbc_interior_pt_template()
        call ini_cpt_coords_table(cpt_coords_tab)
        do j=1,7
           do i=1,7
              cpt_coords = cpt_coords_tab(i,j,:)
              nbc_template = interface_used%create_nbc_interior_pt_template(
     $             cpt_coords)
              call print_nbc_template(i,j,nbc_template)
           end do
        end do

        
        contains

        !< print the content of the matrix nbc_template
        !> on an output file
        subroutine print_nbc_template(i,j,nbc_template)

          implicit none

          integer                , intent(in) :: i
          integer                , intent(in) :: j
          integer, dimension(3,3), intent(in) :: nbc_template


          character(len=20) :: nbc_filename
          integer :: ios


          write(nbc_filename,
     $         '(''nbc_template'',I1,I1,''.dat'')') i,j

          open(unit=1,
     $          file=nbc_filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nbc_template
              close(unit=1)
           else
              stop 'file opening pb'
           end if


        end subroutine print_nbc_template


        !< initialize the coordinates tested as cpt_coords
        subroutine ini_cpt_coords_table(cpt_coords)
        
          implicit none

          integer(ikind), dimension(7,7,2), intent(out) :: cpt_coords

          integer(ikind), dimension(7) :: x_coords
          integer(ikind), dimension(7) :: y_coords
          integer :: i,j

          x_coords = [ 1, 2, 3, 4, nx-2, nx-1, nx]
          y_coords = [ 1, 2, 3, 4, ny-2, ny-1, ny]

          do j=1,7
             do i=1,7
                cpt_coords(i,j,1) = x_coords(i)
                cpt_coords(i,j,2) = y_coords(j)
             end do
          end do

        end subroutine ini_cpt_coords_table


        !< initialize the nodes using a constant variable
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
