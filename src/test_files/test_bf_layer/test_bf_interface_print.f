      program test_bf_interface_print

        use bf_interface_print_class, only :
     $     bf_interface_print

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       interior_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        call test_ini()
        print '(''test_ini: check interior_grdpts_id.nc'')'
        print '()'


        call test_print_netcdf()
        print '(''test_print_netcdf: check N_1_10.nc...'')'
        print '()'

        
        print '(''test_validated: '',L1)', .true.

        contains


        subroutine test_ini()

          implicit none

          type(bf_interface_print)   :: bf_interface_used
          real(rkind), dimension(nx) :: interior_x_map
          real(rkind), dimension(ny) :: interior_y_map

          call bf_interface_used%ini(interior_x_map,interior_y_map)

        end subroutine test_ini


        subroutine test_print_netcdf()

          implicit none

          type(bf_interface_print)         :: bf_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: added_sublayer

          integer(ikind) :: i,j
          integer        :: k


          !input
          interior_nodes = reshape((/
     $         (((200*(k-1) + 20*(j-1) + (i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          call bf_interface_used%ini(interior_x_map,interior_y_map)
          

          !two North layers
          call bf_interface_used%mainlayer_pointers(N)%ini_mainlayer(N)
          added_sublayer => bf_interface_used%mainlayer_pointers(N)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_N, align_W+3, align_N+1/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_N-3+j-1) + (align_W-7+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          added_sublayer => bf_interface_used%mainlayer_pointers(N)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E-3, align_N, align_E+4, align_N+1/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_N-3+j-1) + (align_E-6+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))


          !two south layers
          call bf_interface_used%mainlayer_pointers(S)%ini_mainlayer(S)
          added_sublayer => bf_interface_used%mainlayer_pointers(S)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_S-1, align_W+3, align_S/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_S-4+j-1) + (align_W-7+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          added_sublayer => bf_interface_used%mainlayer_pointers(S)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E-3, align_S-1, align_E+4, align_S/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_S-4+j-1) + (align_E-6+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))


          !two east layers
          call bf_interface_used%mainlayer_pointers(E)%ini_mainlayer(E)
          added_sublayer => bf_interface_used%mainlayer_pointers(E)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E, align_S+1, align_E+4, align_S+2/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_S-2+j-1) + (align_E-3+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          added_sublayer => bf_interface_used%mainlayer_pointers(E)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E, align_N-2, align_E+4, align_N-1/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_N-5+j-1) + (align_E-3+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))


          !two west layers
          call bf_interface_used%mainlayer_pointers(W)%ini_mainlayer(W)
          added_sublayer => bf_interface_used%mainlayer_pointers(W)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_S+1, align_W, align_S+2/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_S-2+j-1) + (align_W-7+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          added_sublayer => bf_interface_used%mainlayer_pointers(W)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_N-2, align_W, align_N-1/),
     $            (/2,2/)))
          added_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          added_sublayer%nodes = reshape((/
     $         (((200*(k-1) + 20*(align_N-5+j-1) + (align_W-7+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))


          call bf_interface_used%print_netcdf(
     $         10,
     $         ['md','qx','qy','Ed'],
     $         ['mass density','momentum-x density','momemtum-y density','total energy density'],
     $         ['(kg.m-3)/(kg.m-3)','(kg.m-2.s-1)/(kg.m-2.s-1)','(kg.m-2.s-1)/(kg.m-2.s-1)','(J.m-3)/(J.m-3)'],
     $         0.01d0)


        end subroutine test_print_netcdf

      end program test_bf_interface_print
