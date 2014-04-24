      module interface_print_module

        use bf_mainlayer_class      , only : bf_mainlayer
        use bf_sublayer_class       , only : bf_sublayer
        use interface_abstract_class, only : interface_abstract

        implicit none

        private
        public :: print_interface


        contains


        subroutine print_interface(interface_used, interface_id)

          
          implicit none
          
          class(interface_abstract), intent(in) :: interface_used
          integer                  , intent(in) :: interface_id

          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr
          
          integer :: i,j
          integer :: nb_sublayers_max

          character(len=18) :: sublayers_nb_filename

          nb_sublayers_max = 1

          !go through the different main layers of the interface
          do i=1, size(interface_used%mainlayer_pointers,1)

             
             !get the main layer if it exists
             if(associated(interface_used%mainlayer_pointers(i)%ptr)) then
                
                mainlayer_ptr => interface_used%mainlayer_pointers(i)%ptr

                !get the sublayers inside the mainlayer
                if(associated(mainlayer_ptr%head_sublayer)) then

                   j=1
                   sublayer_ptr => mainlayer_ptr%head_sublayer                
                   call print_sublayer(sublayer_ptr,i,j,interface_id)

                   do while(associated(sublayer_ptr%next))
                      sublayer_ptr => sublayer_ptr%next
                      j = j+1
                      call print_sublayer(sublayer_ptr,i,j,interface_id)
                   end do
                   nb_sublayers_max = max(nb_sublayers_max,j)

                end if
             end if

          end do

          write(sublayers_nb_filename,'(''sublayers_nb_'',I1,''.dat'')')
     $         interface_id

          call print_nb_sublayers(
     $         sublayers_nb_filename,
     $         nb_sublayers_max)

        end subroutine print_interface


        subroutine print_sublayer(sublayer_ptr,i,j,interface_id)

          implicit none

          type(bf_sublayer), pointer, intent(in) :: sublayer_ptr
          integer                   , intent(in) :: i
          integer                   , intent(in) :: j
          integer                   , intent(in) :: interface_id

          character(2), dimension(8) :: bf_layer_char

          character(len=22) :: sizes_filename
          character(len=22) :: nodes_filename
          character(len=22) :: grdid_filename

          bf_layer_char = ['N_','S_','E_','W_','NE','NW','SE','SW']

          write(sizes_filename,'(A2,I1,''_sizes_'',I1,''.dat'')')
     $         bf_layer_char(i), j, interface_id
          write(nodes_filename,'(A2,I1,''_nodes_'',I1,''.dat'')')
     $         bf_layer_char(i), j, interface_id
          write(grdid_filename,'(A2,I1,''_grdpt_id_'',I1,''.dat'')')
     $         bf_layer_char(i), j, interface_id

          call sublayer_ptr%element%print_sizes(sizes_filename)
          call sublayer_ptr%element%print_nodes(nodes_filename)
          call sublayer_ptr%element%print_grdpts_id(grdid_filename)

        end subroutine print_sublayer


        subroutine print_nb_sublayers(filename,nb_sublayers)

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

        end subroutine print_nb_sublayers

      end module interface_print_module
