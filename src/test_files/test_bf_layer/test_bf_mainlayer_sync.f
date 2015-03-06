      program test_bf_mainlayer_sync

        use bf_mainlayer_sync_class, only :
     $     bf_mainlayer_sync

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_matrix3D_validated,
     $       is_int_vector_validated

        use parameters_bf_layer, only :
     $       align_N,align_S,
     $       align_E,align_W

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_update_interior_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_interior_bc_sections: '',L1)', test_loc
        print '()'


        test_loc = test_sync_nodes_with_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes_with_interior: '',L1)', test_loc
        print '()'


        contains


        function test_update_interior_bc_sections(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc

          
          test_validated = .true.


          do k=1,14

             test_loc = perform_test_update_interior_bc_sections(
     $            k,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_update_interior_bc_sections


        function perform_test_update_interior_bc_sections(
     $       test_id,
     $       detailled)
     $       result(test_validated)

           implicit none

           integer, intent(in) :: test_id
           logical, intent(in) :: detailled
           logical             :: test_validated

           type(bf_mainlayer_sync) :: bf_mainlayer_used

           integer                   :: mainlayer_id
           integer                   :: nb_sublayers
           integer, dimension(2,2,3) :: bf_alignment
           integer                   :: nb_sections_test
           integer, dimension(2,2)   :: bc_sections_test

           type(bf_sublayer), pointer       :: added_sublayer
           real(rkind), dimension(nx)       :: interior_x_map
           real(rkind), dimension(ny)       :: interior_y_map
           real(rkind), dimension(nx,ny,ne) :: interior_nodes

           integer, dimension(:,:), allocatable :: bc_sections


           integer :: k
           integer :: test_loc


           test_validated = .true.


           !input
           call get_param_test_update_interior_bc_sections(
     $          test_id,
     $          mainlayer_id,
     $          nb_sublayers,
     $          bf_alignment,
     $          nb_sections_test,
     $          bc_sections_test)

           call bf_mainlayer_used%ini(mainlayer_id)

           do k=1, nb_sublayers
              added_sublayer => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,k))
           end do


           !output
           call bf_mainlayer_used%update_interior_bc_sections(bc_sections)
           
           !validation
           if(allocated(bc_sections)) then
              
              test_loc = nb_sections_test.eq.size(bc_sections,2)
              test_validated = test_validated.and.test_loc
              if(detailled.and.(.not.test_loc)) then
                 print '(''nb_sections failed'')'
              end if

              do k=1, nb_sections_test

                 test_loc = is_int_vector_validated(
     $                bc_sections(:,k),
     $                bc_sections_test(:,k),
     $                detailled)
                 test_validated = test_validated.and.test_loc

              end do

           else

              test_loc = nb_sections_test.eq.0
              test_validated = test_validated.and.test_loc
              if(detailled.and.(.not.test_loc)) then
                 print '(''nb_sections failed'')'
              end if

           end if

        end function perform_test_update_interior_bc_sections


        subroutine get_param_test_update_interior_bc_sections(
     $     k,
     $     mainlayer_id,
     $     nb_sublayers,
     $     bf_alignment,
     $     nb_sections,
     $     bc_sections)

           implicit none

           integer                  , intent(in)  :: k
           integer                  , intent(out) :: mainlayer_id
           integer                  , intent(out) :: nb_sublayers
           integer, dimension(2,2,3), intent(out) :: bf_alignment
           integer                  , intent(out) :: nb_sections
           integer, dimension(2,2)  , intent(out) :: bc_sections


           select case(k)

             case(1,7)
                nb_sublayers = 0
                if(k.eq.1) then
                   mainlayer_id = N
                else
                   mainlayer_id = S
                end if
                nb_sections = 1
                bc_sections(:,1) = [1,nx]
                
             case(2,8)
                nb_sublayers = 1

                if(k.eq.2) then
                   mainlayer_id = N
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+1, align_N  , align_W+3, align_N+1/),
     $                  (/2,2/))

                else
                   mainlayer_id = S
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+1, align_S-1, align_W+3, align_S/),
     $                  (/2,2/))

                end if
                nb_sections = 1
                bc_sections(:,1) = [align_W+6,nx]

             case(3,9)
                nb_sublayers = 1

                if(k.eq.3) then
                   mainlayer_id = N
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+10,   align_N, align_E-10, align_N+1/),
     $                  (/2,2/))
                else
                   mainlayer_id = S
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+10, align_S-1, align_E-10, align_S/),
     $                  (/2,2/))
                end if
                nb_sections = 2
                bc_sections(:,1) = [1,align_W+7]
                bc_sections(:,2) = [align_E-7,nx]

             case(4,10)
                nb_sublayers = 1

                if(k.eq.4) then
                   mainlayer_id = N
                   bf_alignment(:,:,1) = reshape((/
     $                  align_E-4,   align_N, align_E, align_N+1/),
     $                  (/2,2/))
                else
                   mainlayer_id = S
                   bf_alignment(:,:,1) = reshape((/
     $                  align_E-4, align_S-1, align_E, align_S/),
     $                  (/2,2/))
                end if
                nb_sections = 1
                bc_sections(:,1) = [1,align_E-7]

             case(5,11)
                nb_sublayers = 3

                if(k.eq.5) then
                   mainlayer_id = N
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W-2, align_N, align_W+2, align_N+1/),
     $                  (/2,2/))

                   bf_alignment(:,:,2) = reshape((/
     $                  align_W+9, align_N, align_E-9, align_N+1/),
     $                  (/2,2/))

                   bf_alignment(:,:,3) = reshape((/
     $                  align_E-2, align_N, align_E+2, align_N+1/),
     $                  (/2,2/))

                else
                   mainlayer_id = S
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W-2, align_S-1, align_W+2, align_S/),
     $                  (/2,2/))

                   bf_alignment(:,:,2) = reshape((/
     $                  align_W+9, align_S-1, align_E-9, align_S/),
     $                  (/2,2/))

                   bf_alignment(:,:,3) = reshape((/
     $                  align_E-2, align_S-1, align_E+2, align_S/),
     $                  (/2,2/))
                end if
                nb_sections = 2
                bc_sections(:,1) = [align_W+5,align_W+6]
                bc_sections(:,2) = [align_E-6,align_E-5]

             case(6,12)
                nb_sublayers = 1

                if(k.eq.6) then
                   mainlayer_id = N
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+1, align_N, align_E-1, align_N+1/),
     $                  (/2,2/))

                else
                   mainlayer_id = S
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W+1, align_S-1, align_E-1, align_S/),
     $                  (/2,2/))
                end if
                nb_sections = 0

             case(13,15)
                nb_sublayers = 1
                
                if(k.eq.13) then
                   mainlayer_id = W
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W-3, align_S+5, align_W, align_N-5/),
     $                  (/2,2/))
                else
                   mainlayer_id = E
                   bf_alignment(:,:,1) = reshape((/
     $                  align_E, align_S+5, align_E+3, align_N-5/),
     $                  (/2,2/))
                end if
                nb_sections = 2
                bc_sections(:,1) = [align_S+1,align_S+2]
                bc_sections(:,2) = [align_N-2,align_N-1]

             case(14,16)
                nb_sublayers = 1
                
                if(k.eq.14) then
                   mainlayer_id = W
                   bf_alignment(:,:,1) = reshape((/
     $                  align_W-3, align_S+1, align_W, align_N-1/),
     $                  (/2,2/))
                else
                   mainlayer_id = E
                   bf_alignment(:,:,1) = reshape((/
     $                  align_E, align_S+1, align_E+3, align_N-1/),
     $                  (/2,2/))
                end if
                nb_sections = 0

             end select

        end subroutine get_param_test_update_interior_bc_sections


        function test_sync_nodes_with_interior(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_mainlayer_sync)          :: bf_mainlayer_used
          integer, dimension(2,2,3)        :: bf_alignment
          type(bf_sublayer), pointer       :: added_sublayer1
          type(bf_sublayer), pointer       :: added_sublayer2
          type(bf_sublayer), pointer       :: added_sublayer3
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map

          real(rkind), dimension(9    ,6,ne) :: bf_nodes1
          real(rkind), dimension(nx-16,6,ne) :: bf_nodes2
          real(rkind), dimension(9    ,6,ne) :: bf_nodes3
          real(rkind), dimension(nx,ny,ne)   :: interior_nodes

          real(rkind), dimension(9    ,6,ne) :: bf_nodes1_test
          real(rkind), dimension(nx-16,6,ne) :: bf_nodes2_test
          real(rkind), dimension(9    ,6,ne) :: bf_nodes3_test
          real(rkind), dimension(nx,ny,ne)   :: interior_nodes_test

          integer :: i,j,k
          integer :: test_loc


          test_validated = .true.

          
          !input
          bf_nodes1_test = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_W-5+i-1),
     $         i=1,size(bf_nodes1_test,1)),j=1,size(bf_nodes1_test,2)),k=1,ne)/),
     $         (/size(bf_nodes1_test,1),size(bf_nodes1_test,2),ne/))

          bf_nodes2_test = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_W+6+i-1),
     $         i=1,size(bf_nodes2_test,1)),j=1,size(bf_nodes2_test,2)),k=1,ne)/),
     $         (/size(bf_nodes2_test,1),size(bf_nodes2_test,2),ne/))

          bf_nodes3_test = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_E-5+i-1),
     $         i=1,size(bf_nodes3_test,1)),j=1,size(bf_nodes3_test,2)),k=1,ne)/),
     $         (/size(bf_nodes3_test,1),size(bf_nodes3_test,2),ne/))
          
          interior_nodes_test = reshape((/
     $         (((500*(k-1) + 50*(j-1) + 5*(i-1),
     $         i=1,size(interior_nodes_test,1)),j=1,size(interior_nodes_test,2)),k=1,ne)/),
     $         (/size(interior_nodes_test,1),size(interior_nodes_test,2),ne/))
          interior_nodes_test(7:8,ny-1:ny,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,2),k=1,ne)/),
     $         (/2,2,ne/))
          interior_nodes_test(13:14,ny-1:ny,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,2),k=1,ne)/),
     $         (/2,2,ne/))


          bf_nodes1 = bf_nodes1_test
          bf_nodes1(4:size(bf_nodes1,1),1:2,:) = reshape((/
     $         (((-99.0d0,i=1,size(bf_nodes1,1)),j=1,2),k=1,ne)/),
     $         (/size(bf_nodes1,1)-3,2,ne/))

          bf_nodes2 = bf_nodes2_test
          bf_nodes2(1:size(bf_nodes2,1),1:2,:) = reshape((/
     $         (((-99.0d0,i=1,size(bf_nodes2,1)),j=1,2),k=1,ne)/),
     $         (/size(bf_nodes2,1),2,ne/))

          bf_nodes3 = bf_nodes3_test
          bf_nodes3(1:size(bf_nodes3,1)-3,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,size(bf_nodes3,1)-3),j=1,2),k=1,ne)/),
     $         (/size(bf_nodes3,1)-3,2,ne/))

          interior_nodes = interior_nodes_test
          interior_nodes(:,ny-1:ny,:) = reshape((/
     $         (((-99.0d0,
     $         i=1,size(interior_nodes,1)),j=1,2),k=1,ne)/),
     $         (/size(interior_nodes,1),2,ne/))


          !input
          call bf_mainlayer_used%ini(N)

          bf_alignment(:,:,1) = reshape((/
     $         align_W-2, align_N, align_W+2, align_N+1/),
     $         (/2,2/))
          
          bf_alignment(:,:,2) = reshape((/
     $         align_W+9, align_N, align_E-9, align_N+1/),
     $         (/2,2/))

          bf_alignment(:,:,3) = reshape((/
     $         align_E-2, align_N, align_E+2, align_N+1/),
     $         (/2,2/))

          added_sublayer1 => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,1))
          added_sublayer1%nodes = bf_nodes1

          added_sublayer2 => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,2))
          added_sublayer2%nodes = bf_nodes2

          added_sublayer3 => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,3))
          added_sublayer3%nodes = bf_nodes3


          !output
          call bf_mainlayer_used%sync_nodes_with_interior(
     $         interior_nodes)


          !validation
          test_loc = is_real_matrix3D_validated(
     $         added_sublayer1%nodes,
     $         bf_nodes1_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''sync bf_sublayer1 failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         added_sublayer2%nodes,
     $         bf_nodes2_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''sync bf_sublayer2 failed'')'
          end if
          
          test_loc = is_real_matrix3D_validated(
     $         added_sublayer3%nodes,
     $         bf_nodes3_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''sync bf_sublayer3 failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         interior_nodes,
     $         interior_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''sync interior_nodes failed'')'
          end if


        end function test_sync_nodes_with_interior


        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (nx.eq.20).and.(ny.eq.25))) then

             print '(''the test requires: '')'
             print '(''nx=20'')'
             print '(''ny=25'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_mainlayer_sync
