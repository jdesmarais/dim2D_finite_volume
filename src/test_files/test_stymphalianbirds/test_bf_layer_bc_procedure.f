      program test_bf_layer_bc_procedure

        use bf_layer_bc_procedure_module, only :
     $     SW_corner_type,
     $     SE_corner_type,
     $     NW_corner_type,
     $     NE_corner_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     N_edge_type,
     $     SE_edge_type,
     $     SW_edge_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     get_bc_interior_pt_procedure,
     $     get_bc_pt_procedure,
     $     is_bc_pt_procedure_of_edge_type,
     $     does_corner_procedure_1_match,
     $     does_corner_procedure_2_match,
     $     does_corner_procedure_3_match,
     $     does_corner_procedure_4_match,
     $     does_corner_procedure_5_1_match,
     $     does_corner_procedure_5_2_match,
     $     does_corner_procedure_6_1_match,
     $     does_corner_procedure_6_2_match

        use parameters_bf_layer, only :
     $     interior_pt,
     $     bc_interior_pt,
     $     bc_pt,
     $     no_pt

        implicit none


        integer                              :: test_id        
        integer, dimension(:,:), allocatable :: grdpts_id
        integer                              :: procedure_type
        integer                              :: test_procedure_type

        integer                              :: test_i
        integer                              :: test_j
        logical                              :: edge_type
        logical                              :: test_edge_type

        logical                              :: match
        integer                              :: nb_pts_x
        integer                              :: nb_pts_y
        integer                              :: test_nb_pts_x
        integer                              :: test_nb_pts_y
        logical                              :: test_match
                                             
        logical                              :: test_validated



        !test get_bc_interior_pt_procedure()
        print '(''test_get_bc_interior_pt_procedure'')'
        allocate(grdpts_id(3,3))
        do test_id=1,16

           call make_test_bc_interior_pt_procedure(
     $          test_id,
     $          grdpts_id,
     $          test_procedure_type)

           procedure_type = get_bc_interior_pt_procedure(
     $          2,
     $          2,
     $          grdpts_id)

           test_validated = procedure_type.eq.test_procedure_type

           print '(''test '',I2,'': '',L1)', test_id, test_validated
           
        end do
        deallocate(grdpts_id)
        print '()'


        !test is_bc_pt_procedure_of_edge_type
        print '(''test_is_bc_pt_procedure_of_edge_type'')'
        do test_id=1,12

           call make_test_is_bc_pt_procedure_of_edge_type(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_procedure_type,
     $          test_edge_type)

           edge_type = is_bc_pt_procedure_of_edge_type(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type)

           if(test_edge_type) then
              test_validated = edge_type.eqv.test_edge_type
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)

           else
              test_validated = edge_type.eqv.test_edge_type
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  edge_type: '',L1,2X,L1)', test_edge_type, edge_type
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
           end if              

           deallocate(grdpts_id)

        end do
        print '()'


        !test does_corner_procedure_1_match
        print '(''test_does_corner_procedure_1_match'')'
        do test_id=1,8
           
           call make_test_does_corner_procedure_1_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_1_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'        


        !test does_corner_procedure_2_match
        print '(''test_does_corner_procedure_2_match'')'
        do test_id=1,8
           
           call make_test_does_corner_procedure_2_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_2_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test does_corner_procedure_3_match
        print '(''test_does_corner_procedure_3_match'')'
        do test_id=1,6
           
           call make_test_does_corner_procedure_3_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_3_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test does_corner_procedure_4_match
        print '(''test_does_corner_procedure_4_match'')'
        do test_id=1,6
           
           call make_test_does_corner_procedure_4_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_4_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'

        
        !test does_corner_procedure_5_1_match
        print '(''test_does_corner_procedure_5_1_match'')'
        do test_id=1,2
           
           call make_test_does_corner_procedure_5_1_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_5_1_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test does_corner_procedure_5_2_match
        print '(''test_does_corner_procedure_5_2_match'')'
        do test_id=1,2
           
           call make_test_does_corner_procedure_5_2_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_5_2_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test does_corner_procedure_6_1_match
        print '(''test_does_corner_procedure_6_1_match'')'
        do test_id=1,2
           
           call make_test_does_corner_procedure_6_1_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_6_1_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test does_corner_procedure_6_2_match
        print '(''test_does_corner_procedure_6_2_match'')'
        do test_id=1,2
           
           call make_test_does_corner_procedure_6_2_match(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_match,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           match = does_corner_procedure_6_2_match(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           if(test_match) then
              test_validated = match.eqv.test_match
              test_validated = test_validated.and.
     $             (procedure_type.eq.test_procedure_type)
              test_validated = test_validated.and.
     $             (nb_pts_x.eq.test_nb_pts_x)
              test_validated = test_validated.and.
     $             (nb_pts_y.eq.test_nb_pts_y)

           else
              test_validated = match.eqv.test_match
           end if

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  match         : '',L1,2X,L1)', test_match, match
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)
        
        end do
        print '()'


        !test get_bc_pt_procedure()
        print '(''test_get_bc_pt_procedure'')'
        do test_id=1,44

           call make_test_bc_pt_procedure(
     $          test_id,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_procedure_type,
     $          test_nb_pts_x,
     $          test_nb_pts_y)

           call get_bc_pt_procedure(
     $          test_i,test_j,
     $          grdpts_id,
     $          procedure_type,
     $          nb_pts_x,
     $          nb_pts_y)

           test_validated =
     $          (procedure_type.eq.test_procedure_type)
           test_validated = test_validated.and.
     $          (nb_pts_x.eq.test_nb_pts_x)
           test_validated = test_validated.and.
     $          (nb_pts_y.eq.test_nb_pts_y)

           print '(''test '',I2,'': '',L1)', test_id, test_validated

           if(.not.test_validated) then
              print '(''  procedure_type: '',I2,2X,I2)', test_procedure_type, procedure_type
              print '(''  nb_pts_x      : '',I2,2X,I2)', test_nb_pts_x, nb_pts_x
              print '(''  nb_pts_y      : '',I2,2X,I2)', test_nb_pts_y, nb_pts_y
           end if

           deallocate(grdpts_id)

        end do
        print '()'


        contains

        subroutine make_test_bc_interior_pt_procedure(
     $       test_id,
     $       grdpts_id,
     $       test_procedure_type)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(3,3), intent(out) :: grdpts_id
          integer                , intent(out) :: test_procedure_type


          select case(test_id)

            !  -------
            ! | 1 1 1 |
            ! | 0 0 0 |
            ! |       |
            !  -------
            case(1)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=N_edge_type
               
            !  -------
            ! |       |
            ! | 0 0 0 |
            ! | 1 1 1 |
            !  -------
            case(2)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=S_edge_type

            !  -------
            ! |   0 1 |
            ! |   0 1 |
            ! |   0 1 |
            !  -------
            case(3)
               grdpts_id(1,1) = interior_pt
               grdpts_id(1,2) = interior_pt
               grdpts_id(1,3) = interior_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = bc_pt
               grdpts_id(3,2) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=E_edge_type

            !  -------
            ! | 1 0   |
            ! | 1 0   |
            ! | 1 0   |
            !  -------
            case(4)
               grdpts_id(1,1) = bc_pt
               grdpts_id(1,2) = bc_pt
               grdpts_id(1,3) = bc_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = interior_pt
               grdpts_id(3,2) = interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=W_edge_type


            !  -------
            ! | 1 1 1 |
            ! | 0 0 0 |
            ! | 0     |
            !  -------
            case(5)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(2,1) = interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=N_edge_type
               
            !  -------
            ! |     0 |
            ! | 0 0 0 |
            ! | 1 1 1 |
            !  -------
            case(6)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=S_edge_type

            !  -------
            ! |   0 1 |
            ! |   0 1 |
            ! | 0 0 1 |
            !  -------
            case(7)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(1,2) = interior_pt
               grdpts_id(1,3) = interior_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = bc_pt
               grdpts_id(3,2) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=E_edge_type

            !  -------
            ! | 1 0 0 |
            ! | 1 0   |
            ! | 1 0   |
            !  -------
            case(8)
               grdpts_id(1,1) = bc_pt
               grdpts_id(1,2) = bc_pt
               grdpts_id(1,3) = bc_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = interior_pt
               grdpts_id(3,2) = interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=W_edge_type

            !  -------
            ! | 1 1 1 |
            ! | 1 0 0 |
            ! | 1 0   |
            !  -------
            case(9)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NW_corner_type

            !  -------
            ! | 1 1 1 |
            ! | 0 0 1 |
            ! |   0 1 |
            !  -------
            case(10)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NE_corner_type

            !  -------
            ! | 1 0   |
            ! | 1 0 0 |
            ! | 1 1 1 |
            !  -------
            case(11)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=SW_corner_type

            !  -------
            ! |   0 1 |
            ! | 0 0 1 |
            ! | 1 1 1 |
            !  -------
            case(12)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=SE_corner_type

            !  -------
            ! | 1 1 1 |
            ! | 1 0 0 |
            ! | 0 0   |
            !  -------
            case(13)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NW_corner_type

            !  -------
            ! | 0 1 1 |
            ! | 0 0 1 |
            ! |   0 1 |
            !  -------
            case(14)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = bc_interior_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NE_corner_type

            !  -------
            ! | 0 0   |
            ! | 1 0 0 |
            ! | 1 1 1 |
            !  -------
            case(15)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=SW_corner_type

            !  -------
            ! |   0 0 |
            ! | 0 0 1 |
            ! | 1 1 1 |
            !  -------
            case(16)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=SE_corner_type

            case default
               print '(''test_bf_layer_bc_procedure'')'
               print '(''get_test_bc_interior_pt_procedure'')'
               print '(''test case not implemented'')'
               stop 'choose [1,16]'

          end select

        end subroutine make_test_bc_interior_pt_procedure


        subroutine make_test_is_bc_pt_procedure_of_edge_type(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_procedure_type,
     $     test_edge_type)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer                             , intent(out) :: test_procedure_type
          logical                             , intent(out) :: test_edge_type

          integer :: k

          
          select case(test_id)

            !  -----------
            ! | 2 2 2*2 2 |
            ! | 1 1 1 1 1 |
            ! | 0 0 0 0 0 |
            !  ----------- 
            case(1)
               allocate(grdpts_id(5,3))
               do k=1,5
                  grdpts_id(k,1) = interior_pt
                  grdpts_id(k,2) = bc_interior_pt
                  grdpts_id(k,3) = bc_pt
               end do
               test_i = 3
               test_j = 3
               test_procedure_type = N_edge_type
               test_edge_type      = .true.

            !  ----------- 
            ! |           |
            ! |           |
            ! | 2 2 2*2 2 |
            ! | 2 1 1 1 1 |
            ! | 1 1 0 0 0 |
            !  ----------- 
            case(2)
               allocate(grdpts_id(5,5))
               do k=1,5
                  grdpts_id(k,1) = interior_pt
                  grdpts_id(k,2) = bc_interior_pt
                  grdpts_id(k,3) = bc_pt
                  grdpts_id(k,4) = no_pt
                  grdpts_id(k,5) = no_pt
               end do
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(1,2) = bc_pt               
               test_i = 3
               test_j = 3
               test_procedure_type = N_edge_type
               test_edge_type      = .true.

            !  -----------
            ! | 0 0 0 1 1 |
            ! | 1 1 1 1 2 |
            ! | 2 2 2*2 2 |
            !  ----------- 
            case(3)
               allocate(grdpts_id(5,3))
               do k=1,5
                  grdpts_id(k,1) = bc_pt
                  grdpts_id(k,2) = bc_interior_pt
                  grdpts_id(k,3) = interior_pt
               end do
               grdpts_id(5,2) = bc_pt
               grdpts_id(4,3) = bc_interior_pt
               grdpts_id(5,3) = bc_interior_pt
               test_i = 3
               test_j = 1
               test_procedure_type = S_edge_type
               test_edge_type      = .true.

            !  ----------- 
            ! | 0 0 0 0 0 |
            ! | 1 1 1 1 1 |
            ! | 2 2 2*2 2 |
            ! |         2 |
            ! |           |   
            !  -----------
            case(4)
               allocate(grdpts_id(5,5))
               do k=1,5
                  grdpts_id(k,1) = no_pt
                  grdpts_id(k,2) = no_pt
                  grdpts_id(k,3) = bc_pt
                  grdpts_id(k,4) = bc_interior_pt
                  grdpts_id(k,5) = interior_pt
               end do
               grdpts_id(5,2) = bc_pt
               test_i = 3
               test_j = 3
               test_procedure_type = S_edge_type
               test_edge_type      = .true.

            !  -------
            ! | 2 1 0 |
            ! | 2 1 0 |
            ! | 2*1 0 |
            ! | 2 1 1 |
            ! | 2 2 1 |
            !  -------
            case(5)
               allocate(grdpts_id(3,5))
               do k=1,5
                  grdpts_id(1,k) = bc_pt
                  grdpts_id(2,k) = bc_interior_pt
                  grdpts_id(3,k) = interior_pt
               end do
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt
               test_i = 1
               test_j = 3
               test_procedure_type = W_edge_type
               test_edge_type      = .true.

            !  -----------    
            ! |     2 1 0 |
            ! |     2 1 0 |
            ! |     2*1 0 |
            ! |     2 1 0 |
            ! |   2 2 1 0 |
            !  -----------
            case(6)
               allocate(grdpts_id(5,5))
               do k=1,5
                  grdpts_id(1,k) = no_pt
                  grdpts_id(2,k) = no_pt
                  grdpts_id(3,k) = bc_pt
                  grdpts_id(4,k) = bc_interior_pt
                  grdpts_id(5,k) = interior_pt
               end do
               grdpts_id(2,1)=bc_pt
               test_i = 3
               test_j = 3
               test_procedure_type = W_edge_type
               test_edge_type      = .true.

            !  -------
            ! | 1 2 2 |
            ! | 1 1 2 |
            ! | 0 1 2*|
            ! | 0 1 2 |
            ! | 0 1 2 |
            !  -------
            case(7)
               allocate(grdpts_id(3,5))
               do k=1,5
                  grdpts_id(1,k) = interior_pt
                  grdpts_id(2,k) = bc_interior_pt
                  grdpts_id(3,k) = bc_pt
               end do
               grdpts_id(1,4) = bc_interior_pt
               grdpts_id(1,5) = bc_interior_pt
               grdpts_id(2,5) = bc_pt

               test_i = 3
               test_j = 3
               test_procedure_type = E_edge_type
               test_edge_type      = .true.

            !  -----------    
            ! | 0 1 2 2   |
            ! | 0 1 2     |
            ! | 0 1 2*    |
            ! | 0 1 2     |
            ! | 0 1 2     |
            !  -----------
            case(8)
               allocate(grdpts_id(5,5))
               do k=1,5
                  grdpts_id(1,k) = interior_pt
                  grdpts_id(2,k) = bc_interior_pt
                  grdpts_id(3,k) = bc_pt
                  grdpts_id(4,k) = no_pt
                  grdpts_id(5,k) = no_pt
               end do
               grdpts_id(4,5)=bc_pt
               test_i = 3
               test_j = 3
               test_procedure_type = E_edge_type
               test_edge_type      = .true.

           !  -----------
           ! | 2 2 2*2   |
           ! | 1 1 1 2 2 |
           ! | 0 0 1 1 1 |
           !  ----------- 
           case(9)
              allocate(grdpts_id(5,3))
              do k=1,2
                 grdpts_id(k,1) = interior_pt
              end do
              do k=3,5
                 grdpts_id(k,1) = bc_interior_pt
              end do
              do k=1,3
                 grdpts_id(k,2) = bc_interior_pt
              end do
              do k=4,5
                 grdpts_id(k,2) = bc_pt
              end do
              do k=1,4
                 grdpts_id(k,3) = bc_pt
              end do
              grdpts_id(5,3) = no_pt
              test_i = 3
              test_j = 3
              test_procedure_type = N_edge_type
              test_edge_type      = .false.

            !  ----------- 
            ! | 0 0 0 0 0 |
            ! | 1 1 1 1 1 |
            ! | 1 2 2*2 2 |
            ! | 2 2     2 |
            ! |           |   
            !  -----------
            case(10)
               allocate(grdpts_id(5,5))
               do k=1,5
                  grdpts_id(k,1) = no_pt
               end do

               do k=1,2
                  grdpts_id(k,2) = bc_pt
               end do
               do k=3,4
                  grdpts_id(k,2) = no_pt
               end do
               grdpts_id(5,2) = bc_pt

               grdpts_id(1,3) = bc_interior_pt
               do k=2,5
                  grdpts_id(k,3) = bc_pt
               end do

               do k=1,5
                  grdpts_id(k,4) = bc_interior_pt
               end do

               do k=1,5
                  grdpts_id(k,5) = interior_pt
               end do
               test_i = 3
               test_j = 3
               test_procedure_type = S_edge_type
               test_edge_type      = .false.

            !  -------
            ! | 2 1 0 |
            ! | 2 1 0 |
            ! | 2*1 1 |
            ! | 2 2 2 |
            ! |       |
            !  -------
            case(11)
               allocate(grdpts_id(3,5))
               do k=1,3
                  grdpts_id(k,1) = no_pt
               end do
               do k=1,3
                  grdpts_id(k,2) = bc_pt
               end do
               grdpts_id(1,3) = bc_pt
               do k=2,3
                  grdpts_id(k,3) = bc_interior_pt
               end do
               do k=4,5
                  grdpts_id(1,k) = bc_pt
                  grdpts_id(2,k) = bc_interior_pt
                  grdpts_id(3,k) = interior_pt
               end do
               test_i = 1
               test_j = 3
               test_procedure_type = W_edge_type
               test_edge_type      = .false.

            !  -----------    
            ! | 0 1 2 2   |
            ! | 0 1 2     |
            ! | 1 1 2*    |
            ! | 2 2 2     |
            ! | 2         |
            !  -----------
            case(12)
               allocate(grdpts_id(5,5))
               grdpts_id(1,1) = bc_pt
               do k=4,5
                  grdpts_id(k,1) = no_pt
               end do
               do k=1,3
                  grdpts_id(k,2) = bc_pt
               end do
               do k=4,5
                  grdpts_id(k,2) = no_pt
               end do
               do k=1,2
                  grdpts_id(k,3) = bc_interior_pt
               end do
               grdpts_id(3,3) = bc_pt
               do k=4,5
                  grdpts_id(k,3) = no_pt
               end do
               grdpts_id(1,4) = interior_pt
               grdpts_id(2,4) = bc_interior_pt
               grdpts_id(3,4) = bc_pt
               do k=4,5
                  grdpts_id(k,4) = no_pt
               end do
               grdpts_id(1,5) = interior_pt
               grdpts_id(2,5) = bc_interior_pt
               do k=3,4
                  grdpts_id(k,5) = bc_pt
               end do
               grdpts_id(5,5) = no_pt
               test_i = 3
               test_j = 3
               test_procedure_type = E_edge_type
               test_edge_type      = .false.

            case default
               print '( ''test_bf_layer_bc_procedure'')'
               print '(''make_test_is_bc_pt_procedure_of_edge_type'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_is_bc_pt_procedure_of_edge_type


        subroutine make_test_does_corner_procedure_1_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! | 0 1 1 1 1 |
          ! | 0 1 2 2 2 |
          ! | 1 1 2*    |
          ! | 1 2 2     |
          ! | 1 2       |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! |           |
          ! |           |
          ! |     2*2 2 |
          ! | 2 2 2 1 1 |
          ! | 1 1 1 1 0 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |           |
          ! |     2 2 2 |
          ! |     2*1 1 |
          ! | 2 2 2 1 0 |
          ! | 1 1 1 1 0 |
          !  -----------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! |           |
          ! |     2 2 2 |
          ! |     2*1 1 |
          ! | 2 2 2 1 0 |
          ! | 1 1 1 1 0 |
          !  -----------
            case(4)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! |       2 1 |
          ! |       2 1 |
          ! |   2 2*2 1 |
          ! |   2 1 1 1 |
          ! |   2 1 0 0 |
          !  -----------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,interior_pt,interior_pt,
     $              no_pt,bc_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 1 2     |
          ! | 1 1 2     |
          ! | 1 2 2*    |
          ! | 1 2       |
          ! | 1 2       |
          !  -----------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 1 1 1 |
          ! | 0 1 1 2 2 |
          ! | 1 2 2*2   |
          ! | 2 2       |
          ! |           |
          !  -----------
            case(7)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,bc_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,
     $              interior_pt,interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 0 1 2 |
          ! | 1 1 1 1 2 |
          ! | 1 2 2*2 2 |
          ! | 1 2       |
          ! | 2 2       |
          !  -----------
            case(8)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,interior_pt,interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = S_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_1_match'')'
               print '(''make_test_does_corner_procedure_1_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_1_match


        subroutine make_test_does_corner_procedure_2_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! | 1 1 1 1 0 |
          ! | 2 2 2 1 0 |
          ! |     2*1 1 |
          ! |     2 2 2 |
          ! |         2 |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! |           |
          ! | 2         |
          ! | 2 2 2*    |
          ! | 1 1 2 2   |
          ! | 0 1 1 2 2 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |           |
          ! | 2 2 2     |
          ! | 1 1 2*    |
          ! | 0 1 2 2   |
          ! | 0 1 1 2   |
          !  -----------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! | 0 1 2 2 2 |
          ! | 0 1 2     |
          ! | 0 1 2*    |
          ! | 0 1 2 2   |
          ! | 0 1 1 2   |
          !  -----------
            case(4)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = E_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |           |
          ! | 2 2       |
          ! | 1 2 2*2   |
          ! | 1 1 1 2 2 |
          ! | 0 0 1 1 1 |
          !  -----------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt,bc_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! | 1 1 1 1 0 |
          ! | 2 2 2 1 1 |
          ! |     2*2 1 |
          ! |       2 1 |
          ! |       2 1 |
          !  -----------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              no_pt, no_pt, no_pt, bc_pt, bc_interior_pt,
     $              no_pt, no_pt, no_pt, bc_pt, bc_interior_pt,
     $              no_pt, no_pt, bc_pt, bc_pt, bc_interior_pt,
     $              bc_pt, bc_pt, bc_pt, bc_interior_pt, bc_interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt, interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 1 1 1 0 0 |
          ! | 2 2 1 1 1 |
          ! |   2 2*2 1 |
          ! |       2 1 |
          ! |       2 1 |
          !  -----------
            case(7)
               grdpts_id = reshape(
     $              (/
     $              no_pt, no_pt, no_pt, bc_pt, bc_interior_pt,
     $              no_pt, no_pt, no_pt, bc_pt, bc_interior_pt,
     $              no_pt, bc_pt, bc_pt, bc_pt, bc_interior_pt,
     $              bc_pt, bc_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, interior_pt, interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 0 0 0 |
          ! | 1 1 1 1 1 |
          ! | 2 2 2*2 1 |
          ! |       2 1 |
          ! |       2 2 |
          !  -----------
            case(8)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = S_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_2_match'')'
               print '(''make_test_does_corner_procedure_2_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_2_match


        subroutine make_test_does_corner_procedure_3_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! | 1 1 1 0 0 |
          ! | 2 2 1 1 0 |
          ! |   2 2*1 0 |
          ! |     2 1 0 |
          ! |     2 1 0 |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! | 2 2       |
          ! | 1 2       |
          ! | 1 2 2*    |
          ! | 1 1 2 2 2 |
          ! | 0 1 1 1 1 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 1 2       |
          ! | 1 2       |
          ! | 1 2 2*2 2 |
          ! | 1 1 1 1 1 |
          ! | 0 0 0 0 0 |
          !  -----------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = N_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

c$$$          !  -----------    
c$$$          ! | 2 2       |
c$$$          ! | 1 2 2     |
c$$$          ! | 1 1 2*2   |
c$$$          ! | 0 1 1 2   |
c$$$          ! | 0 0 1 2   |
c$$$          !  -----------
          !  -----------    
          ! | 2 2       |
          ! | 1 2 2     |
          ! | 1 1 2*    |
          ! | 0 1 2     |
          ! | 0 1 2     |
          !  -----------
            case(4)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,no_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_pt,bc_pt,no_pt,no_pt,
c$$$     $              bc_pt,bc_pt,no_pt,no_pt,no_pt
c$$$     $              /),
c$$$     $              (/5,5/))
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! | 1 1 1 1 0 |
          ! | 2 2 2 1 1 |
          ! |     2*2 2 |
          ! |           |
          ! |           |
          !  -----------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 1 1 1 1 0 |
          ! | 2 2 2 1 0 |
          ! |     2*1 0 |
          ! |     2 1 0 |
          ! |     2 1 0 |
          !  -----------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              no_pt, no_pt, bc_pt, bc_interior_pt, interior_pt,
     $              no_pt, no_pt, bc_pt, bc_interior_pt, interior_pt,
     $              no_pt, no_pt, bc_pt, bc_interior_pt, interior_pt,
     $              bc_pt, bc_pt, bc_pt, bc_interior_pt, interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt, interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = W_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0


            case default
               print '( ''test_does_corner_procedure_3_match'')'
               print '(''make_test_does_corner_procedure_3_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_3_match


        subroutine make_test_does_corner_procedure_4_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! | 0 0 1 1 1 |
          ! | 0 1 1 2 2 |
          ! | 0 1 2*2   |
          ! | 0 1 2     |
          ! | 0 1 2     |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,
     $              interior_pt,interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_corner_type
               test_nb_pts_x       = 1
               test_nb_pts_y       = 0

          !  -----------    
          ! |       2 2 |
          ! |       2 1 |
          ! |     2*2 1 |
          ! |     2 1 1 |
          ! |     2 1 0 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |           |
          ! |       2 2 |
          ! | 2 2 2*2 1 |
          ! | 1 1 1 1 1 |
          ! | 0 0 0 0 0 |
          !  -----------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = N_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |       2 1 |
          ! |     2 2 1 |
          ! |     2*1 1 |
          ! |     2 1 0 |
          ! |     2 1 0 |
          !  -----------
            case(4)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,no_pt,no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 1

          !  -----------    
          ! | 0 1 1 2   |
          ! | 1 1 2 2   |
          ! | 2 2 2*    |
          ! |           |
          ! |           |
          !  -----------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,no_pt,
     $              interior_pt, bc_interior_pt, bc_interior_pt, bc_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = Se_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 1 1 1 2 |
          ! | 0 1 2 2 2 |
          ! | 0 1 2*    |
          ! | 0 1 2     |
          ! | 0 1 2     |
          !  -----------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = E_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0


            case default
               print '( ''test_does_corner_procedure_4_match'')'
               print '(''make_test_does_corner_procedure_4_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_4_match


        subroutine make_test_does_corner_procedure_5_1_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! |           |
          ! |           |
          ! | 2 2 2*    |
          ! | 1 1 2 2 2 |
          ! | 0 1 1 1 1 |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 0 0 0 |
          ! | 1 1 1 1 0 |
          ! | 2 2 2*1 0 |
          ! |     2 1 0 |
          ! |     2 1 0 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_5_1_match'')'
               print '(''make_test_does_corner_procedure_5_1_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_5_1_match


        subroutine make_test_does_corner_procedure_5_2_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! |         2 |
          ! |         2 |
          ! |     2*2 2 |
          ! |   2 2 1 1 |
          ! |   2 1 1 0 |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              no_pt,bc_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,bc_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 0 0 0 |
          ! | 0 1 1 1 1 |
          ! | 0 1 2*2 2 |
          ! | 0 1 2     |
          ! | 0 1 2     |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_5_2_match'')'
               print '(''make_test_does_corner_procedure_5_2_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_5_2_match

      
        subroutine make_test_does_corner_procedure_6_1_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! | 0 1 2     |
          ! | 1 1 2     |
          ! | 2 2 2*    |
          ! |           |
          ! |           |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SE_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! |     2 1 0 |
          ! |     2 1 0 |
          ! | 2 2 2*1 0 |
          ! | 1 1 1 1 0 |
          ! | 0 0 0 0 0 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NW_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_6_1_match'')'
               print '(''make_test_does_corner_procedure_6_1_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_6_1_match


        subroutine make_test_does_corner_procedure_6_2_match(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_match,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          logical                             , intent(out) :: test_match
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y


          test_i = 3
          test_j = 3

          allocate(grdpts_id(5,5))

          
          select case(test_id)

          !  -----------    
          ! |     2 1 0 |
          ! |     2 1 1 |
          ! |     2*2 2 |
          ! |           |
          ! |           |
          !  -----------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = SW_corner_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

          !  -----------    
          ! | 0 0 2     |
          ! | 0 1 2     |
          ! | 0 1 2*2 2 |
          ! | 0 1 1 1 1 |
          ! | 0 0 0 0 0 |
          !  -----------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt,bc_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt
     $              /),
     $              (/5,5/))
               test_match          = .true.
               test_procedure_type = NE_edge_type
               test_nb_pts_x       = 0
               test_nb_pts_y       = 0

            case default
               print '( ''test_does_corner_procedure_6_2_match'')'
               print '(''make_test_does_corner_procedure_6_2_match'')'
               print '(''test case not implemented: '',I2)', test_id
               stop 'not implemented'

          end select

        end subroutine make_test_does_corner_procedure_6_2_match


        subroutine make_test_bc_pt_procedure(
     $     test_id,
     $     grdpts_id,
     $     test_i,
     $     test_j,
     $     test_procedure_type,
     $     test_nb_pts_x,
     $     test_nb_pts_y)

          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer                             , intent(out) :: test_procedure_type
          integer                             , intent(out) :: test_nb_pts_x
          integer                             , intent(out) :: test_nb_pts_y

          logical :: test_match
          logical :: test_edge_type
          integer :: test_id2

          
          select case(test_id)


            !edge tests
            case(1:8)

               test_id2 = test_id
               
               call make_test_is_bc_pt_procedure_of_edge_type(
     $              test_id,
     $              grdpts_id,
     $              test_i,
     $              test_j,
     $              test_procedure_type,
     $              test_edge_type)
               
               test_nb_pts_x = 0
               test_nb_pts_y = 0

            !corner 1 tests
            case(9:16)

                test_id2 = test_id-8
                
                call make_test_does_corner_procedure_1_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 2 tests
            case(17:24)
               
               test_id2 = test_id-16

               call make_test_does_corner_procedure_2_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 3 tests
            case(25:30)
               
               test_id2 = test_id-24

               call make_test_does_corner_procedure_3_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 4 tests
            case(31:36)
               
               test_id2 = test_id-30

               call make_test_does_corner_procedure_4_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 5_1 tests
            case(37:38)
               
               test_id2 = test_id-36

               call make_test_does_corner_procedure_5_1_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 5_2 tests
            case(39:40)
               
               test_id2 = test_id-38

               call make_test_does_corner_procedure_5_2_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)


            !corner 6_1 tests
            case(41:42)
               
               test_id2 = test_id-40

               call make_test_does_corner_procedure_6_1_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            !corner 5_2 tests
            case(43:44)
               
               test_id2 = test_id-42

               call make_test_does_corner_procedure_6_2_match(
     $               test_id2,
     $               grdpts_id,
     $               test_i,
     $               test_j,
     $               test_match,
     $               test_procedure_type,
     $               test_nb_pts_x,
     $               test_nb_pts_y)

            case default
               
               print '(''test_bf_layer_bc_procedure'')'
               print '(''make_test_bc_pt_procedure'')'
               print '(''case not recognized: '',I2)',test_id
               stop 'case not implemented'

           end select

        end subroutine make_test_bc_pt_procedure

      end program test_bf_layer_bc_procedure
