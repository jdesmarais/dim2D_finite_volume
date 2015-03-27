      program test_geometry_update

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_validated

        use field_extended_class, only :
     $       field_extended

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,
     $       x_min,x_max,
     $       y_min,y_max,
     $       debug_adapt_computational_domain,
     $       debug_geometry_update

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        integer, parameter :: spot_transported=1
        integer, parameter :: bubble_expansion=2

        integer, parameter :: main_case_choice = bubble_expansion
        integer, parameter :: sub_case_choice = 1 

        call test_geometry_update_on_case()


        contains


        subroutine test_geometry_update_on_case()

          implicit none

          type(field_extended) :: f_simulated !  field simulated
          integer(ikind)       :: nt
          integer(ikind)       :: t
          real(rkind ) :: dt
          real(rkind), dimension(nx,ny,ne) :: nodes0


          call check_inputs()


          ! initialize the field
          call f_simulated%ini()
          call overwrite_nodes(f_simulated,0)
          call f_simulated%write_data()
  
  
          nt = 35
  
          ! adapt the field
          do t=1, nt
  
             ! overwrite the nodes
             call overwrite_nodes(f_simulated,t)

             ! adapt the domain
             call f_simulated%domain_extension%initialize_before_timeInt(
     $            f_simulated%interior_bc_sections)
             call f_simulated%adapt_domain(dt,nodes0)
             call f_simulated%domain_extension%finalize_after_timeInt()
  
             !  write the output data
             call f_simulated%write_data()
  
          end do

        end subroutine test_geometry_update_on_case


        subroutine overwrite_nodes(f_simulated,t)

          implicit none

          class(field_extended), intent(inout) :: f_simulated
          integer              , intent(in)    :: t

          select case(main_case_choice)
            case(spot_transported)
               call overwrite_nodes_with_spot_transported(
     $              f_simulated,t)
            case(bubble_expansion)
               call overwrite_nodes_with_bubble_expansion(
     $              f_simulated,t)
          end select

        end subroutine overwrite_nodes


        subroutine overwrite_nodes_with_spot_transported(f_simulated,t)

          implicit none

          class(field_extended), intent(inout) :: f_simulated
          integer              , intent(in)    :: t

          

          real(rkind) :: x_center
          real(rkind) :: y_center

          integer :: k
          integer :: m
          integer :: nb_sublayers
          type(bf_sublayer), pointer :: bf_sublayer_ptr

          real(rkind) :: radius
          real(rkind) :: ux
          real(rkind) :: uy


          radius = 5.0d0

          
          select case(sub_case_choice)
            !right
            case(1)
               ux = 1.0d0
               uy = 0.0d0

            !left
            case(2)
               ux =-1.0d0
               uy = 0.0d0

            !top
            case(3)
               ux = 0.0d0
               uy = 1.0d0

            !bottom
            case(4)
               ux = 0.0d0
               uy =-1.0d0

            !NE_corner
            case(5)
               ux = 1.0d0
               uy = 1.0d0

            !NW corner
            case(6)
               ux =-1.0d0
               uy = 1.0d0

            !SE_corner
            case(7)
               ux = 1.0d0
               uy =-1.0d0

            !SW corner
            case(8)
               ux =-1.0d0
               uy =-1.0d0
          end select


          x_center = 0.0d0 + ux*t
          y_center = 0.0d0 + uy*t

          
          !interior nodes
          call apply_spot_transported(
     $         f_simulated%x_map,
     $         f_simulated%y_map,
     $         f_simulated%nodes,
     $         x_center,
     $         y_center,
     $         radius)

          !buffer layers
          do k=1,4
             
             nb_sublayers = f_simulated%domain_extension%mainlayer_pointers(k)%get_nb_sublayers()

             if(nb_sublayers.gt.0) then

                bf_sublayer_ptr => f_simulated%domain_extension%mainlayer_pointers(k)%get_head_sublayer()

                do m=1, nb_sublayers
                   
                   call apply_spot_transported(
     $                  bf_sublayer_ptr%x_map,
     $                  bf_sublayer_ptr%y_map,
     $                  bf_sublayer_ptr%nodes,
     $                  x_center,
     $                  y_center,
     $                  radius)

                   bf_sublayer_ptr => bf_sublayer_ptr%get_next()
                   
                end do

             end if

          end do

        end subroutine overwrite_nodes_with_spot_transported


        subroutine overwrite_nodes_with_bubble_expansion(f_simulated,t)

          implicit none

          class(field_extended), intent(inout) :: f_simulated
          integer              , intent(in)    :: t

          

          real(rkind) :: x_center
          real(rkind) :: y_center

          integer :: k
          integer :: m
          integer :: nb_sublayers
          type(bf_sublayer), pointer :: bf_sublayer_ptr

          real(rkind) :: radius1
          real(rkind) :: radius2
          real(rkind) :: ux
          real(rkind) :: uy


          radius1 = 1.0d0 + t
          radius2 = 5.0d0 + t

          x_center = 0.0d0
          y_center = 0.0d0

          
          !interior nodes
          call apply_bubble_expansion(
     $         f_simulated%x_map,
     $         f_simulated%y_map,
     $         f_simulated%nodes,
     $         x_center,
     $         y_center,
     $         radius1,
     $         radius2)

          !buffer layers
          do k=1,4
             
             nb_sublayers = f_simulated%domain_extension%mainlayer_pointers(k)%get_nb_sublayers()

             if(nb_sublayers.gt.0) then

                bf_sublayer_ptr => f_simulated%domain_extension%mainlayer_pointers(k)%get_head_sublayer()

                do m=1, nb_sublayers
                   
                   call apply_bubble_expansion(
     $                  bf_sublayer_ptr%x_map,
     $                  bf_sublayer_ptr%y_map,
     $                  bf_sublayer_ptr%nodes,
     $                  x_center,
     $                  y_center,
     $                  radius1,
     $                  radius2)

                   bf_sublayer_ptr => bf_sublayer_ptr%get_next()
                   
                end do

             end if

          end do

        end subroutine overwrite_nodes_with_bubble_expansion


        subroutine apply_spot_transported(
     $     x_map,
     $     y_map,
     $     nodes,
     $     x_center,
     $     y_center,
     $     radius)

           implicit none

           real(rkind), dimension(:)    , intent(in)  :: x_map
           real(rkind), dimension(:)    , intent(in)  :: y_map
           real(rkind), dimension(:,:,:), intent(out) :: nodes
           real(rkind)                  , intent(in)  :: x_center
           real(rkind)                  , intent(in)  :: y_center
           real(rkind)                  , intent(in)  :: radius
           
           
           integer(ikind) :: i,j
           real(rkind)    :: distance
           

           do j=1, size(y_map,1)
              do i=1, size(x_map,1)

                 distance = SQRT((x_map(i)-x_center)**2+(y_map(j)-y_center)**2)

                 if(distance.gt.radius) then
                    nodes(i,j,1) =  1.0d0
                 else
                    nodes(i,j,1) = -1.0d0
                 end if

              end do
           end do

        end subroutine apply_spot_transported


        subroutine apply_bubble_expansion(
     $     x_map,
     $     y_map,
     $     nodes,
     $     x_center,
     $     y_center,
     $     radius1,
     $     radius2)

           implicit none

           real(rkind), dimension(:)    , intent(in)  :: x_map
           real(rkind), dimension(:)    , intent(in)  :: y_map
           real(rkind), dimension(:,:,:), intent(out) :: nodes
           real(rkind)                  , intent(in)  :: x_center
           real(rkind)                  , intent(in)  :: y_center
           real(rkind)                  , intent(in)  :: radius1
           real(rkind)                  , intent(in)  :: radius2
           
           
           integer(ikind) :: i,j
           real(rkind)    :: distance
           

           do j=1, size(y_map,1)
              do i=1, size(x_map,1)

                 distance = SQRT((x_map(i)-x_center)**2+(y_map(j)-y_center)**2)

                 if(distance.gt.radius2) then
                    nodes(i,j,1) =  1.0d0
                 else
                    if(distance.gt.radius1) then
                       nodes(i,j,1) = -1.0d0
                    else
                       nodes(i,j,1) = 1.0d0
                    end if
                 end if

              end do
           end do

        end subroutine apply_bubble_expansion


        subroutine check_inputs()

          implicit none


          if(.not.(
     $         (is_real_validated(x_min,-10.0d0,.false.)).and.
     $         (is_real_validated(x_max, 10.0d0,.false.)).and.
     $         (is_real_validated(y_min,-15.0d0,.false.)).and.
     $         (is_real_validated(y_max, 15.0d0,.false.)).and.
     $         (nx.eq.25).and.
     $         (ny.eq.35).and.
     $         (debug_adapt_computational_domain.eqv.(.true.)).and.
     $         (debug_geometry_update.eqv.(.true.)))) then

             print '(''the test requires:'')'
             print '(''  - x_min = -10.0d0'')'
             print '(''  - x_max =  10.0d0'')'
             print '(''  - y_min = -10.0d0'')'
             print '(''  - y_max =  15.0d0'')'
             print '(''  - nx    =  25'')'
             print '(''  - ny    =  35'')'
             print '(''  - debug_adapt_computational_domain=.true.'')'
             print '(''  - debug_geometry_update=.true.'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_geometry_update
