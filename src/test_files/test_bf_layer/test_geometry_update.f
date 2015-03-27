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
        integer            :: case_choice


        case_choice = spot_transported

        call test_geometry_update_on_case()


        contains


        subroutine test_geometry_update_on_case()

          implicit none

          type(field_extended) :: f_simulated !  field simulated
          integer(ikind)       :: nt
          integer(ikind)       :: t
          

          ! initialize the field
          call overwrite_nodes(f_simulated,0)
          call f_simulated%write_data()
  
  
          nt = 20
  
          ! adapt the field
          do t=1, nt
  
             ! overwrite the nodes
             call overwrite_nodes(f_simulated,t)
  
             ! adapt the domain
             call f_simulated%integrate(0.1)
  
             !  write the output data
             call f_simulated%write_data()
  
          end do

        end subroutine test_geometry_update_on_case


        subroutine overwrite_nodes(f_simulated,t)

          implicit none

          class(field_extended), intent(inout) :: f_simulated
          integer              , intent(in)    :: t

          select case(case_choice)
            case(spot_transported)
               call overwrite_nodes_with_spot_transported(
     $              f_simulated,t)
          end select

        end subroutine overwrite_nodes


        subroutine overwrite_nodes_with_spot_transported(f_simulated,t)

          implicit none

          class(field_extended), intent(inout) :: f_simulated
          integer              , intent(in)    :: t

          real(rkind) :: radius
          real(rkind) :: ux
          real(rkind) :: uy

          real(rkind) :: x_center
          real(rkind) :: y_center

          integer :: k
          integer :: m
          integer :: nb_sublayers
          type(bf_sublayer), pointer :: bf_sublayer_ptr


          radius = 5.0d0
          ux = 1.0d0
          uy = 0.0d0

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

             bf_sublayer_ptr => f_simulated%domain_extension%mainlayer_pointers(k)%get_head_sublayer()

             do m=1, nb_sublayers

                call apply_spot_transported(
     $               bf_sublayer_ptr%x_map,
     $               bf_sublayer_ptr%y_map,
     $               bf_sublayer_ptr%nodes,
     $               x_center,
     $               y_center,
     $               radius)

                bf_sublayer_ptr => bf_sublayer_ptr%get_next()
                
             end do

          end do

        end subroutine overwrite_nodes_with_spot_transported


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


        subroutine check_inputs()

          implicit none


          if(.not.(
     $         (is_real_validated(x_min,-10.0d0,.false.)).and.
     $         (is_real_validated(x_max, 10.0d0,.false.)).and.
     $         (is_real_validated(y_min,-15.0d0,.false.)).and.
     $         (is_real_validated(y_max, 15.0d0,.false.)).and.
     $         (nx.eq.25).and.
     $         (ny.eq.35).and.
     $         (debug_adapt_computational_domain.eqv.(.false.)).and.
     $         (debug_geometry_update.eqv.(.true.)))) then

             print '(''the test requires:'')'
             print '(''  - x_min = -10.0d0'')'
             print '(''  - x_max =  10.0d0'')'
             print '(''  - y_min = -15.0d0'')'
             print '(''  - y_max =  15.0d0'')'
             print '(''  - nx    =  25'')'
             print '(''  - ny    =  35'')'
             print '(''  - debug_adapt_computational_domain=.true.'')'
             print '(''  - debug_geometry_update=.true.'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_geometry_update
