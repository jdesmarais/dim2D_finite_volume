      program test_bf_sublayers_merge_prog

        use bf_sublayer_class          , only : bf_sublayer

        use bf_sublayers_merge_module  , only : merge_sublayers_N,
     $                                          merge_sublayers_S,
     $                                          merge_sublayers_E,
     $                                          merge_sublayers_W
        use interface_abstract_class   , only : interface_abstract
        
        use parameters_constant        , only : N,S,E,W
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : ikind, rkind
                                       
        use test_bf_layer_module       , only : print_nodes,
     $                                          print_sizes,
     $                                          ini_nodes
        use test_cases_interface_module, only : ini_interface
        !use test_cases_path_module     , only : ini_path
        use interface_print_module     , only : print_interface

        implicit none

        !parameters for the test
        integer, parameter :: no_alignment = 0
        integer, parameter :: test_case_id = 3
        integer, parameter :: mainlayer_id = 4
        integer, parameter :: nb_sublayers = 2
        integer, parameter :: size_case = 4
        integer, parameter :: distance_case = 4
        integer, parameter :: merge_id = 1
        integer, parameter :: merge_inverse = 1
        integer, parameter :: random_seed = 86456
        integer, parameter :: neighbor_case = 4

        integer, parameter :: interface_before = 0
        integer, parameter :: interface_after = 1


        !local variables
        real(rkind), dimension(nx,ny,ne)   :: nodes
        type(interface_abstract)           :: interface_used
        integer, dimension(2,nb_sublayers) :: relative_sizes
        integer, dimension(nb_sublayers+1) :: relative_distance
        integer(ikind), dimension(2,2)     :: alignment
        logical                            :: neighbor_log1
        logical                            :: neighbor_log2

        integer :: i

        !variables for the test
        type(bf_sublayer), pointer :: sublayer1
        type(bf_sublayer), pointer :: sublayer2
        type(bf_sublayer), pointer :: merged_sublayer


        !verify the test parameters
        call verify_test_parameters()


        !initialize the nodes for the test
        call ini_nodes(nodes)
        call print_sizes(nodes,'interior_sizes.dat')
        call print_nodes(nodes,'interior_nodes.dat')


        !initialize the relative_size and relative_distance
        call ini_relative_sizes(
     $       size_case,
     $       random_seed,
     $       relative_sizes)

        call ini_relative_distance(
     $       distance_case,
     $       random_seed,
     $       relative_distance)

        call ini_alignment(
     $       alignment,
     $       mainlayer_id,
     $       merge_id,
     $       relative_sizes,
     $       relative_distance)

        call ini_neighbors(
     $       neighbor_case,
     $       neighbor_log1, neighbor_log2)


        !initialize the interface with two sublayers
        !in the same mainlayer
        call ini_interface(
     $       interface_used,
     $       test_case_id,
     $       nodes,
     $       mainlayer_id=mainlayer_id,
     $       nb_sublayers=nb_sublayers,
     $       relative_distance=relative_distance,
     $       relative_sizes=relative_sizes)
        

        !select the two sublayers for the merge
        if(merge_inverse.eq.0) then
           sublayer1 => interface_used%mainlayer_pointers(mainlayer_id)%ptr%head_sublayer
           do i=1, merge_id-1
              sublayer1 => sublayer1%next
           end do
           sublayer2 => sublayer1%next
        else
           sublayer2 => interface_used%mainlayer_pointers(mainlayer_id)%ptr%head_sublayer
           do i=1, merge_id-1
              sublayer2 => sublayer2%next
           end do
           sublayer1 => sublayer2%next
        end if


        !print interface before
        call print_interface(interface_used, interface_before)


        !merge the sublayers
        select case(mainlayer_id)

          !north main layer
          case(N)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_N(
     $               interface_used%mainlayer_pointers(N)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_N(
     $               interface_used%mainlayer_pointers(N)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(N)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(N)%ptr%nb_sublayers-1

          !south main layer
          case(S)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_S(
     $               interface_used%mainlayer_pointers(S)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_S(
     $               interface_used%mainlayer_pointers(S)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(S)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(S)%ptr%nb_sublayers-1

          !east main layer
          case(E)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_E(
     $               interface_used%mainlayer_pointers(E)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_E(
     $               interface_used%mainlayer_pointers(E)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(E)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(E)%ptr%nb_sublayers-1

          !west main layer
          case(W)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_W(
     $               interface_used%mainlayer_pointers(W)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_W(
     $               interface_used%mainlayer_pointers(W)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(W)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(W)%ptr%nb_sublayers-1

          case default
             print '(''mainlayer not recognized'')'
             print '(''mainlayer_id: '',I2)', mainlayer_id
             stop 'test not implemented yet'
        end select



        !print interface after
        call print_interface(interface_used, interface_after)


        contains

        subroutine verify_test_parameters()

          if((no_alignment.ne.0).and.(no_alignment.ne.1)) then
             print '(''no_alignment should equal 0 or 1'')'
             print '(''it determines whether the alignment should'')'
             print '(''be taken into account'')'
             stop 'change no_alignment'
          end if
          
          if(test_case_id.ne.3) then
             print '(''only one initialization is implemented'')'
             print '(''for this test'')'
             print '(''please select test_case_id = 3'')'
             stop 'change test_case_id'
          end if

          if(  (mainlayer_id.ne.N).and.
     $         (mainlayer_id.ne.S).and.
     $         (mainlayer_id.ne.E).and.
     $         (mainlayer_id.ne.W)) then
             print '(''mainlayer_id not recognized'')'
             print '(''please select N, S, E, or W'')'
             stop 'change mainlayer_id'
          end if

          if((size_case.le.0).or.(size_case.ge.5)) then
             print '(''size_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change size_case'
          end if

          if((distance_case.le.0).or.(distance_case.ge.5)) then
             print '(''distance_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change distance_case'
          end if

          if(merge_id.gt.(nb_sublayers-1)) then
             print '(''merge_id not recognized'')'
             print '(''please select merge_id<nb_sublayers'')'
             stop 'change merge_id'
          end if

          if((merge_inverse.ne.0).and.(merge_inverse.ne.1)) then
             print '(''merge_inverse not recognized'')'
             print '(''please select 0 or 1'')'
             stop 'change merge_inverse'
          end if

          if((neighbor_case.le.0).or.(neighbor_case.ge.5)) then
             print '(''neighbor_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change neighbor_case'
          end if

        end subroutine verify_test_parameters


        subroutine ini_relative_sizes(
     $     size_case,
     $     seed,
     $     relative_sizes)

          implicit none
          
          integer                , intent(in)  :: size_case
          integer                , intent(in)  :: seed
          integer, dimension(:,:), intent(out) :: relative_sizes


          integer :: random_value
          integer :: i,j


          select case(size_case)

            !they all have the same size=1
            case(1)
               do j=1, size(relative_sizes,2)
                  do i=1,2                   
                     relative_sizes(i,j) = 0
                  end do
               end do

            !they all have the same size
            !determined by a random generator
            case(2)
               call srand(seed)
               random_value = nint(5.0*RAND())

               do j=1, size(relative_sizes,2)
                  relative_sizes(1,j) = random_value
                  relative_sizes(2,j) = random_value-1
               end do

            !they have different sizes in increasing order
            !starting from a random value
            case(3)
               call srand(seed)
               random_value = nint(5.0*RAND())

               do j=1, size(relative_sizes,2)
                  relative_sizes(1,j) = random_value+(j-1)
                  relative_sizes(2,j) = random_value+(j-1)-1
               end do

            !they have different sizes in decreasing order
            !starting from a random value
            case(4)
               call srand(seed)
               random_value = nint(5.0*RAND())

               do j=1, size(relative_sizes,2)
                  relative_sizes(1,j) = random_value-(j-1)
                  relative_sizes(2,j) = random_value-(j-1)-1
               end do

          end select

        end subroutine ini_relative_sizes


        subroutine ini_relative_distance(
     $     distance_case,
     $     seed,
     $     relative_distance)

          implicit none
          
          integer              , intent(in)  :: distance_case
          integer              , intent(in)  :: seed
          integer, dimension(:), intent(out) :: relative_distance

          integer :: random_value

          select case(distance_case)

            !the distance between the different
            !sublayers is fixed to 0
            case(1)
               do i=1, size(relative_distance,1)
                  relative_distance(i) = 0
               end do
               
            !the distance between the different
            !sublayers is fixed to 1
            case(2)
               do i=1, size(relative_distance,1)
                  relative_distance(i) = 1
               end do

            !the distance between the different
            !sublayers is fixed to 2
            case(3)
               do i=1, size(relative_distance,1)
                  relative_distance(i) = 2
               end do

            !the distance between the different
            !sublayers is fixed by a random value
            case(4)
               call srand(seed)
               random_value = nint(5.0*RAND())

               do i=1, size(relative_distance,1)
                  relative_distance(i) = random_value
               end do

          end select

          !if(no_alignment.eq.0) then
          !   relative_distance(1)=3
          !end if

c$$$          relative_distance(merge_id)=0
c$$$          relative_distance(merge_id+1)=1
c$$$          relative_distance(merge_id+2)=0

        end subroutine ini_relative_distance


        subroutine ini_alignment(
     $     alignment,
     $     mainlayer_id,
     $     merge_id,
     $     relative_sizes,
     $     relative_distance)
    
          implicit none

          integer(ikind), dimension(2,2), intent(out) :: alignment
          integer                       , intent(in)  :: mainlayer_id
          integer                       , intent(in)  :: merge_id
          integer, dimension(:,:)       , intent(in)  :: relative_sizes
          integer, dimension(:)         , intent(in)  :: relative_distance


          select case(mainlayer_id)
            case(N,S)
               alignment(1,1) = bc_size+1
               alignment(1,2) = alignment(1,1)+
     $                          relative_distance(1)+
     $                          relative_sizes(1,1)
               do i=2, merge_id+1
                  alignment(1,2) =
     $                 alignment(1,2) +
     $                 2*bc_size+1 +
     $                 relative_sizes(1,i) +
     $                 relative_distance(i)
               end do
               alignment(1,2) =
     $              alignment(1,2)+
     $              relative_distance(merge_id+2)

            case(E,W)
               alignment(2,1) = bc_size+1
               alignment(2,2) = alignment(2,1)+
     $                          relative_distance(1)+
     $                          relative_sizes(1,1)
               do i=2, merge_id+1
                  alignment(2,2) =
     $                 alignment(2,2) +
     $                 2*bc_size+1 +
     $                 relative_sizes(1,i) +
     $                 relative_distance(i)
               end do
               alignment(2,2) =
     $              alignment(2,2)+
     $              relative_distance(merge_id+2)

            case default
               print '(''test_bf_sublayers_merge_prog'')'
               print '(''ini_alignment'')'
               print '(''mainlayer_id not recognized'')'
               print '(''mainlayer_id: '',I2)', mainlayer_id
               stop 'not implemented yet'
          end select

        end subroutine ini_alignment

        subroutine ini_neighbors(
     $     neighbor_case,
     $     neighbor_log1, neighbor_log2)

          implicit none

          integer, intent(in)  :: neighbor_case
          logical, intent(out) :: neighbor_log1
          logical, intent(out) :: neighbor_log2


          select case(neighbor_case)
            case(1)
               neighbor_log1=.false.
               neighbor_log2=.false.
            case(2)
               neighbor_log1=.true.
               neighbor_log2=.false.
            case(3)
               neighbor_log1=.false.
               neighbor_log2=.true.
            case(4)
               neighbor_log1=.true.
               neighbor_log2=.true.
          end select

        end subroutine ini_neighbors

      end program test_bf_sublayers_merge_prog
