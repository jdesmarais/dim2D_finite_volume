      !module for the fast initialization 
      !of sublayers for the test
      module sublayers_ini_module
      
        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : bc_size
        use parameters_kind    , only : ikind


        implicit none

        private
        public  :: ini_relative_sizes,
     $             ini_relative_distance,
     $             ini_alignment,
     $             ini_neighbors

        contains


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

          integer :: random_value,i

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

          integer :: i          

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

      end module sublayers_ini_module
