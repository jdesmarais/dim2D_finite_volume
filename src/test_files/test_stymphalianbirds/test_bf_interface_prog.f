      !> @file
      !> test file for the time integration procedures
      !> on the object 'bf_interface'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the computation of the time derivatives
      !> using the finite volume method by comparing 
      !> with expected data
      !
      !> @date
      ! 17_07_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_bf_interface_prog

        use bf_interface_class , only : bf_interface
        use bf_mainlayer_class , only : bf_mainlayer
        use bf_sublayer_class  , only : bf_sublayer
                               
        use parameters_bf_layer, only : interior_pt
                               
        use parameters_constant, only : periodic_xy_choice,
     $                                  bc_nodes_choice,
     $                                  N,S,E,W
        use parameters_input   , only : nx,ny,ne,
     $                                  bc_choice,bc_size,
     $                                  bcx_type_choice,
     $                                  bcy_type_choice
        use parameters_kind    , only : ikind, rkind

        implicit none

        
        !<operators tested
        integer    , parameter                     :: nxt = 10
        integer    , parameter                     :: nyt = 6
        real(rkind), dimension(nx,ny,ne)           :: interior_nodes
        real(rkind), dimension(:,:,:), allocatable :: nodes
        integer    , dimension(:,:)  , allocatable :: grdpts_id
        real(rkind)                                :: dx
        real(rkind)                                :: dy

        type(bf_interface)                         :: bf_interface_used
        type(bf_mainlayer), pointer                :: bf_mainlayer_used
        integer    , parameter                     :: nb_sublayers=3
        type(bf_sublayer), pointer                 :: added_sublayer
        integer(ikind), dimension(2,2)             :: alignment

        type(bf_sublayer), pointer                 :: current_sublayer
        real(rkind), dimension(:,:,:), allocatable :: time_dev


        !<CPU recorded times
        real :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.true.
        integer(ikind)             :: i,j
        logical                    :: test_parameter


        !<if nx<4, ny<4 then the test cannot be done
        test_parameter=.true.
        test_parameter=test_parameter.and.(nxt.eq.10)
        test_parameter=test_parameter.and.(nyt.eq.6)
        test_parameter=test_parameter.and.(ne.eq.1)
        test_parameter=test_parameter.and.(bc_choice.eq.periodic_xy_choice)
        test_parameter=test_parameter.and.(bcx_type_choice.eq.bc_nodes_choice)
        test_parameter=test_parameter.and.(bcy_type_choice.eq.bc_nodes_choice)    
        if(.not.test_parameter) then
           print *, 'the test requires several parameters'
           print *, 'test designed for simpletest eq'
           print *, 'nx=10'
           print *, 'ny=6'
           print *, 'ne=1'
           print *, 'bc_choice=periodic_xy_choice'
           print *, 'bcx_type_choice=bc_nodes_choice'
           print *, 'bcy_type_choice=bc_nodes_choice'
           stop ''
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        dx=1.0
        dy=1.0


        !initialization of the bf_interface
        call bf_interface_used%ini()

        do j=1,4
           do i=1, nb_sublayers

              alignment = get_alignment(i,j)
              
              added_sublayer => bf_interface_used%allocate_sublayer(
     $             j, interior_nodes, alignment, dx, dy)
              
              call initialize_nodes_and_grdpts(nodes,grdpts_id)
              
              call added_sublayer%set_nodes(nodes)
              call added_sublayer%set_grdpts_id(grdpts_id)

           end do
        end do


        !computation of the time derivatives + check
        do j=1,4

           print '(''mainlayer:'', I1)', j
           
           bf_mainlayer_used => bf_interface_used%get_mainlayer(j)

           current_sublayer => bf_mainlayer_used%get_head_sublayer()
        
           do i=1, nb_sublayers
              
              call current_sublayer%allocate_before_timeInt()
              call current_sublayer%compute_time_dev()
              call current_sublayer%get_time_dev(time_dev)
              call check_time_dev(detailled, time_dev)
              call current_sublayer%deallocate_after_timeInt()

              current_sublayer => current_sublayer%get_next()
              
           end do

           print '()'

        end do   

        call CPU_TIME(time2)
        print '(''time elapsed:'', F6.2)', time2-time1


        contains


        subroutine check_time_dev(detailled, time_dev)

          implicit none

          logical                      , intent(in) :: detailled
          real(rkind), dimension(:,:,:), intent(in) :: time_dev

          logical        :: global,local
          integer(ikind) :: i,j

          !<print the time derivatives
          global=.true.
          do j=1+bc_size, nyt-bc_size
             do i=1+bc_size, nxt-bc_size
                local=(time_dev(i,j,1).eq.(-20.0d0))
                global=global.and.local
                if(detailled) then
                   print *, i,j, local
                end if
             end do
          end do
          print '(''test_validated: '', L1)', global          

        end subroutine check_time_dev        


        subroutine initialize_nodes_and_grdpts(nodes,grdpts_id)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(out) :: nodes
          integer    , dimension(:,:)  , allocatable, intent(out) :: grdpts_id

          integer :: i,j

          allocate(nodes(nxt,nyt,ne))
          allocate(grdpts_id(nxt,nyt))
          
          do j=1, nyt
             do i=1, nxt
                nodes(i,j,1) = i + (j-1)*nxt
                grdpts_id(i,j) = interior_pt
             end do
          end do

        end subroutine initialize_nodes_and_grdpts


        function get_alignment(i,j) result(alignment)

          implicit none

          integer           , intent(in) :: i,j
          integer(ikind), dimension(2,2) :: alignment


          select case(j)
            case(N)
               alignment = get_alignment_north(i)
            case(S)
               alignment = get_alignment_south(i)
            case(E)
               alignment = get_alignment_east(i)
            case(W)
               alignment = get_alignment_west(i)
            case default
               alignment(1,1) = 0
               print '(''test_bf_interface'')'
               print '(''get_alignment'')'
               stop 'case not recognized'
          end select          

        end function get_alignment


        function get_alignment_north(i) result(alignment)

          implicit none

          integer           , intent(in) :: i
          integer(ikind), dimension(2,2) :: alignment

          select case(i)
            case(1)
               alignment(1,1) = bc_size+1
            case(2)
               alignment(1,1) = nxt+bc_size
            case(3)
               alignment(1,1) = 2*nxt+bc_size
            case default
               alignment(1,1) = 0
               print '(''test_bf_mainlayer_prog'')'
               print '(''get_alignment_north'')'
               stop 'case not recognized'
          end select

          alignment(2,1) = nx-bc_size+1
          alignment(1,2) = alignment(1,1)-2*bc_size+nxt-1
          alignment(2,2) = alignment(2,1)-2*bc_size+nyt-1

        end function get_alignment_north


        function get_alignment_south(i) result(alignment)

          implicit none

          integer           , intent(in) :: i
          integer(ikind), dimension(2,2) :: alignment

          select case(i)
            case(1)
               alignment(1,1) = bc_size+1
            case(2)
               alignment(1,1) = nxt+bc_size
            case(3)
               alignment(1,1) = 2*nxt+bc_size
            case default
               alignment(1,1) = 0
               alignment(2,2) = 0
               print '(''test_bf_mainlayer_prog'')'
               print '(''get_alignment_south'')'
               stop 'case not recognized'
          end select

          alignment(2,2) = bc_size
          alignment(1,2) = alignment(1,1)-2*bc_size+nxt-1
          alignment(2,1) = alignment(2,2)+2*bc_size-nyt+1

        end function get_alignment_south


        function get_alignment_west(i) result(alignment)

          implicit none

          integer           , intent(in) :: i
          integer(ikind), dimension(2,2) :: alignment

          select case(i)
            case(1)
               alignment(2,1) = bc_size+1
            case(2)
               alignment(2,1) = nyt+bc_size
            case(3)
               alignment(2,1) = bc_size
            case default
               alignment(2,1) = 0
               print '(''test_bf_mainlayer_prog'')'
               print '(''get_alignment_west'')'
               stop 'case not recognized'
          end select

          alignment(1,2) = bc_size
          alignment(1,1) = alignment(1,2)+2*bc_size-nxt+1
          alignment(2,2) = alignment(2,1)-2*bc_size+nyt-1

        end function get_alignment_west


        function get_alignment_east(i) result(alignment)

          implicit none

          integer           , intent(in) :: i
          integer(ikind), dimension(2,2) :: alignment

          select case(i)
            case(1)
               alignment(2,1) = bc_size+1
            case(2)
               alignment(2,1) = nyt+bc_size
            case(3)
               alignment(2,1) = bc_size
            case default
               alignment(2,1) = 0
               print '(''test_bf_mainlayer_prog'')'
               print '(''get_alignment_east'')'
               stop 'case not recognized'
          end select

          alignment(1,1) = nx-bc_size+1
          alignment(1,2) = alignment(1,1)-2*bc_size+nxt-1
          alignment(2,2) = alignment(2,1)-2*bc_size+nyt-1

        end function get_alignment_east

      end program test_bf_interface_prog
