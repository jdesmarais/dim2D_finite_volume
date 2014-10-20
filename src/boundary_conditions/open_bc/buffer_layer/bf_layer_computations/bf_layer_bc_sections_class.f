      module bf_layer_bc_sections_class

        use parameters_bf_layer, only :
     $     bc_interior_pt,
     $     bc_pt

        use bf_layer_bc_procedure_module, only : 
     $     N_edge_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     NE_corner_type,
     $     NW_corner_type,
     $     SE_corner_type,
     $     SW_corner_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     SE_edge_type,
     $     SW_edge_type,
     $     get_bc_interior_pt_procedure

        implicit none

        private
        public :: bf_layer_bc_sections


        integer, parameter :: max_bc_sections_temp = 6


        type :: bf_layer_bc_sections

          integer                              :: nb_ele_temp
          integer, dimension(:,:), allocatable :: bc_sections_temp
          integer, dimension(:,:), allocatable :: bc_sections_buffer

          integer                              :: nb_ele_final
          integer, dimension(:,:), allocatable :: bc_sections_final

          contains

          procedure,   pass :: ini
          procedure,   pass :: deallocate_tables
          procedure,   pass :: add_to_temporary_bc_sections
          procedure,   pass :: add_to_final_bc_sections
          procedure,   pass :: add_to_bc_sections
          procedure,   pass :: remove_from_bc_sections_temp
          procedure,   pass :: remove_from_bc_sections_buffer
          procedure, nopass :: get_bc_section
          procedure, nopass :: analyse_grdpt_with_bc_section
          procedure,   pass :: analyse_grdpt
          procedure,   pass :: sort_bc_sections

          !only for tests
          procedure,   pass :: get_nb_ele_temp
          procedure,   pass :: get_nb_ele_final
          procedure,   pass :: get_bc_sections_temp
          procedure,   pass :: get_bc_sections_buffer
          procedure,   pass :: get_bc_sections_final
          procedure,   pass :: print_bc_sections

        end type bf_layer_bc_sections


        contains


        !initialize the number of elements in the object
        subroutine ini(this)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          
          this%nb_ele_temp=0
          this%nb_ele_final=0

        end subroutine ini


        !deallocate the allocatable atributes
        subroutine deallocate_tables(this)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this


          if(allocated(this%bc_sections_temp)) then
             deallocate(this%bc_sections_temp)
          end if

          if(allocated(this%bc_sections_buffer)) then
             deallocate(this%bc_sections_buffer)
          end if

          if(allocated(this%bc_sections_final)) then
             deallocate(this%bc_sections_final)
          end if

        end subroutine deallocate_tables


        !add a boundary section to list of boundary layers
        !bc_section, integer, dimension(5)
        ![procedure_type,edge_min,edge_max,coord,match_nb]
        subroutine add_to_temporary_bc_sections(this,bc_section)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer, dimension(5)      , intent(in)    :: bc_section

          integer, dimension(:,:), allocatable :: temp


          !increment the total number of elements saved in 
          this%nb_ele_temp=this%nb_ele_temp+1

          
          if(.not.allocated(this%bc_sections_temp)) then

             !first save in this%bc_sections_temp
             allocate(this%bc_sections_temp(5,max_bc_sections_temp))
             this%bc_sections_temp(:,1) = bc_section

          else

             !save in this%bc_sections_buffer
             if(this%nb_ele_temp.gt.max_bc_sections_temp) then
                
                if(.not.allocated(this%bc_sections_buffer)) then
                   allocate(this%bc_sections_buffer(5,1))
                   this%bc_sections_buffer(:,1) = bc_section

                else
                   if(this%nb_ele_temp.gt.(max_bc_sections_temp+size(this%bc_sections_buffer,2))) then

                      allocate(temp(5,size(this%bc_sections_buffer,2)+1))
                      temp(:,1:size(this%bc_sections_buffer,2)) = this%bc_sections_buffer(:,:)
                      temp(:,size(temp,2)) = bc_section(:)
                      call MOVE_ALLOC(temp,this%bc_sections_buffer)
                      
                   else
                      this%bc_sections_buffer(:,this%nb_ele_temp-max_bc_sections_temp) = bc_section
                   end if

                end if

             !save in this%bc_sections_temp
             else
                this%bc_sections_temp(:,this%nb_ele_temp) = bc_section
             end if

          end if

        end subroutine add_to_temporary_bc_sections


        !add the boundary section to the list of final boundary
        !sections
        subroutine add_to_final_bc_sections(this,bc_section)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer, dimension(5)      , intent(in)    :: bc_section

          integer, dimension(:,:), allocatable :: temp

          this%nb_ele_final = this%nb_ele_final+1

          if(.not.allocated(this%bc_sections_final)) then
             allocate(this%bc_sections_final(5,max_bc_sections_temp))
             this%bc_sections_final(:,1) = bc_section

          else
             if(this%nb_ele_final.gt.size(this%bc_sections_final,2)) then
                allocate(temp(5,size(this%bc_sections_final,2)+1))
                temp(:,1:size(temp,2)-1) = this%bc_sections_final(:,:)
                temp(:,size(temp,2))     = bc_section
                call MOVE_ALLOC(temp,this%bc_sections_final)

             else
                this%bc_sections_final(:,this%nb_ele_final) = bc_section

             end if

          end if

        end subroutine add_to_final_bc_sections


        !add the boundary section to the list of final boundary
        !sections
        subroutine add_to_bc_sections(this,bc_section)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer, dimension(5)      , intent(in)    :: bc_section

          select case(bc_section(1))

            !corner type boundary procedures are directly added
            !to the final bc_sections while the other are added
            !to the temporary boundary procedures
            case(
     $         NE_corner_type,
     $         NW_corner_type,
     $         SE_corner_type,
     $         SW_corner_type)

              call add_to_final_bc_sections(this,bc_section)

            case default

              call add_to_temporary_bc_sections(this,bc_section)

          end select

        end subroutine add_to_bc_sections


        subroutine remove_from_bc_sections_temp(this,k)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer                    , intent(in)    :: k

          integer :: l

          !before:           ------------------------
          !bf_sections_temp |  |  |  | k |  |  |  |  |
          !                  ------------------------
          !                  \______/    \__________/
          !                     |           __|
          !                    \|/        \|/
          !                   _ .___  _____.________
          !                  /      \/              \
          !after:            ------------------------
          !bf_sections_temp |  |  |  |  |  |  |  |  |
          !                  ------------------------
          !
          do l=k+1, min(this%nb_ele_temp,size(this%bc_sections_temp,2))
             this%bc_sections_temp(:,l-1) = 
     $            this%bc_sections_temp(:,l)
          end do

          if(this%nb_ele_temp.gt.size(this%bc_sections_temp,2)) then
             this%bc_sections_temp(:,size(this%bc_sections_temp,2)) =
     $            this%bc_sections_buffer(:,1)
             
             do l=2,this%nb_ele_temp-size(this%bc_sections_temp,2)
                this%bc_sections_buffer(:,l-1) =
     $               this%bc_sections_buffer(:,l)
             end do

          end if

          this%nb_ele_temp = this%nb_ele_temp-1

        end subroutine remove_from_bc_sections_temp


        subroutine remove_from_bc_sections_buffer(this,k)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer                    , intent(in)    :: k

          integer :: l

          !before:             ------------------------
          !bf_sections_buffer |  |  |  | k |  |  |  |  |
          !                    ------------------------
          !                    \______/    \__________/
          !                       |           __|
          !                      \|/        \|/
          !                     _ .___  _____.________
          !                    /      \/              \
          !after:              ------------------------
          !bf_sections_buffer |  |  |  |  |  |  |  |  |
          !                    ------------------------
          !
          do l=k+1,
     $         min(this%nb_ele_temp,
     $             this%nb_ele_temp-size(this%bc_sections_temp,2))
             this%bc_sections_buffer(:,l-1) =
     $            this%bc_sections_buffer(:,l)
          end do

          this%nb_ele_temp = this%nb_ele_temp-1

        end subroutine remove_from_bc_sections_buffer


        !using the boundary procedure given by bf_layer_bc_prcoedure
        !one can get the bc_section corresponding to the grid point(i,j)
        function get_bc_section(i,j,grdpts_id) result(bc_section)

          implicit none

          integer                , intent(in) :: i
          integer                , intent(in) :: j
          integer, dimension(:,:), intent(in) :: grdpts_id
          integer, dimension(5)               :: bc_section

          integer :: procedure_type
          integer :: i_proc
          integer :: j_proc
                    

          call get_bc_interior_pt_procedure(
     $         i,j,
     $         grdpts_id,
     $         procedure_type,
     $         i_proc,
     $         j_proc)

          bc_section(1)=procedure_type

          select case(procedure_type)
            case(N_edge_type,S_edge_type)
               bc_section(2) = i_proc
               bc_section(3) = i_proc
               bc_section(4) = j_proc

            case(E_edge_type,W_edge_type)
               bc_section(2) = j_proc
               bc_section(3) = j_proc
               bc_section(4) = i_proc

            case(SE_edge_type,SW_edge_type,NE_edge_type,NW_edge_type)
               bc_section(2)=i_proc
               bc_section(3)=j_proc
               bc_section(5)=1

            case(SE_corner_type,SW_corner_type,NE_corner_type,NW_corner_type)
               bc_section(2)=i_proc
               bc_section(3)=j_proc

            case default
               print '(''bf_layer_bc_sections'')'
               print '(''get_bc_section'')'
               print '(''procedure type: '',I2)', procedure_type
               stop 'procedure type not recognized'
               
          end select

        end function get_bc_section


        !check if the gridpoint tested (i,j) is comaptible with an existing
        !boundary layer bc_section
        !teh compatibility depends on the type of boundary layer
        function analyse_grdpt_with_bc_section(
     $     i,j,grdpts_id,bc_section,remove_ele)
     $     result(compatible)

          implicit none

          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer, dimension(5)  , intent(inout) :: bc_section
          logical                , intent(out)   :: remove_ele
          logical                                :: compatible

          remove_ele = .false.


          !type of procedure
          select case(bc_section(1))

            !bc_section(2):i_min
            !bc_section(3):i_max
            !bc_section(4):j
            case(N_edge_type)
               compatible =
     $              ((j-bc_section(4)).eq.0).and.
     $              ((i-bc_section(3)).eq.1).and.
     $              (grdpts_id(i,j+1).eq.bc_pt).and.
     $              (grdpts_id(i+1,j+1).eq.bc_pt).and.
     $              (grdpts_id(i+1,j).eq.bc_interior_pt)

               if(compatible) then
                  bc_section(3)=i
               end if

               remove_ele = j.gt.bc_section(4)
               

            !bc_section(2):i_min
            !bc_section(3):i_max
            !bc_section(4):j
            case(S_edge_type)
               compatible = 
     $              ((j-bc_section(4)).eq.0).and.
     $              ((i-bc_section(3)).eq.1).and.
     $              (grdpts_id(i,j-1).eq.bc_pt).and.
     $              (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $              (grdpts_id(i+1,j).eq.bc_interior_pt)

               if(compatible) then
                  bc_section(3)=i
               end if

               remove_ele = j.gt.bc_section(4)

            !bc_section(2): j_min
            !bc_section(3): j_max
            !bc_section(4): i
            case(E_edge_type)
               compatible =
     $              ((i-bc_section(4)).eq.0).and.
     $              ((j-bc_section(3)).eq.1).and.
     $              (grdpts_id(i+1,j).eq.bc_pt).and.
     $              (grdpts_id(i,j+1).eq.bc_interior_pt).and.
     $              (grdpts_id(i+1,j+1).eq.bc_pt)

               if(compatible) then
                  bc_section(3)=j
               end if

               remove_ele=j.gt.(bc_section(3)+1)


            !bc_section(2): j_min
            !bc_section(3): j_max
            !bc_section(4): i
            case(W_edge_type)
               compatible = 
     $              ((i-bc_section(4)).eq.0).and.
     $              ((j-bc_section(3)).eq.1).and.
     $              (grdpts_id(i-1,j).eq.bc_pt).and.
     $              (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $              (grdpts_id(i,j+1).eq.bc_interior_pt)

               if(compatible) then
                  bc_section(3)=j
               end if

               remove_ele=j.gt.(bc_section(3)+1)


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NE_edge_type)
               compatible = 
     $              ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $              (((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).and.(grdpts_id(i+1,j).eq.bc_interior_pt)).or.
     $              (((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i,j+1).eq.bc_interior_pt))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               else
                  remove_ele = j.gt.(bc_section(3)+1)
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SE_edge_type)
               compatible = 
     $              (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(grdpts_id(i,j-1).eq.bc_interior_pt)).or.
     $              ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).or.
     $              (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i+1,j).eq.bc_interior_pt))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               else
                  remove_ele = j.gt.(bc_section(3)+1)
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SW_edge_type)
               compatible = 
     $              (((i.eq.bc_section(2)+1).and.(j.eq.(bc_section(3)))).and.(grdpts_id(i,j-1).eq.bc_interior_pt)).or.
     $              (((i.eq.(bc_section(2))).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i-1,j).eq.bc_interior_pt)).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3)+1))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               else
                  remove_ele = j.gt.(bc_section(3)+1)
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NW_edge_type)
               compatible = 
     $              (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(grdpts_id(i-1,j).eq.bc_interior_pt)).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).or.
     $              (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i,j+1).eq.bc_interior_pt))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               else
                  remove_ele = j.gt.(bc_section(3)+1)
               end if


          end select

        end function analyse_grdpt_with_bc_section


        !analyse the grid point and decide whether it is part of an existing
        !boundary layer or whether it is the starting point of another boundary
        !layer
        subroutine analyse_grdpt(this,i,j,grdpts_id)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer                    , intent(in)    :: i
          integer                    , intent(in)    :: j
          integer, dimension(:,:)    , intent(in)    :: grdpts_id

          integer               :: k,k_buffer
          integer, dimension(5) :: new_bc_section
          logical               :: compatible
          logical               :: remove_ele


          !if there are already boundary layers, check whether the
          !grid point analysed is compatible with one of them
          if(this%nb_ele_temp.gt.0) then

             compatible = .false.

             !is the grid point compatible with the boundary layers
             !saved in this%bc_sections_temp
             k=1
             do while (k.le.min(this%nb_ele_temp,size(this%bc_sections_temp,2)))

                compatible = this%analyse_grdpt_with_bc_section(
     $               i,j,grdpts_id,
     $               this%bc_sections_temp(:,k),
     $               remove_ele)

                if(remove_ele) then
                   call this%add_to_final_bc_sections(this%bc_sections_temp(:,k))
                   call this%remove_from_bc_sections_temp(k)
                   k = k-1
                end if

                if(compatible) then
                   exit
                end if

                k=k+1

             end do

             if(.not.compatible) then

               !is the grid point compatible with the boundary layers
               !saved in this%bc_sections_buffer
                k=size(this%bc_sections_temp,2)+1
                do while(k.le.this%nb_ele_temp)
                   
                   k_buffer = k-size(this%bc_sections_temp,2)

                   compatible = this%analyse_grdpt_with_bc_section(
     $                  i,j,grdpts_id,
     $                  this%bc_sections_buffer(:,k_buffer),
     $                  remove_ele)

                   if(remove_ele) then
                      call this%add_to_final_bc_sections(this%bc_sections_buffer(:,k_buffer))
                      call this%remove_from_bc_sections_buffer(k)
                      k = k-1
                   end if

                   if(compatible) then
                      exit
                   end if

                   k=k+1

                end do


                !if no boundary layer matches the current grid point
                !the grid point is used as the starting point of a new
                !boundary layer
                if(.not.compatible) then

                   new_bc_section = this%get_bc_section(i,j,grdpts_id)
                   call this%add_to_bc_sections(new_bc_section)

                end if

             end if

          !if there is no existing boundary layers, the grid point
          !is analysed to know what type of starting point for a
          !boundary layer it is and then added to the boundary
          !layers of the object
          else

             new_bc_section = this%get_bc_section(i,j,grdpts_id)
             call this%add_to_bc_sections(new_bc_section)

          end if

        end subroutine analyse_grdpt


        subroutine sort_bc_sections(this,sorted_bc_sections)

          implicit none

          class(bf_layer_bc_sections)        , intent(in)  :: this
          integer, dimension(:,:),allocatable, intent(out) :: sorted_bc_sections

          integer :: nb_ele_tot
          integer :: k
          integer :: k_buffer
          integer :: k_final


          !1) create an intermediate table with the i_min,j_min of
          !   each boundary layer and the extent
          !   sorted_ele = [procedure_type,i_min,j_min,extent]
          !   extent = i_max for N_edge and S_edge
          !   extent = j_max for E_edge and W_edge
          !   extent = corner for other
          nb_ele_tot = this%nb_ele_temp + this%nb_ele_final
          allocate(sorted_bc_sections(4,nb_ele_tot))
          
          if(allocated(this%bc_sections_temp)) then
             do k=1, min(this%nb_ele_temp,size(this%bc_sections_temp,2))
                sorted_bc_sections(:,k) = get_sorted_ele(this%bc_sections_temp(:,k))
             end do
          end if

          if(allocated(this%bc_sections_buffer)) then
             do k=size(this%bc_sections_temp,2)+1,this%nb_ele_temp
                k_buffer = k-size(this%bc_sections_temp,2)
                sorted_bc_sections(:,k) = get_sorted_ele(this%bc_sections_buffer(:,k_buffer))
             end do
          end if

          if(allocated(this%bc_sections_final)) then
             do k=this%nb_ele_temp+1, this%nb_ele_temp+this%nb_ele_final
                k_final = k-this%nb_ele_temp
                sorted_bc_sections(:,k) = get_sorted_ele(this%bc_sections_final(:,k_final))
             end do
          end if


          !2) sort the table with increasing i and j
          call bubble_sort(sorted_bc_sections)

        end subroutine sort_bc_sections


        function order_bc_sections(p,q)
     $     result(p_larger_than_q)
        
          implicit none

          integer, dimension(4), intent(in) :: p
          integer, dimension(4), intent(in) :: q
          logical                           :: p_larger_than_q

          p_larger_than_q = .false.

          !if(p(j_min)>q(j_min))
          if (p(3)>q(3)) then
             p_larger_than_q = .true.
          else
             if(p(3).eq.q(3)) then

                !if(p(i_min)>q(i_min))
                if(p(2)>q(2)) then
                   p_larger_than_q = .true.
                end if

             end if
          end if

        end function order_bc_sections


        subroutine exchange_bc_sections(p,q)

          implicit none

          integer, dimension(4), intent(inout) :: p
          integer, dimension(4), intent(inout) :: q
          integer, dimension(4)                :: temp

          temp = p
          p    = q
          q    = temp

        end subroutine exchange_bc_sections

        
        subroutine bubble_sort(a)

          implicit none

          integer, dimension(:,:), intent(inout) :: a

          integer :: i,j
          integer :: n
          integer :: max_j

          n = size(a,2)

          do i=1, n

             max_j = 1

             do j=2, n-i+1
                if(order_bc_sections(a(:,j),a(:,max_j))) then
                   max_j = j
                end if
             end do

             call exchange_bc_sections(a(:,max_j),a(:,n-i+1))

          end do

        end subroutine bubble_sort


        function get_sorted_ele(bc_section) result(sorted_ele)

          implicit none

          integer, dimension(5), intent(in) :: bc_section
          integer, dimension(4)             :: sorted_ele

          select case(bc_section(1))

            case(N_edge_type)
               sorted_ele = [bc_section(1),bc_section(2),bc_section(4),bc_section(3)]

            case(S_edge_type)
               sorted_ele = [bc_section(1),bc_section(2),bc_section(4)-1,bc_section(3)]

            case(E_edge_type)
               sorted_ele = [bc_section(1),bc_section(4),bc_section(2),bc_section(3)]

            case(W_edge_type)
               sorted_ele = [bc_section(1),bc_section(4)-1,bc_section(2),bc_section(3)]

            case(
     $              NW_corner_type,
     $              NE_corner_type,
     $              SE_corner_type,
     $              SW_corner_type,
     $              NE_edge_type,
     $              NW_edge_type,
     $              SE_edge_type,
     $              SW_edge_type)
               sorted_ele = [bc_section(1),bc_section(2),bc_section(3),0]

            case default
               print '(''bf_layer_bc_sections_class'')'
               print '(''get_sorted_ele'')'
               print '(''case not recognized: '',I2)', bc_section(1)
               stop

          end select

        end function get_sorted_ele


        function get_nb_ele_temp(this) result(nb_ele_temp)

          implicit none

          class(bf_layer_bc_sections), intent(in) :: this
          integer                                 :: nb_ele_temp

          nb_ele_temp = this%nb_ele_temp

        end function get_nb_ele_temp


        function get_nb_ele_final(this) result(nb_ele_final)

          implicit none

          class(bf_layer_bc_sections), intent(in) :: this
          integer                                 :: nb_ele_final

          nb_ele_final = this%nb_ele_final

        end function get_nb_ele_final


        !print the boundary layers saved in the bc_sections
        subroutine print_bc_sections(this)
        
          implicit none

          class(bf_layer_bc_sections), intent(in) :: this

          integer :: k


          !print the bc_sections saved in this%bc_sections_temp
          !and this%bc_sections_buffer
          print '()'
          print '(''temporary bc_sections'')'
          do k=1, min(size(this%bc_sections_temp,2),this%nb_ele_temp)
             call print_bc_procedure(
     $            k,
     $            this%bc_sections_temp(:,k))
          end do

          do k=1, min(size(this%bc_sections_buffer,2),
     $         this%nb_ele_temp-max_bc_sections_temp)
             call print_bc_procedure(
     $            max_bc_sections_temp+k,
     $            this%bc_sections_buffer(:,k))
          end do
          print '()'


          !print the bc_sections saved in this%bc_sections_final
          print '()'
          print '(''final bc_sections'')'
          do k=1, min(this%nb_ele_final,size(this%bc_sections_final,2))
             call print_bc_procedure(
     $            k,
     $            this%bc_sections_final(:,k))
          end do
          print '()'

        end subroutine print_bc_sections


        subroutine get_bc_sections_temp(this,bc_sections_temp)

          implicit none

          class(bf_layer_bc_sections)         , intent(in)  :: this
          integer, dimension(:,:), allocatable, intent(out) :: bc_sections_temp

          integer :: k

          if(allocated(this%bc_sections_temp)) then
             allocate(bc_sections_temp(5,size(this%bc_sections_temp,2)))
             do k=1, size(this%bc_sections_temp,2)
                bc_sections_temp(:,k) = this%bc_sections_temp(:,k) 
             end do
          end if

        end subroutine get_bc_sections_temp


        subroutine get_bc_sections_buffer(this,bc_sections_buffer)

          implicit none

          class(bf_layer_bc_sections)         , intent(in)  :: this
          integer, dimension(:,:), allocatable, intent(out) :: bc_sections_buffer

          integer :: k

          if(allocated(this%bc_sections_buffer)) then
             allocate(bc_sections_buffer(5,size(this%bc_sections_buffer,2)))
             do k=1, size(this%bc_sections_buffer,2)
                bc_sections_buffer(:,k) = this%bc_sections_buffer(:,k) 
             end do
          end if

        end subroutine get_bc_sections_buffer


        subroutine get_bc_sections_final(this,bc_sections_final)

          implicit none

          class(bf_layer_bc_sections)         , intent(in)  :: this
          integer, dimension(:,:), allocatable, intent(out) :: bc_sections_final

          integer :: k

          if(allocated(this%bc_sections_final)) then
             allocate(bc_sections_final(5,size(this%bc_sections_final,2)))
             do k=1, size(this%bc_sections_final,2)
                bc_sections_final(:,k) = this%bc_sections_final(:,k) 
             end do
          end if

        end subroutine get_bc_sections_final


        subroutine print_bc_procedure(k,bc_section)

          implicit none
          
          integer              , intent(in) :: k
          integer, dimension(5), intent(in) :: bc_section

          select case(bc_section(1))
            case(N_edge_type)
               print '(I2,'' N_edge    i=['', 2I3, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(S_edge_type)
               print '(I2,'' S_edge    i=['', 2I3, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(E_edge_type)  
               print '(I2,'' E_edge    i=['', 2I3, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(W_edge_type)
               print '(I2,'' W_edge    i=['', 2I3, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(NE_corner_type)
               print '(I2,'' NE_corner [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(NW_corner_type)
               print '(I2,'' NW_corner [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(SE_corner_type)
               print '(I2,'' SE_corner [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(SW_corner_type)
               print '(I2,'' SW_corner [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(NE_edge_type)
               print '(I2,'' NE_edge   [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(NW_edge_type)
               print '(I2,'' NW_edge   [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(SE_edge_type) 
               print '(I2,'' SE_edge   [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

            case(SW_edge_type)
               print '(I2,'' SW_edge   [i_min,j_min]: '', 2I3)',
     $              k, bc_section(2:3)

          end select

        end subroutine print_bc_procedure

      end module bf_layer_bc_sections_class
