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
          procedure, nopass :: get_bc_section
          procedure, nopass :: analyse_grdpt_with_bc_section

          !only for tests
          procedure,   pass :: get_nb_ele_temp
          procedure,   pass :: get_nb_ele_final
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


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NE_edge_type)
               compatible = 
     $              ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).or.
     $              ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1)))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SE_edge_type)
               compatible = 
     $              ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $              ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1)))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SW_edge_type)
               compatible = 
     $              ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3)))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               end if


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NW_edge_type)
               compatible = 
     $              ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).or.
     $              ((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1)))

               if(compatible) then
                  bc_section(5) = bc_section(5)+1
                  remove_ele = bc_section(5).eq.3
               end if


          end select

        end function analyse_grdpt_with_bc_section


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


        subroutine print_bc_procedure(k,bc_section)

          implicit none
          
          integer              , intent(in) :: k
          integer, dimension(5), intent(in) :: bc_section

          select case(bc_section(1))
            case(N_edge_type)
               print '(I2,'' N_edge    i=['', 2I2, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(S_edge_type)
               print '(I2,'' S_edge    i=['', 2I2, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(E_edge_type)  
               print '(I2,'' E_edge    i=['', 2I2, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(W_edge_type)
               print '(I2,'' W_edge    i=['', 2I2, ''] j: '',I2)',
     $              k, bc_section(2:4)

            case(NE_corner_type)
               print '(I2,'' NE_corner [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(NW_corner_type)
               print '(I2,'' NW_corner [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(SE_corner_type)
               print '(I2,'' SE_corner [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(SW_corner_type)
               print '(I2,'' SW_corner [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(NE_edge_type)
               print '(I2,'' NE_edge   [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(NW_edge_type)
               print '(I2,'' NW_edge   [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(SE_edge_type) 
               print '(I2,'' SE_edge   [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

            case(SW_edge_type)
               print '(I2,'' SW_edge   [i_min,j_min]: '', 2I2)',
     $              k, bc_section(2:3)

          end select

        end subroutine print_bc_procedure

      end module bf_layer_bc_sections_class
