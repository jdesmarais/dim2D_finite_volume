      !> @file
      !> module encapsulating the bf_layer_bc_sections object. It
      !> encapsulates the functions needed to analyze the grdpts_id
      !> and determine where the boundary layers are located. The
      !> procedures needed to compute the boundary points are then
      !> deduced
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the bf_layer_bc_sections object. It
      !> encapsulates the functions needed to analyze the grdpts_id
      !> and determine where the boundary layers are located. The
      !> procedures needed to compute the boundary points are then
      !> deduced
      !
      !> @date
      ! 26_01_2015 - documentation update  - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_bc_sections_class

        use parameters_bf_layer, only :
     $     bc_interior_pt,
     $     bc_pt,
     $     BF_SUCCESS

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

        use bf_layer_errors_module, only :
     $       error_overlap_index,
     $       error_overlap_incompatible

        implicit none

        private
        public ::
     $       bf_layer_bc_sections,
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap,
     $       determine_edge_points_computed


        integer, parameter :: max_bc_sections_temp = 6

        integer, parameter :: no_overlap = 0
        integer, parameter :: N_overlap  = 1
        integer, parameter :: S_overlap  = 2
        integer, parameter :: E_overlap  = 3
        integer, parameter :: W_overlap  = 4
        integer, parameter :: NE_overlap = 5
        integer, parameter :: NW_overlap = 6
        integer, parameter :: SE_overlap = 7
        integer, parameter :: SW_overlap = 8


        !> @class bf_layer_bc_sections
        !> class encapsulating the bf_layer_bc_sections object. It
        !> encapsulates the functions needed to analyze the grdpts_id
        !> and determine where the boundary layers are located. The
        !> procedures needed to compute the boundary points are then
        !> deduced. When analyzing the grdpts_id, to prevent numerous
        !> re-allocation, is used
        !
        !> @param nb_ele_temp
        !> number of bc_sections stored in the temporary
        !> bc_sections_temp attribute waiting to be completed by
        !> other grdpts_id before being finalized in bc_sections_final
        !
        !> @param bc_sections_temp
        !> the bc_setions waiting to be completed by other grdpts_id
        !> are stored in an array made of two array: bc_sections_temp
        !> and bc_sctions_buffer. The size of bc_sections_temp is
        !> fixed by max_bc_sections_temp whike bc_sections_buffer can
        !> be reallocated depending on the needs
        !
        !> @param bc_sections_buffer
        !> array that can be reallocated depending on the needs to
        !> store the temporary bc_sections
        !
        !> @param nb_ele_final
        !> number of elements stored in the bc_sections_final table
        !
        !> @param bc_sections_final
        !> array where the bc_scetions once completed are stored
        !
        !> @param ini
        !> initialize the main attributes of bf_bc_sections_class
        !
        !> @param deallocate_tables
        !> deallocate the array attributes
        !
        !> @param add_to_temporary_bc_sections
        !> add the bc_section to the bc_sections_temp or
        !> bc_sections_buffer attribute depending on the number of 
        !> temporary bc_sections already stored
        !
        !> @param add_to_final_bc_sections
        !> add the bc_section to the bc_sections_final array
        !
        !> @param add_to_bc_sections
        !> add the bc_sections iether to the temporary bc_sections
        !> or the final bc_sections
        !
        !> @param remove_from_bc_sections_temp
        !> remove a bc_section from the bc_sections_temp array
        !
        !> @param remove_from_bc_sections_buffer
        !> remove a bc_section from the bc_sections_buffer array
        !
        !> @param get_bc_section
        !> analyze the grid point ID and determine the corresponding
        !> bc_section (alone if there is no bc_section stoerd in the
        !> temporary array or with the existing bc_section)
        !
        !> @param analyse_grdpt_with_bc_section
        !> analyze the grid point ID and compare it to previous bc_section
        !> to know whether the bc_section can be completed
        !
        !> @param analyze_grdpt
        !> analyze the grid point ID and determine the corresponding
        !> bc_section w/o comparing it to the existing bc_sections
        !
        !> @param sort_bc_sections
        !> gather the bc_sections saved in final and temporary arrays
        !> and order them in increasing j and increasing i
        !
        !> @param add_overlap_between_corners_and_anti_corners
        !> analyze the ordered bc_sections and mark the overlapping
        !> gridpoints b/w the corner and anti-corner boundary layers
        !
        !> @param get_nb_ele_temp
        !> get the nb_ele_temp attribute
        !
        !> @param get_nb_ele_final
        !> get the nb_ele_final attribute
        !
        !> @param get_bc_sections_temp
        !> get the bc_sections_temp attribute
        !
        !> @param get_bc_sections_buffer
        !> get the bc_sections_buffer attribute
        !
        !> @param print_bc_sections
        !> display the bc_sections in a graphical form
        !------------------------------------------------------------
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
          procedure, nopass :: add_overlap_between_corners_and_anti_corners
          procedure,   pass :: finalize_bc_sections

          !only for tests
          procedure,   pass :: get_nb_ele_temp
          procedure,   pass :: get_nb_ele_final
          procedure,   pass :: get_bc_sections_temp
          procedure,   pass :: get_bc_sections_buffer
          procedure,   pass :: get_bc_sections_final
          procedure,   pass :: print_bc_sections

        end type bf_layer_bc_sections


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the number of elements in the object
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          
          this%nb_ele_temp=0
          this%nb_ele_final=0

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the allocatable atributes
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add a boundary section to list of boundary layers
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param bc_section
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the boundary section to the list of final boundary
        !> sections
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param bc_section
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the boundary section to the list of final boundary
        !> sections
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param bc_section
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the element (k) from the temporary bc_section
        !> array
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param k
        !> index identifying the element in the temporary bc_section
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the element (k) from the bc_sections_buffer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param k
        !> index identifying the element in the temporary bc_section
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> using the boundary procedure given by bf_layer_bc_procedure
        !> one can get the bc_section corresponding to the grid-
        !> point(i,j)
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param i
        !> integer identifying the x-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param j
        !> integer identifying the y-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param grdpts_id
        !> identity of the grid-point
        !
        !> @param ierror
        !> integer identifying whether the identification of the 
        !> bc_section was successful or not
        !--------------------------------------------------------------
        function get_bc_section(i,j,grdpts_id,ierror) result(bc_section)

          implicit none

          integer                , intent(in)  :: i
          integer                , intent(in)  :: j
          integer, dimension(:,:), intent(in)  :: grdpts_id
          logical                , intent(out) :: ierror
          integer, dimension(5)                :: bc_section

          integer :: procedure_type
          integer :: i_proc
          integer :: j_proc

          ierror = BF_SUCCESS                    

          call get_bc_interior_pt_procedure(
     $         i,j,
     $         grdpts_id,
     $         procedure_type,
     $         i_proc,
     $         j_proc,ierror)

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
               print '(''procedure not recognized'')'
               print '(''****************************************'')'
               print '()'
               ierror = .not.BF_SUCCESS
               
          end select

        end function get_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the gridpoint tested (i,j) is compatible with
        !> an existing boundary layer bc_section, the compatibility
        !> critrion depends on the type of boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param i
        !> integer identifying the x-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param j
        !> integer identifying the y-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param grdpts_id
        !> identity of the grid-point
        !
        !> @param bc_section
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !
        !> @param remove_ele
        !> integer identifying whether the bc_section analyzed is not
        !> complete and should be transfered to the bc_sections_final
        !> array
        !
        !> @return compatible
        !> logical identifying whether the grdpts_id(i,j) was compatible
        !> with the bc_section
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyse the grid point and decide whether it is part of
        !> an existing boundary layer or whether it is the starting
        !> point of another boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param i
        !> integer identifying the x-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param j
        !> integer identifying the y-coordinate of the grdpts_id
        !> analyzed
        !
        !> @param grdpts_id
        !> identity of the grid-point
        !
        !> @param ierror
        !> integer identifying whether the analyze of the grdpt was
        !> successful
        !--------------------------------------------------------------
        subroutine analyse_grdpt(this,i,j,grdpts_id,ierror)

          implicit none

          class(bf_layer_bc_sections), intent(inout) :: this
          integer                    , intent(in)    :: i
          integer                    , intent(in)    :: j
          integer, dimension(:,:)    , intent(in)    :: grdpts_id
          logical                    , intent(out)   :: ierror

          integer               :: k,k_buffer
          integer, dimension(5) :: new_bc_section
          logical               :: compatible
          logical               :: remove_ele


          !this condition prevents to test the grid points
          !at the very edge of the grdpts_id array while it
          !should not be possible since bc_interior_pt are
          !always one grid point away from the border
          !this is one exception to this rule: when the buffer
          !layer is attached to the interior and both borders are 
          !strictly inside the interior
          !
          !                          ____ buffer layer
          !         |---------------|
          !         | 3 3 3 3 3 3 3 |
          ! 3 3 3 3 | 3 2 2 2 2 2 3 | 3 3 3 3
          ! 2 2 2 2 | 2 2 1 1 1 2 2 | 2 2 2 2
          ! ------- | ------------- | -------
          !         | 1 1 1 1 1 1 1 |
          !         | 1 1 1 1 1 1 1 |
          !         -----------------  
          !          interior
          !
          !the left bc_interior_pt which is at the edge
          !belongs to an NW_edge. As this edge is treated
          !as a corner, the grid points outside the buffer
          !layer are not needed for the computation and as
          !the grid point next to it can used to determine
          !that it is a corner, the condition below is
          !justified
          !the only restriction is if the NW_edge should be 
          !treated in a way different from the corner
          if((i.ne.1).and.(i.ne.size(grdpts_id,1)).and.(j.ne.1).and.(j.ne.size(grdpts_id,2))) then
             

             !if there are already boundary layers, check whether the
             !grid point analysed is compatible with one of them
             if(this%nb_ele_temp.gt.0) then
             
                compatible = .false.
             
                !is the grid point compatible with the boundary layers
                !saved in this%bc_sections_temp
                k=1
                do while (k.le.min(this%nb_ele_temp,size(this%bc_sections_temp,2)))
             
                   compatible = this%analyse_grdpt_with_bc_section(
     $                  i,j,grdpts_id,
     $                  this%bc_sections_temp(:,k),
     $                  remove_ele)
             
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
     $                     i,j,grdpts_id,
     $                     this%bc_sections_buffer(:,k_buffer),
     $                     remove_ele)
             
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
             
                      new_bc_section = this%get_bc_section(i,j,grdpts_id,ierror)
                      call this%add_to_bc_sections(new_bc_section)
             
                   end if
             
                end if
             
             !if there is no existing boundary layers, the grid point
             !is analysed to know what type of starting point for a
             !boundary layer it is and then added to the boundary
             !layers of the object
             else
             
                new_bc_section = this%get_bc_section(i,j,grdpts_id,ierror)
                call this%add_to_bc_sections(new_bc_section)
             
             end if

          end if

        end subroutine analyse_grdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> finalize the bc_sections analyzed by the
        !> bf_layer_bc_sections object:
        !> 1) keep only the information needed to compute the
        !>    boundary grid-points
        !> 2) sort the bc_sections in increasing j and i to
        !>    improve cache efficieny when computing the
        !>    boundary layer gridpoints
        !> 3) determine the overlap b/w the corner and anti-corner
        !>    boundary procedures to prevent interactions and symetry
        !>    violation when applying the boundary procedures
        !> 4) deallocate the intermediate attributes used to analyze
        !>    the boundary layers
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param sorted_bc_sections
        !> array with the boundary sections sorted in increasing j
        !> and increasing i
        !--------------------------------------------------------------
        subroutine finalize_bc_sections(this,sorted_bc_sections)

          implicit none

          class(bf_layer_bc_sections)        , intent(inout) :: this
          integer, dimension(:,:),allocatable, intent(out)   :: sorted_bc_sections


          !> 1) keep only the information needed to compute the
          !>    boundary grid-points
          !> 2) sort the bc_sections in increasing j and i to
          !>    improve cache efficieny when computing the
          !>    boundary layer gridpoints
          call this%sort_bc_sections(sorted_bc_sections)


          !4) deallocate the intermediate attributes used to analyze
          !>  the boundary layers
          call this%deallocate_tables()


          !3) determine the overlap b/w the corner and anti-corner
          !>    boundary procedures to prevent interactions and symetry
          !>    violation when applying the boundary procedures
          call this%add_overlap_between_corners_and_anti_corners(
     $         sorted_bc_sections)
          
        end subroutine finalize_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> sort the boundary layers with increasing j and
        !> increasing i
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param sorted_bc_sections
        !> array with the boundary sections sorted in increasing j
        !> and increasing i
        !--------------------------------------------------------------
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
          if(nb_ele_tot.gt.0) then
             allocate(sorted_bc_sections(4,nb_ele_tot))
          end if
          
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> operator ordering two bc_sections
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param p
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !
        !> @param q
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !
        !> @return p_larger_than_q
        !> logical indicating whether p>q
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> permutation of two bc_sections
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param p
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !
        !> @param q
        !> representation of a boundary layer: integer, dimension(5)
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !--------------------------------------------------------------
        subroutine exchange_bc_sections(p,q)

          implicit none

          integer, dimension(4), intent(inout) :: p
          integer, dimension(4), intent(inout) :: q
          integer, dimension(4)                :: temp

          temp = p
          p    = q
          q    = temp

        end subroutine exchange_bc_sections

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> use bubble sorting to order the bc_section array
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param a
        !> bc_section array sorted
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> turn a boundary layer represented with 
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !> into a reduced element where only the information needed
        !> to compute the boundary points are kept:
        !> [procedure_type,i_min,j_min,extent]
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary layer represented as
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !
        !> @param sorted_ele
        !> boundary layer represented as
        !> [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
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

               sorted_ele = [
     $           bc_section(1),
     $           bc_section(2),
     $           bc_section(3),
     $           no_overlap]

            case default
               print '(''bf_layer_bc_sections_class'')'
               print '(''get_sorted_ele'')'
               print '(''case not recognized: '',I2)', bc_section(1)
               stop

          end select

        end function get_sorted_ele


        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyse the sorted elements of boundary layers to identify
        !> whether some boundary elements overlap (corners and
        !> anti-corners)
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_sections_sorted
        !> boundary layer represented as
        !> [procedure_type,i_min,j_min,extent]
        !> and ordered with increasing i and increasing j
        !--------------------------------------------------------------
        subroutine add_overlap_between_corners_and_anti_corners(
     $     bc_sections_sorted)

          implicit none

          integer, dimension(:,:), intent(inout) :: bc_sections_sorted

          integer :: k_prev_stage    !index where the bc_section for j-1 begins
          integer :: k_current_stage !index where the bc_section for j begins

          integer :: j_current_stage !index for the y-coordinate of the current stage
          integer :: j_stage         !index for the y-coordinate of the bc_section analyzed

          integer :: k               !index for the bc_section analyzed


          k_prev_stage    = 1
          k_current_stage = 1
          j_current_stage = get_j_stage(bc_sections_sorted(:,1))
          
          ! loop over the bc_sections stored in
          ! bc_sections_sorted
          do k=1, size(bc_sections_sorted,2)

             ! get the j_stage identifying the
             ! j_min component of the bc_section
             ! this way, we can determine how far
             ! in the bc_sections_sorted list we
             ! need to look for to get potential
             ! anti-corner bc_sections overlaped
             ! by the corner bc_sections
             j_stage = get_j_stage(bc_sections_sorted(:,k))
             
             ! update the indices identifying the
             ! prev and next j-stages
             if(j_stage.gt.j_current_stage) then
                k_prev_stage    = k_current_stage
                k_current_stage = k
                j_current_stage = j_stage
             end if

             ! if the bc_section analyzed is a corner,
             ! it should be compared to the bc_sections
             ! of the previous and the next stages
             ! corresponding to the stages where an
             ! overlap is possible
             if(is_a_corner(bc_sections_sorted(:,k))) then
                
                call compare_corner_to_previous_stage_bc_sections(
     $               bc_sections_sorted(:,k),
     $               bc_sections_sorted,
     $               k_prev_stage,
     $               k-1)

                call compare_corner_to_next_stage_bc_sections(
     $               bc_sections_sorted(:,k),
     $               bc_sections_sorted,
     $               k+1,
     $               j_stage)

             end if

          end do

        end subroutine add_overlap_between_corners_and_anti_corners


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compare the corner located at (i,j) with the boundary layer
        !> elements located b/w k_min an k_max in bc_sections_sorted
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param corner
        !> boundary layer represented as [corner_type,i_min,j_min,extent]
        !
        !> @param bc_sections_sorted
        !> boundary layers represented as
        !> [procedure_type,i_min,j_min,extent]
        !> and ordered with increasing i and increasing j
        !
        !> @param k_min
        !> index identifying the first element compared in
        !> bc_sections_sorted
        !
        !> @param k_max
        !> index identifying the last element compared in
        !> bc_sections_sorted
        !--------------------------------------------------------------
        subroutine compare_corner_to_previous_stage_bc_sections(
     $     corner,
     $     bc_sections_sorted,
     $     k_min,
     $     k_max)

          implicit none

          integer, dimension(4)  , intent(in)    :: corner
          integer, dimension(:,:), intent(inout) :: bc_sections_sorted
          integer                , intent(in)    :: k_min
          integer                , intent(in)    :: k_max

          integer :: k

          do k=k_min,k_max

             if(is_an_anti_corner(bc_sections_sorted(:,k))) then

                call overlap(corner,bc_sections_sorted(:,k))

             end if

          end do

        end subroutine compare_corner_to_previous_stage_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compare the corner located at (i,j) with the boundary layer
        !> elements located from k_min and whose y-position j is such
        !> that j.le.(j_stage+1) to allow an overlap b/w the corner
        !> element and the bc_section element
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param corner
        !> boundary layer represented as [corner_type,i_min,j_min,extent]
        !
        !> @param bc_sections_sorted
        !> boundary layers represented as
        !> [procedure_type,i_min,j_min,extent]
        !> and ordered with increasing i and increasing j
        !
        !> @param k_min
        !> index identifying the first element compared in
        !> bc_sections_sorted
        !
        !> @param j_stage
        !> index identifying the y-position of the corner
        !--------------------------------------------------------------
        subroutine compare_corner_to_next_stage_bc_sections(
     $     corner,
     $     bc_sections_sorted,
     $     k_min,
     $     j_stage)

          implicit none

          integer, dimension(4)  , intent(in)    :: corner
          integer, dimension(:,:), intent(inout) :: bc_sections_sorted
          integer                , intent(in)    :: k_min
          integer                , intent(in)    :: j_stage

          integer :: k
          integer :: j

          k=k_min
          j=j_stage

          do while((k.le.size(bc_sections_sorted,2)).and.(j.le.(j_stage+1)))

             j = get_j_stage(bc_sections_sorted(:,k))

             if(is_an_anti_corner(bc_sections_sorted(:,k))) then

                call overlap(corner,bc_sections_sorted(:,k))

             end if

             k=k+1

          end do

        end subroutine compare_corner_to_next_stage_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of the grid-points in common
        !> with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param corner
        !> boundary layer represented as [corner_type,i_min,j_min,extent]
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        subroutine overlap(corner,anti_corner)

          implicit none

          integer, dimension(4), intent(in)    :: corner
          integer, dimension(4), intent(inout) :: anti_corner


          if(corner(3).eq.anti_corner(3)) then

          !i_corner = i_anti_corner+1 AND j_corner = j_anti_corner
             if(corner(2).eq.(anti_corner(2)+1)) then
                call overlap_E(anti_corner)
             else
             
          !i_corner = i_anti_corner-1 AND j_corner = j_anti_corner
                if(corner(2).eq.(anti_corner(2)-1)) then
                   call overlap_W(anti_corner)
                end if
             end if
          end if


          if(corner(2).eq.anti_corner(2)) then

          !j_corner = j_anti_corner+1 AND i_corner = i_anti_corner
             if(corner(3).eq.(anti_corner(3)+1)) then
                call overlap_N(anti_corner)
             else
             
          !j_corner = j_anti_corner-1 AND i_corner = i_anti_corner
                if(corner(3).eq.(anti_corner(3)-1)) then
                   call overlap_S(anti_corner)
                end if
             end if
          end if

        end subroutine overlap


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of North the grid-points in
        !> common with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        subroutine overlap_N(anti_corner)

          implicit none

          integer, dimension(4), intent(inout) :: anti_corner

          
          select case(anti_corner(4))

            case(no_overlap)
               anti_corner(4) = N_overlap

            case(N_overlap,NE_overlap,NW_overlap)
               anti_corner(4) = anti_corner(4)

            case(S_overlap)
               call error_overlap_incompatible(
     $              'bf_layer_bc_sections_class',
     $              'overlap_N',
     $              anti_corner(4),
     $              N_overlap)

            case(E_overlap)
               anti_corner(4) = NE_overlap

            case(W_overlap)
               anti_corner(4) = NW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_N',
     $              anti_corner(4))

          end select

        end subroutine overlap_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of South the grid-points in
        !> common with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        subroutine overlap_S(anti_corner)

          implicit none

          integer, dimension(4), intent(inout) :: anti_corner

          
          select case(anti_corner(4))

            case(no_overlap)
               anti_corner(4) = S_overlap

            case(S_overlap,SE_overlap,SW_overlap)
               anti_corner(4) = anti_corner(4)

            case(N_overlap)
               call error_overlap_incompatible(
     $              'bf_layer_bc_sections_class',
     $              'overlap_S',
     $              anti_corner(4),
     $              S_overlap)

            case(E_overlap)
               anti_corner(4) = SE_overlap

            case(W_overlap)
               anti_corner(4) = SW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_S',
     $              anti_corner(4))

          end select

        end subroutine overlap_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of East the grid-points in
        !> common with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        subroutine overlap_E(anti_corner)

          implicit none

          integer, dimension(4), intent(inout) :: anti_corner

          
          select case(anti_corner(4))

            case(no_overlap)
               anti_corner(4) = E_overlap

            case(E_overlap,SE_overlap,NE_overlap)
               anti_corner(4) = anti_corner(4)
               
            case(W_overlap)
               call error_overlap_incompatible(
     $              'bf_layer_bc_sections_class',
     $              'overlap_E',
     $              anti_corner(4),
     $              E_overlap)

            case(N_overlap)
               anti_corner(4) = NE_overlap

            case(S_overlap)
               anti_corner(4) = SE_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_E',
     $              anti_corner(4))

          end select

        end subroutine overlap_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> modify the properties of the anti_corner boundary layer
        !> to remove the computation of West the grid-points in
        !> common with the corner boundary layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param anti-corner
        !> boundary layer represented as [anti_corner_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        subroutine overlap_W(anti_corner)

          implicit none

          integer, dimension(4), intent(inout) :: anti_corner

          
          select case(anti_corner(4))

            case(no_overlap)
               anti_corner(4) = W_overlap

            case(W_overlap,SW_overlap,NW_overlap)
               anti_corner(4) = anti_corner(4)
               
            case(E_overlap)
               call error_overlap_incompatible(
     $              'bf_layer_bc_sections_class',
     $              'overlap_W',
     $              anti_corner(4),
     $              W_overlap)

            case(N_overlap)
               anti_corner(4) = NW_overlap

            case(S_overlap)
               anti_corner(4) = SW_overlap

            case default
               call error_overlap_index(
     $              'bf_layer_bc_sections_class',
     $              'overlap_W',
     $              anti_corner(4))

          end select

        end subroutine overlap_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the boundary layer is an anti-corner
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary layer represented as [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        function is_an_anti_corner(bc_section)

          implicit none

          integer, dimension(4), intent(in) :: bc_section
          logical                           :: is_an_anti_corner

          integer :: procedure_type

          procedure_type = bc_section(1)

          is_an_anti_corner = (procedure_type.eq.NE_edge_type).or.
     $                        (procedure_type.eq.NW_edge_type).or.
     $                        (procedure_type.eq.SE_edge_type).or.
     $                        (procedure_type.eq.SW_edge_type)

        end function is_an_anti_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the boundary layer is a corner
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary layer represented as [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        function is_a_corner(bc_section)

          implicit none

          integer, dimension(4), intent(in) :: bc_section
          logical                           :: is_a_corner

          integer :: procedure_type

          procedure_type = bc_section(1)

          is_a_corner = (procedure_type.eq.NE_corner_type).or.
     $                  (procedure_type.eq.NW_corner_type).or.
     $                  (procedure_type.eq.SE_corner_type).or.
     $                  (procedure_type.eq.SW_corner_type)

        end function is_a_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the index identifying the y-position of the boundary
        !> layer
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section_sorted
        !> boundary layer represented as [procedure_type,i_min,j_min,extent]
        !--------------------------------------------------------------
        function get_j_stage(bc_section_sorted)

          implicit none

          integer, dimension(4), intent(in) :: bc_section_sorted
          integer                           :: get_j_stage
          
          get_j_stage = bc_section_sorted(3)

        end function get_j_stage


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 'nb_ele_temp' attribute
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @return nb_ele_temp
        !> 'nb_ele_temp' attribute
        !--------------------------------------------------------------
        function get_nb_ele_temp(this) result(nb_ele_temp)

          implicit none

          class(bf_layer_bc_sections), intent(in) :: this
          integer                                 :: nb_ele_temp

          nb_ele_temp = this%nb_ele_temp

        end function get_nb_ele_temp


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 'nb_ele_final' attribute
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @return nb_ele_temp
        !> 'nb_ele_final' attribute
        !--------------------------------------------------------------
        function get_nb_ele_final(this) result(nb_ele_final)

          implicit none

          class(bf_layer_bc_sections), intent(in) :: this
          integer                                 :: nb_ele_final

          nb_ele_final = this%nb_ele_final

        end function get_nb_ele_final


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the boundary layers saved in the bc_sections
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 'bc_sections_temp' attribute
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !>@return bc_sections_temp
        !> 'bc_sections_temp' attribute
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 'bc_sections_buffer' attribute
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !>@return bc_sections_buffer
        !> 'bc_sections_buffer' attribute
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 'bc_sections_final' attribute
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !  
        !> @return bc_sections_final
        !> 'bc_sections_final' attribute
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the boundary layer information in a human
        !> readable format
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param k
        !> index identifying the boundary layer printed
        !
        !> @param bc_section represented as
        !> [procedure_type,edge_min,edge_max,coord,match_nb]
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine which grid-points are computed when
        !> the anti-corner bc_section is overlap by corner
        !> bc_section
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param overlap_type
        !> index identifying how the anti-corner bc_section
        !> is overlaped by a neighboring corner bc_section
        !
        !> @param compute_point1
        !> logical indicating whether the SW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point2
        !> logical indicating whether the SE grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point3
        !> logical indicating whether the NW grid-point of
        !> the bc_section should be computed
        !
        !> @param compute_point4
        !> logical indicating whether the NE grid-point of
        !> the bc_section should be computed
        !--------------------------------------------------------------
        subroutine determine_edge_points_computed(
     $     overlap_type,
     $     compute_point1,
     $     compute_point2,
     $     compute_point3,
     $     compute_point4)

          implicit none

          integer, intent(in)  :: overlap_type
          logical, intent(out) :: compute_point1
          logical, intent(out) :: compute_point2
          logical, intent(out) :: compute_point3
          logical, intent(out) :: compute_point4


          compute_point1 = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.SW_overlap))

          compute_point2 = .not.(
     $         (overlap_type.eq.S_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.SE_overlap))

          compute_point3 = .not.(
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.W_overlap).or.
     $         (overlap_type.eq.NW_overlap))

          compute_point4 = .not.(
     $         (overlap_type.eq.N_overlap).or.
     $         (overlap_type.eq.E_overlap).or.
     $         (overlap_type.eq.NE_overlap))

        end subroutine determine_edge_points_computed

      end module bf_layer_bc_sections_class
