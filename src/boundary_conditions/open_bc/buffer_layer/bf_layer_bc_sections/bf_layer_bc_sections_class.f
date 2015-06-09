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

        use bf_layer_bc_procedure_module, only : 
     $       get_bc_interior_pt_procedure

        use bf_layer_bc_sections_overlap_module, only :
     $       overlap_square_bc_sections,
     $       overlap_bc_section_by_integration_borders

        use bf_layer_errors_module, only :
     $       error_bc_section_type,
     $       error_overlap_index,
     $       error_overlap_incompatible

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table,
     $       get_indices_to_extract_interior_data,
     $       get_indices_to_extract_bf_layer_data,
     $       get_grdpts_id_from_interior

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       bc_pt,
     $       BF_SUCCESS,
     $       no_overlap
        
        use parameters_constant, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind

        implicit none

        private
        public :: bf_layer_bc_sections


        integer, parameter :: max_bc_sections_temp = 6
        integer, parameter :: reallocation_ele_nb  = 5
        


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
        !> @param check_square_overlap
        !> analyze the ordered bc_sections and mark the overlapping
        !> gridpoints b/w corner and anti-corner boundary layers
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
          procedure, nopass :: create_tmp_grdpts_id_for_analyse
          procedure,   pass :: analyse_grdpt

          procedure, nopass :: bubble_sort
          procedure,   pass :: sort_bc_sections
          procedure, nopass :: check_square_overlap
          procedure, nopass :: check_integration_overlap
          procedure, nopass :: check_overlaps
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
          integer :: l_buffer

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
          do l=k+1, this%nb_ele_temp

             l_buffer = l-size(this%bc_sections_temp,2)

             this%bc_sections_buffer(:,l_buffer-1) =
     $            this%bc_sections_buffer(:,l_buffer)
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
        function get_bc_section(
     $     i,j,grdpts_id,
     $     tmp_grdpts_id, use_tmp_grdpts_id,
     $     ierror)
     $     result(bc_section)

          implicit none

          integer(ikind)                , intent(in)  :: i
          integer(ikind)                , intent(in)  :: j
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer       , dimension(:,:), intent(in)  :: tmp_grdpts_id
          logical                       , intent(in)  :: use_tmp_grdpts_id
          logical                       , intent(out) :: ierror
          integer       , dimension(5)                :: bc_section

          integer :: procedure_type
          integer :: i_proc
          integer :: j_proc

          ierror = BF_SUCCESS                    


          ! analyse the procedure from the grdpts_id
          if(use_tmp_grdpts_id) then

             call get_bc_interior_pt_procedure(
     $            2,2,
     $            tmp_grdpts_id,
     $            procedure_type,
     $            i_proc,
     $            j_proc,
     $            ierror)

             ! the i_proc and j_proc extracted from
             ! tmp_grdpts_id are expressed in a reference
             ! frame where (2,2) <-> (i,j)
             ! we need to turn (i_proc,j_proc) bask to the
             ! general frame
             ! i_gen = (i_loc-2) + i
             ! j_gen = (j_loc-2) + j
             !------------------------------------------------------------
             i_proc = i_proc-2+i
             j_proc = j_proc-2+j

          else

             call get_bc_interior_pt_procedure(
     $            i,j,
     $            grdpts_id,
     $            procedure_type,
     $            i_proc,
     $            j_proc,
     $            ierror)

          end if


          ! determine the bc_section from the procedure
          if(ierror.eqv.BF_SUCCESS) then
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

                case(SE_corner_type, SW_corner_type,
     $               NE_corner_type, NW_corner_type,
     $               SE_edge_type  , SW_edge_type,
     $               NE_edge_type  , NW_edge_type)
                  bc_section(2)=i_proc
                  bc_section(3)=j_proc

                case default
                   print '(''bf_layer_bc_sections'')'
                   print '(''get_bc_section'')'
                   print '(''procedure type: '',I3)', procedure_type
                   print '(''procedure not recognized'')'
                   print '(''****************************************'')'
                   print '()'
                   ierror = .not.BF_SUCCESS
             
              end select

           else
              print '(''bf_layer_bc_sections'')'
              print '(''get_bc_section'')'
              print '(''get_bc_interior_pt_procedure failed'')'
              print '(''****************************************'')'
              print '()'
              ierror = .not.BF_SUCCESS
           end if

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
     $     i,j,grdpts_id,
     $     tmp_grdpts_id,use_tmp_grdpts_id,
     $     bc_section,remove_ele)
     $     result(compatible)

          implicit none

          integer                , intent(in)    :: i
          integer                , intent(in)    :: j
          integer, dimension(:,:), intent(in)    :: grdpts_id
          integer, dimension(:,:), intent(in)    :: tmp_grdpts_id
          logical                , intent(in)    :: use_tmp_grdpts_id
          integer, dimension(5)  , intent(inout) :: bc_section
          logical                , intent(out)   :: remove_ele
          logical                                :: compatible

          remove_ele = .false.
          compatible = .false.


          !type of procedure
          select case(bc_section(1))

            !bc_section(2):i_min
            !bc_section(3):i_max
            !bc_section(4):j
            case(N_edge_type)
               compatible =
     $              ((j-bc_section(4)).eq.0).and.
     $              ((i-bc_section(3)).eq.1)

               if(compatible) then
                  if(use_tmp_grdpts_id) then
                     compatible = 
     $                    (tmp_grdpts_id(2,2+1).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2+1,2+1).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2+1,2).eq.bc_interior_pt)
                  else
                     compatible = 
     $                    (grdpts_id(i,j+1).eq.bc_pt).and.
     $                    (grdpts_id(i+1,j+1).eq.bc_pt).and.
     $                    (grdpts_id(i+1,j).eq.bc_interior_pt)
                  end if
                     
                  if(compatible) then
                     bc_section(3)=i
                  end if

               end if

               remove_ele = j.gt.bc_section(4)
               

            !bc_section(2):i_min
            !bc_section(3):i_max
            !bc_section(4):j
            case(S_edge_type)
               compatible = 
     $              ((j-bc_section(4)).eq.0).and.
     $              ((i-bc_section(3)).eq.1)

               if(compatible) then
                  if(use_tmp_grdpts_id) then
                     compatible = 
     $                    (tmp_grdpts_id(2,2-1).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2+1,2-1).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2+1,2).eq.bc_interior_pt)
                  else
                     compatible = 
     $                    (grdpts_id(i,j-1).eq.bc_pt).and.
     $                    (grdpts_id(i+1,j-1).eq.bc_pt).and.
     $                    (grdpts_id(i+1,j).eq.bc_interior_pt)
                  end if

                  if(compatible) then
                     bc_section(3)=i
                  end if
               end if

               remove_ele = j.gt.bc_section(4)

            !bc_section(2): j_min
            !bc_section(3): j_max
            !bc_section(4): i
            case(E_edge_type)
               compatible =
     $              ((i-bc_section(4)).eq.0).and.
     $              ((j-bc_section(3)).eq.1)

               if(compatible) then
                  if(use_tmp_grdpts_id) then
                     compatible =
     $                    (tmp_grdpts_id(2+1,2).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2,2+1).eq.bc_interior_pt).and.
     $                    (tmp_grdpts_id(2+1,2+1).eq.bc_pt)
                  else
                     compatible = 
     $                    (grdpts_id(i+1,j).eq.bc_pt).and.
     $                    (grdpts_id(i,j+1).eq.bc_interior_pt).and.
     $                    (grdpts_id(i+1,j+1).eq.bc_pt)
                  end if
                     
                  if(compatible) then
                     bc_section(3)=j
                  end if
               end if

               remove_ele=j.gt.(bc_section(3)+1)


            !bc_section(2): j_min
            !bc_section(3): j_max
            !bc_section(4): i
            case(W_edge_type)
               compatible = 
     $              ((i-bc_section(4)).eq.0).and.
     $              ((j-bc_section(3)).eq.1)

               if(compatible) then
                  if(use_tmp_grdpts_id) then
                     compatible = 
     $                    (tmp_grdpts_id(2-1,2).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2-1,2+1).eq.bc_pt).and.
     $                    (tmp_grdpts_id(2,2+1).eq.bc_interior_pt)
                  else
                     compatible = 
     $                    (grdpts_id(i-1,j).eq.bc_pt).and.
     $                    (grdpts_id(i-1,j+1).eq.bc_pt).and.
     $                    (grdpts_id(i,j+1).eq.bc_interior_pt)
                  end if

                  if(compatible) then
                     bc_section(3)=j
                  end if
               end if

               remove_ele=j.gt.(bc_section(3)+1)


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NE_edge_type)
               if(use_tmp_grdpts_id) then
                  compatible = 
     $                 ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).and.(tmp_grdpts_id(2+1,2).eq.bc_interior_pt)).or.
     $                 (((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).and.(tmp_grdpts_id(2,2+1).eq.bc_interior_pt))
               else
                  compatible = 
     $                 ((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).and.(grdpts_id(i+1,j).eq.bc_interior_pt)).or.
     $                 (((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i,j+1).eq.bc_interior_pt))
               end if
                  
               remove_ele = j.gt.(bc_section(3)+1)


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SE_edge_type)
               if(use_tmp_grdpts_id) then
                  compatible = 
     $                 (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(tmp_grdpts_id(2,2-1).eq.bc_interior_pt)).or.
     $                 ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(tmp_grdpts_id(2+1,2).eq.bc_interior_pt))
               else
                  compatible = 
     $                 (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(grdpts_id(i,j-1).eq.bc_interior_pt)).or.
     $                 ((i.eq.bc_section(2)).and.(j.eq.(bc_section(3)+1))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i+1,j).eq.bc_interior_pt))
               end if

               remove_ele = j.gt.(bc_section(3)+1)


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(SW_edge_type)
               if(use_tmp_grdpts_id) then
                  compatible = 
     $                 (((i.eq.bc_section(2)+1).and.(j.eq.(bc_section(3)))).and.(tmp_grdpts_id(2,2-1).eq.bc_interior_pt)).or.
     $                 (((i.eq.(bc_section(2))).and.(j.eq.(bc_section(3)+1))).and.(tmp_grdpts_id(2-1,2).eq.bc_interior_pt)).or.
     $                 ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3)+1))
               else
                  compatible = 
     $                 (((i.eq.bc_section(2)+1).and.(j.eq.(bc_section(3)))).and.(grdpts_id(i,j-1).eq.bc_interior_pt)).or.
     $                 (((i.eq.(bc_section(2))).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i-1,j).eq.bc_interior_pt)).or.
     $                 ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3)+1))
               end if

               remove_ele = j.gt.(bc_section(3)+1)


            !bc_section(2): i_min
            !bc_section(3): j_min
            !bc_section(5): match_nb
            case(NW_edge_type)
               if(use_tmp_grdpts_id) then
                  compatible = 
     $                 (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(tmp_grdpts_id(2-1,2).eq.bc_interior_pt)).or.
     $                 ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(tmp_grdpts_id(2,2+1).eq.bc_interior_pt))
               else
                  compatible = 
     $                 (((i.eq.bc_section(2)).and.(j.eq.bc_section(3))).and.(grdpts_id(i-1,j).eq.bc_interior_pt)).or.
     $                 ((i.eq.(bc_section(2)+1)).and.(j.eq.bc_section(3))).or.
     $                 (((i.eq.(bc_section(2)+1)).and.(j.eq.(bc_section(3)+1))).and.(grdpts_id(i,j+1).eq.bc_interior_pt))
               end if

               remove_ele = j.gt.(bc_section(3)+1)

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
        subroutine analyse_grdpt(
     $     this,
     $     bf_alignment,
     $     i,j,grdpts_id,
     $     ierror)

          implicit none

          class(bf_layer_bc_sections)   , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: bf_alignment
          integer(ikind)                , intent(in)    :: i
          integer(ikind)                , intent(in)    :: j
          integer       , dimension(:,:), intent(in)    :: grdpts_id
          logical                       , intent(out)   :: ierror

          integer               :: k,k_buffer
          integer, dimension(5) :: new_bc_section
          logical               :: compatible
          logical               :: remove_ele

          integer, dimension(3,3) :: tmp_grdpts_id
          logical                 :: use_tmp_grdpts_id


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
          

          ! determine whether the grid-point analyzed is at the edge
          ! of the grdpts_id array. In this case, for the analysis,
          ! we need to create a temporary array containing the grdpts_id
          ! around the grid-point asked
          use_tmp_grdpts_id = .false.

          if((i.eq.1).or.(i.eq.size(grdpts_id,1)).or.(j.eq.1).or.(j.eq.size(grdpts_id,2))) then

             tmp_grdpts_id = create_tmp_grdpts_id_for_analyse(
     $            bf_alignment,
     $            [i,j],grdpts_id)

             use_tmp_grdpts_id = .true.

          end if


          !if there are already boundary layers, check whether the
          !grid point analysed is compatible with one of them
          if(this%nb_ele_temp.gt.0) then
             
             compatible = .false.
             
             ! is the grid point compatible with the boundary layers
             ! saved in this%bc_sections_temp ?
             k=1
             do while (k.le.min(this%nb_ele_temp,size(this%bc_sections_temp,2)))
                
                compatible = this%analyse_grdpt_with_bc_section(
     $               i,j,grdpts_id,
     $               tmp_grdpts_id, use_tmp_grdpts_id,
     $               this%bc_sections_temp(:,k),
     $               remove_ele)
             
                if(remove_ele) then
                   call this%add_to_final_bc_sections(this%bc_sections_temp(:,k))
                   call this%remove_from_bc_sections_temp(k)
                   k = k-1
                end if
             
                if(compatible) then
                   ierror = BF_SUCCESS
                   exit
                end if
                
                k=k+1
             
             end do
             
             if(.not.compatible) then
             
                ! is the grid point compatible with the boundary layers
                ! saved in this%bc_sections_buffer ?
                k=size(this%bc_sections_temp,2)+1
                do while(k.le.this%nb_ele_temp)
                      
                   k_buffer = k-size(this%bc_sections_temp,2)
             
                   compatible = this%analyse_grdpt_with_bc_section(
     $                  i,j,grdpts_id,
     $                  tmp_grdpts_id, use_tmp_grdpts_id,
     $                  this%bc_sections_buffer(:,k_buffer),
     $                  remove_ele)
             
                   if(remove_ele) then
                      call this%add_to_final_bc_sections(this%bc_sections_buffer(:,k_buffer))
                      call this%remove_from_bc_sections_buffer(k)
                      k = k-1
                   end if
             
                   if(compatible) then
                      ierror = BF_SUCCESS
                      exit
                   end if
                   
                   k=k+1
                      
                end do
             
             
                ! if no boundary layer matches the current grid point
                ! the grid point is used as the starting point of a new
                ! boundary layer
                if(.not.compatible) then
             
                   new_bc_section = this%get_bc_section(
     $                  i,j,grdpts_id,
     $                  tmp_grdpts_id, use_tmp_grdpts_id,
     $                  ierror)

                   call this%add_to_bc_sections(new_bc_section)
             
                end if
             
             end if
             
          !if there is no existing boundary layers, the grid point
          !is analysed to know what type of starting point for a
          !boundary layer it is and then added to the boundary
          !layers of the object
          else
             
             new_bc_section = this%get_bc_section(
     $            i,j,grdpts_id,
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            ierror)

             call this%add_to_bc_sections(new_bc_section)
             
          end if

        end subroutine analyse_grdpt


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
             allocate(sorted_bc_sections(5,nb_ele_tot))
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

          integer, dimension(5), intent(in) :: p
          integer, dimension(5), intent(in) :: q
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

          integer, dimension(5), intent(inout) :: p
          integer, dimension(5), intent(inout) :: q
          integer, dimension(5)                :: temp

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
          integer, dimension(5)             :: sorted_ele

          select case(bc_section(1))

            case(N_edge_type)
               sorted_ele = [bc_section(1),bc_section(2),bc_section(4),bc_section(3),no_overlap]

            case(S_edge_type)
               sorted_ele = [bc_section(1),bc_section(2),bc_section(4)-1,bc_section(3),no_overlap]

            case(E_edge_type)
               sorted_ele = [bc_section(1),bc_section(4),bc_section(2),bc_section(3),no_overlap]

            case(W_edge_type)
               sorted_ele = [bc_section(1),bc_section(4)-1,bc_section(2),bc_section(3),no_overlap]

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
     $           no_overlap,
     $           no_overlap]

            case default
               print '(''bf_layer_bc_sections_class'')'
               print '(''get_sorted_ele'')'
               print '(''case not recognized: '',I10)', bc_section(1)
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
        !
        !> @param x_borders
        !> x-borders for the time integration of the gridpoints
        !
        !> @param y_borders
        !> y-borders for the time integration of the gridpoints
        !
        !> @param N_bc_sections
        !> x-borders of the North gridpoints [size(2)-bc_size+1,size(2)]
        !> computed with the time integration scheme
        !
        !> @param S_bc_sections
        !> x-borders of the South gridpoints [1,bc_size] computed with
        !> the time integration scheme
        !--------------------------------------------------------------
        subroutine check_overlaps(
     $     bc_sections_sorted,
     $     x_borders,
     $     y_borders)

          implicit none

          integer, dimension(:,:), intent(inout) :: bc_sections_sorted
          integer, dimension(2)  , intent(in)    :: x_borders
          integer, dimension(2)  , intent(in)    :: y_borders

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

             !=============================================
             ! checking the overlap between the corner
             ! b.c. sections and the anti-corner b.c.
             ! sections
             !=============================================
             ! if the bc_section analyzed is a corner,
             ! it should be compared to the bc_sections
             ! of the previous and the next stages
             ! corresponding to the stages where an
             ! overlap is possible
             if(is_a_square(bc_sections_sorted(:,k))) then
                
                call check_square_overlap(
     $               bc_sections_sorted(:,k),
     $               bc_sections_sorted,
     $               k_prev_stage,
     $               k-1)

             end if

             !=============================================
             ! checking the overlap between the bc_section
             ! and the time integration borders
             !=============================================
             call check_integration_overlap(
     $            bc_sections_sorted(:,k),
     $            x_borders,
     $            y_borders)

          end do

        end subroutine check_overlaps


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
        subroutine check_square_overlap(
     $     square,
     $     bc_sections_sorted,
     $     k_min,
     $     k_max)

          implicit none

          integer, dimension(5)  , intent(inout) :: square
          integer, dimension(:,:), intent(inout) :: bc_sections_sorted
          integer                , intent(in)    :: k_min
          integer                , intent(in)    :: k_max

          integer :: k

          do k=k_min,k_max

             if(is_a_square(bc_sections_sorted(:,k))) then

                call overlap_square_bc_sections(square,bc_sections_sorted(:,k))

             end if

          end do

        end subroutine check_square_overlap        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> overlap the bc_section with the time integration borders
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bc_section
        !> boundary section represented as [bc_section_type,i_min,j_min,extent,overlap]
        
        !> @param x_borders
        !> [i_min,i_max] in the time integration
        !
        !> @param y_borders
        !> [j_min,j_max] in the time integration
        !--------------------------------------------------------------
        subroutine check_integration_overlap(
     $     bc_section,
     $     x_borders,
     $     y_borders)

          implicit none

          integer, dimension(5), intent(inout) :: bc_section
          integer, dimension(2), intent(in)    :: x_borders
          integer, dimension(2), intent(in)    :: y_borders

          
          integer, dimension(2,2) :: gen_borders
          integer, dimension(2,2) :: bc_section_borders

          gen_borders(1,1) = x_borders(1)
          gen_borders(1,2) = x_borders(2)
          gen_borders(2,1) = y_borders(1)
          gen_borders(2,2) = y_borders(2)


          select case(bc_section(1))
            case(N_edge_type,S_edge_type)
               bc_section_borders(1,1) = bc_section(2)
               bc_section_borders(2,1) = bc_section(3)
               bc_section_borders(1,2) = bc_section(4)
               bc_section_borders(2,2) = bc_section(3)+1

               call overlap_bc_section_by_integration_borders(
     $              bc_section(1),
     $              bc_section_borders,
     $              bc_section(5),
     $              gen_borders)

               bc_section(2) = bc_section_borders(1,1)
               bc_section(3) = bc_section_borders(2,1)
               bc_section(4) = bc_section_borders(1,2)

            case(E_edge_type,W_edge_type)
               bc_section_borders(1,1) = bc_section(2)
               bc_section_borders(2,1) = bc_section(3)
               bc_section_borders(1,2) = bc_section(2)+1
               bc_section_borders(2,2) = bc_section(4)

               call overlap_bc_section_by_integration_borders(
     $              bc_section(1),
     $              bc_section_borders,
     $              bc_section(5),
     $              gen_borders)

               bc_section(2) = bc_section_borders(1,1)
               bc_section(3) = bc_section_borders(2,1)
               bc_section(4) = bc_section_borders(2,2)

            case(   NE_corner_type, NW_corner_type,
     $              SE_corner_type, SW_corner_type,
     $              NE_edge_type,   NW_edge_type,
     $              SE_edge_type,   SW_edge_type)

               bc_section_borders(1,1) = bc_section(2)
               bc_section_borders(2,1) = bc_section(3)
               bc_section_borders(1,2) = bc_section(2)+1
               bc_section_borders(2,2) = bc_section(3)+1

               call overlap_bc_section_by_integration_borders(
     $              bc_section(1),
     $              bc_section_borders,
     $              bc_section(5),
     $              gen_borders)

            case default
               call error_bc_section_type(
     $              'bf_layer_bc_sections_class',
     $              'check_overlap_integration_borders',
     $              bc_section(1))

          end select

        end subroutine check_integration_overlap


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
        function is_a_square(bc_section)

          implicit none

          integer, dimension(5), intent(in) :: bc_section
          logical                           :: is_a_square

          integer :: procedure_type

          procedure_type = bc_section(1)

          is_a_square = (procedure_type.eq.NE_corner_type).or.
     $                  (procedure_type.eq.NW_corner_type).or.
     $                  (procedure_type.eq.SE_corner_type).or.
     $                  (procedure_type.eq.SW_corner_type).or.
     $                  (procedure_type.eq.NE_edge_type).or.
     $                  (procedure_type.eq.NW_edge_type).or.
     $                  (procedure_type.eq.SE_edge_type).or.
     $                  (procedure_type.eq.SW_edge_type)

        end function is_a_square        


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
        !> finalize the bc_sections analyzed by the
        !> bf_layer_bc_sections object:
        !> 1) keep only the information needed to compute the
        !>    boundary grid-points
        !> 2) sort the bc_sections in increasing j and i to
        !>    improve cache efficiency when computing the
        !>    boundary layer gridpoints
        !> 3) determine the overlap b/w the corner and anti-corner
        !>    boundary procedures to prevent interactions and symetry
        !>    violation when applying the boundary procedures
        !> 4) determine the overlap b/w the anti-corners and the
        !>    time integration borders
        !> 5) deallocate the intermediate attributes used to analyze
        !>    the boundary layers
        !
        !> @date
        !> 26_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer_bc_sections object encapsulating the
        !> localization of the boundary layers
        !
        !> @param x_borders
        !> x-borders for the time integration of the gridpoints
        !
        !> @param y_borders
        !> y-borders for the time integration of the gridpoints
        !
        !> @param N_bc_sections
        !> x-borders of the North gridpoints [size(2)-bc_size+1,size(2)]
        !> computed with the time integration scheme
        !
        !> @param S_bc_sections
        !> x-borders of the South gridpoints [1,bc_size] computed with
        !> the time integration scheme
        !
        !> @param sorted_bc_sections
        !> array with the boundary sections sorted in increasing j
        !> and increasing i
        !--------------------------------------------------------------
        subroutine finalize_bc_sections(
     $     this,
     $     x_borders,
     $     y_borders,
     $     sorted_bc_sections)

          implicit none

          class(bf_layer_bc_sections)                , intent(inout) :: this
          integer(ikind), dimension(2)               , intent(in)    :: x_borders
          integer(ikind), dimension(2)               , intent(in)    :: y_borders
          integer       , dimension(:,:), allocatable, intent(out)   :: sorted_bc_sections
          

          !> 1) keep only the information needed to compute the
          !>    boundary grid-points
          !> 2) sort the bc_sections in increasing j and i to
          !>    improve cache efficieny when computing the
          !>    boundary layer gridpoints
          call sort_bc_sections(this,sorted_bc_sections)


          !4) deallocate the intermediate attributes used to analyze
          !>  the boundary layers
          call this%deallocate_tables()


          !3) determine the overlap b/w the corner and anti-corner
          !>  boundary procedures to prevent interactions and symetry
          !>  violation when applying the boundary procedures
          call check_overlaps(sorted_bc_sections,x_borders,y_borders)

        end subroutine finalize_bc_sections


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
        !> create a temporary array of grdpts_id around the
        !> grid-point (i,j)
        !
        !> @date
        !> 01_04_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !> @param central_coords
        !> indices identifying the central grid-point along the
        !> x- and y- directions
        !
        !> @param grdpts_id
        !> confiugration of the grid-points in the buffer layer
        !--------------------------------------------------------------
        function create_tmp_grdpts_id_for_analyse(
     $     bf_alignment,
     $     central_coords,grdpts_id)
     $     result(tmp_grdpts_id)

           implicit none

           integer(ikind), dimension(2,2), intent(in) :: bf_alignment
           integer(ikind), dimension(2)  , intent(in) :: central_coords
           integer       , dimension(:,:), intent(in) :: grdpts_id
           integer       , dimension(3,3)             :: tmp_grdpts_id

           integer(ikind), dimension(2)   :: match_table
           integer(ikind), dimension(2,2) :: gen_coords
           integer(ikind)                 :: i,j
           integer(ikind)                 :: size_x,size_y
           integer(ikind)                 :: i_recv,i_send,j_recv,j_send          


           ! determination of the borders for the extraction
           !============================================================
           match_table = get_bf_layer_match_table(bf_alignment)

           gen_coords(1,1) = central_coords(1) -1 + match_table(1)
           gen_coords(2,1) = central_coords(2) -1 + match_table(2)
           gen_coords(1,2) = central_coords(1) +1 + match_table(1)
           gen_coords(2,2) = central_coords(2) +1 + match_table(2)


           ! grid-point extraction
           !============================================================
           ! from the interior
           call get_grdpts_id_from_interior(tmp_grdpts_id,gen_coords)

           ! from the buffer layer
           call get_indices_to_extract_bf_layer_data(
     $          bf_alignment,
     $          gen_coords,
     $          size_x, size_y,
     $          i_recv, j_recv,
     $          i_send, j_send)

           do j=1, size_y
              do i=1, size_x
                 
                 tmp_grdpts_id(i_recv+i-1,j_recv+j-1) = 
     $                grdpts_id(i_send+i-1,j_send+j-1)
 
              end do
           end do

        end function create_tmp_grdpts_id_for_analyse

      end module bf_layer_bc_sections_class
