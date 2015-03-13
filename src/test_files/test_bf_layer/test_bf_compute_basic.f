      !test the function td_operators%compute_time_dev_nopt()
      program test_bf_compute_basic

        use bf_compute_basic_class, only :
     $       bf_compute_basic

        use check_data_module, only :
     $       is_real_vector_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_S,align_E,
     $       interior_pt

        use parameters_input, only :
     $       nx,ny,ne, bc_size,
     $       x_min, x_max,
     $       y_min, y_max

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_does_previous_timestep_exist(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_does_previous_timestep_exist: '',L1)', test_loc
        print '()'

        
        test_loc = test_allocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_allocate_tables: '',L1)', test_loc
        print '()'


        test_loc = test_deallocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_deallocate_tables: '',L1)', test_loc
        print '()'

        contains


        function test_does_previous_timestep_exist(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute_basic) :: bf_compute_used          

          test_validated = .not.(bf_compute_used%does_previous_timestep_exist())

          if(detailled.and.(.not.test_validated)) then
             print '(''test failed: '',L1,'' -> '',L1)', .true., .false.
          end if

        end function test_does_previous_timestep_exist


        function test_allocate_tables(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute_basic)                :: bf_compute_used
          integer    , parameter                :: size_x=6
          integer    , parameter                :: size_y=10
          integer    , dimension(2,2)           :: alignment
          real(rkind)                           :: dx
          real(rkind)                           :: dy
          real(rkind), dimension(size_x)        :: x_map
          real(rkind), dimension(size_y)        :: y_map
          integer    , dimension(size_x,size_y) :: grdpts_id

          integer :: i,j

          test_validated = .true.

          !input
          alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+size_x-2*bc_size,
     $         align_S+1+size_y-2*bc_size/),
     $         (/2,2/))

          dx = (x_max-x_min)/size_x
          dy = (y_max-y_min)/size_y

          x_map = (/ (x_min+(i-1)*dx,i=1,size_x) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,size_y) /)
          grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,size_x),j=1,size_y) /),
     $         (/size_x,size_y/))

          !output
          call bf_compute_used%allocate_tables(
     $         size_x,size_y,
     $         alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          !validation
          !alignment
          test_loc = is_int_matrix_validated(
     $         bf_compute_used%alignment_tmp,alignment,detailled)
          test_validated = test_validated.and.test_loc

          !x_map
          test_loc = is_real_vector_validated(
     $         bf_compute_used%x_map_tmp,x_map,detailled)
          test_validated = test_validated.and.test_loc

          !y_map
          test_loc = is_real_vector_validated(
     $         bf_compute_used%y_map_tmp,y_map,detailled)
          test_validated = test_validated.and.test_loc

          !grdpts_id
          test_loc = is_int_matrix_validated(
     $         bf_compute_used%grdpts_id_tmp,grdpts_id,detailled)
          test_validated = test_validated.and.test_loc

          !nodes
          test_loc = size_x.eq.size(bf_compute_used%nodes_tmp,1)
          test_validated = test_validated.and.test_loc

          test_loc = size_y.eq.size(bf_compute_used%nodes_tmp,2)
          test_validated = test_validated.and.test_loc

          test_loc = ne.eq.size(bf_compute_used%nodes_tmp,3)
          test_validated = test_validated.and.test_loc
          
          !timedev
          test_loc = size_x.eq.size(bf_compute_used%time_dev,1)
          test_validated = test_validated.and.test_loc

          test_loc = size_y.eq.size(bf_compute_used%time_dev,2)
          test_validated = test_validated.and.test_loc

          test_loc = ne.eq.size(bf_compute_used%time_dev,3)
          test_validated = test_validated.and.test_loc

        end function test_allocate_tables

      
        function test_deallocate_tables(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute_basic)                :: bf_compute_used
          integer    , parameter                :: size_x=6
          integer    , parameter                :: size_y=10
          integer    , dimension(2,2)           :: alignment
          real(rkind), dimension(size_x)        :: x_map
          real(rkind), dimension(size_y)        :: y_map
          integer    , dimension(size_x,size_y) :: grdpts_id

          logical :: test_loc

          test_validated = .true.

          call bf_compute_used%deallocate_tables()

          test_loc = .not.allocated(bf_compute_used%alignment_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%grdpts_id_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%x_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%y_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%nodes_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%time_dev)
          test_validated = test_validated.and.test_loc
          
          if(detailled.and.(.not.test_validated)) then
             print '(''test_deallocate_tables(1) failed'')'
          end if

          call bf_compute_used%allocate_tables(
     $         size_x,size_y,
     $         alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          call bf_compute_used%deallocate_tables()

          test_loc = .not.allocated(bf_compute_used%alignment_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%grdpts_id_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%x_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%y_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%nodes_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%time_dev)
          test_validated = test_validated.and.test_loc

        end function test_deallocate_tables


        subroutine check_inputs()

          implicit none

          logical :: test_parameters

          test_parameters = test_parameters.and.(nx.eq.6)
          test_parameters = test_parameters.and.(ny.eq.6)
          
        end subroutine check_inputs

      end program test_bf_compute_basic
