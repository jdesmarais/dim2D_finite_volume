      !> @file
      !> test file for the object 'fv_operators'
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
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_td_operators

        use bc_operators_class , only :
     $       bc_operators

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       interior_pt

        use parameters_constant, only :
     $       periodic_xy_choice,
     $       bc_nodes_choice

        use parameters_input, only :
     $       nx,ny,ne,bc_size,bc_choice,
     $       x_min,y_min,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class , only :
     $       sd_operators

        use td_operators_class , only :
     $       td_operators        

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        
        call check_inputs()

        detailled = .true.
        test_validated = .true.


        test_loc = test_compute_time_dev(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_compute_time_dev: '',L1)', test_loc
        print '()'


        test_loc = test_compute_time_dev_nopt(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_compute_time_dev_nopt: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_loc

        contains


        function test_compute_time_dev(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)                      :: t
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind)                      :: dx
          real(rkind)                      :: dy
          type(sd_operators)               :: s
          type(pmodel_eq)                  :: p_model
          type(bc_operators)               :: bc_used
          type(td_operators)               :: t_operator

          real(rkind), dimension(nx,ny,ne) :: time_dev
          real(rkind), dimension(nx,ny,ne) :: time_dev_test

          integer(ikind)                   :: i,j
          logical                          :: test_loc

          test_validated = .true.

          !input
          t=0.0d0

          dx=1.0
          dy=1.0

          x_map = (/ (x_min+(i-1)*dx,i=1,nx) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,ny) /)

          nodes = reshape(
     $         (/ ((i+(j-1)*nx,i=1,nx),j=1,ny) /),
     $         (/nx,ny,ne/))

          time_dev_test = reshape(
     $         (/ ((-20.0d0,i=1,nx),j=1,ny) /),
     $         (/nx,ny,ne/))

          !output
          time_dev = t_operator%compute_time_dev(
     $         t,nodes,x_map,y_map,s,p_model,bc_used)

          !validation
          test_loc = is_real_matrix3D_validated(
     $         time_dev(bc_size+1:nx-bc_size,bc_size+1:ny-bc_size,:),
     $         time_dev_test(bc_size+1:nx-bc_size,bc_size+1:ny-bc_size,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          
        end function test_compute_time_dev


        function test_compute_time_dev_nopt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)                          :: t
          real(rkind), dimension(nx)           :: x_map
          real(rkind), dimension(ny)           :: y_map
          integer    , dimension(nx,ny)        :: grdpts_id
          real(rkind), dimension(nx,ny,ne)     :: nodes
          real(rkind), dimension(nx,ny,ne)     :: interior_nodes
          real(rkind)                          :: dx
          real(rkind)                          :: dy
          type(sd_operators)                   :: s
          type(pmodel_eq)                      :: p_model
          type(bc_operators)                   :: bc_used
          type(td_operators)                   :: t_operator
                                               
          real(rkind), dimension(nx,ny,ne)     :: time_dev
          real(rkind), dimension(nx,ny,ne)     :: time_dev_test
                                               
          integer(ikind)                       :: i,j
          logical                              :: test_loc
          integer, dimension(:,:), allocatable :: bc_sections

          test_validated = .true.

          !input
          t=0.0d0

          dx=1.0
          dy=1.0

          x_map = (/ (x_min+(i-1)*dx,i=1,nx) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,ny) /)

          grdpts_id = reshape(
     $         (/ ((interior_pt,i=1,nx),j=1,ny) /),
     $         (/nx,ny/))

          nodes = reshape(
     $         (/ ((i+(j-1)*nx,i=1,nx),j=1,ny) /),
     $         (/nx,ny,ne/))

          interior_nodes = reshape(
     $         (/ ((-99.0d0,i=1,nx),j=1,ny) /),
     $         (/nx,ny,ne/))

          time_dev_test = reshape(
     $         (/ ((-20.0d0,i=1,nx),j=1,ny) /),
     $         (/nx,ny,ne/))

          !output
          call t_operator%compute_time_dev_nopt(
     $         t,nodes,x_map,y_map,s,p_model,bc_used,time_dev,
     $         reshape((/2,2,2,2/),(/2,2/)),
     $         grdpts_id,
     $         interior_nodes,
     $         bc_sections,
     $         x_borders=[bc_size+1,nx-bc_size],
     $         y_borders=[bc_size+1,ny-bc_size])

          !validation
          test_loc = is_real_matrix3D_validated(
     $         time_dev(bc_size+1:nx-bc_size,bc_size+1:ny-bc_size,:),
     $         time_dev_test(bc_size+1:nx-bc_size,bc_size+1:ny-bc_size,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          
        end function test_compute_time_dev_nopt


        subroutine check_inputs()

          implicit none

          logical :: test_parameter

          test_parameter=.true.
          test_parameter=test_parameter.and.(nx.eq.10)
          test_parameter=test_parameter.and.(ny.eq.6)
          test_parameter=test_parameter.and.(ne.eq.1)
          test_parameter=test_parameter.and.(bc_choice.eq.periodic_xy_choice)
          test_parameter=test_parameter.and.(bc_N_type_choice.eq.bc_nodes_choice)
          test_parameter=test_parameter.and.(bc_S_type_choice.eq.bc_nodes_choice)    
          test_parameter=test_parameter.and.(bc_E_type_choice.eq.bc_nodes_choice)    
          test_parameter=test_parameter.and.(bc_W_type_choice.eq.bc_nodes_choice)    
          if(.not.test_parameter) then
             print *, 'the test requires several parameters'
             print *, 'test designed for simpletest eq'
             print *, 'nx=10'
             print *, 'ny=6'
             print *, 'ne=1'
             print *, 'bc_choice=periodic_xy_choice'
             print *, 'bc_N_type_choice=bc_nodes_choice'
             print *, 'bc_S_type_choice=bc_nodes_choice'
             print *, 'bc_E_type_choice=bc_nodes_choice'
             print *, 'bc_W_type_choice=bc_nodes_choice'
             stop ''
          end if

        end subroutine check_inputs


        function check_around_timedev(
     $     timedev,
     $     gen_coords,
     $     cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: timedev
          integer(ikind), dimension(2,2)  , intent(in) :: gen_coords
          real(rkind)                     , intent(in) :: cst
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          integer(ikind) :: i,j
          integer        :: k

          
          test_validated = .true.


          do k=1, size(timedev,3)
             do j=1, size(timedev,2)
                do i=1, size(timedev,1)

                   test_loc = is_real_validated(timedev(i,j,k),cst,.false.)

                   if(.not.test_loc) then
                      test_loc =
     $                     (i.ge.gen_coords(1,1)).and.
     $                     (i.le.gen_coords(1,2)).and.
     $                     (j.ge.gen_coords(2,1)).and.
     $                     (j.le.gen_coords(2,2))
                      
                   end if

                   if(detailled.and.(.not.test_loc)) then
                      print *, timedev(i,j,k)
                      print *, cst
                      print '(''['',3I4,'']'')', i,j,k
                   end if

                   test_validated = test_validated.and.test_loc

                end do
             end do
          end do

        end function check_around_timedev

      end program test_td_operators
