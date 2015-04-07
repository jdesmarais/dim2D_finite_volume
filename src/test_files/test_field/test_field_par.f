      !  _____________                    _______________
      ! |             |                  |       |       |
      ! |             |                  |       |       |
      ! |             |  compared with   |_______|_______|
      ! |             |                  |       |       |
      ! |             |                  |       |       |
      ! |_____________|                  |_______|_______|
      !
      !      field                           field_par
      !------------------------------------------------------------
      program test_field_par

        use bc_operators_class, only :
     $       bc_operators

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_real_validated,
     $       is_real_matrix3D_validated

        use field_class, only :
     $       field

        use field_par_class, only :
     $       field_par

        use mpi

        use mpi_process_class, only :
     $       mpi_process

        use parameters_constant, only :
     $       N,S,E,W,
     $       sincos,
     $       reflection_xy_choice

        use parameters_input, only :
     $       nx,ny,ne,
     $       x_min,
     $       y_min,
     $       x_max,
     $       y_max,
     $       ic_choice,
     $       bc_choice,
     $       npx,npy,
     $       ntx,nty

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class ,only :
     $       td_operators

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt

        implicit none

        real(rkind), parameter :: t     =  0.0d0
        real(rkind), parameter :: dt    =  0.05d0
        

        logical, parameter :: generate_small_domain = .false.
        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        type(mpi_process) :: mpi_process_used
        integer :: ierror
        integer :: rank


        if(generate_small_domain) then

           call generate_small_domain_results()

        else
           
           call check_inputs()
           
           detailled = .true.
           test_validated = .true.


           test_loc = test_ini(detailled)
           test_validated = test_validated.and.test_loc
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_ini: '',L1)', test_loc
              print '()'
           end if


           test_loc = test_compute_integration_step(detailled)
           test_validated = test_validated.and.test_loc
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_compute_integration_step: '',L1)', test_loc
              print '()'
           end if


           test_loc = test_integrate(detailled)
           test_validated = test_validated.and.test_loc
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_integrate: '',L1)', test_loc
              print '()'
           end if


           call test_write_data()
         
           
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_validated: '',L1)', test_validated
              print '()'
           end if


           call mpi_process_used%finalize_mpi()


        end if


        contains


        subroutine test_write_data()

          implicit none

          type(field_par) :: field_used

          call field_used%ini()
          call field_used%write_data()

        end subroutine test_write_data


        function test_integrate(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(field_par)                   :: field_used
          real(rkind), dimension(100,110,3) :: nodesInt

          logical :: test_loc
          integer :: ios

          integer :: ierror
          integer :: i
          integer :: rank
          logical, dimension(4) :: test_loc_gathered


          test_validated = .true.


          !output
          !------------------------------------------------------------
          call field_used%ini()
          call field_used%integrate(dt)


          !validation
          !------------------------------------------------------------

          !validation: read the results from the small domain
          open(unit=2,
     $         file='field_nodesInt.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodesInt
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          select case(field_used%usr_rank)
      
             !SW tile
             case(0)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodesInt(1:52,1:57,:),
     $               detailled)

             !NW tile
             case(1)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodesInt(1:52,54:110,:),
     $               detailled)
                
             !SE tile
             case(2)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodesInt(49:100,1:57,:),
     $               detailled)
                
             !NE tile
             case(3)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodesInt(49:100,54:110,:),
     $               detailled)

          end select

          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes('',I2,'') failed'')', field_used%usr_rank
          end if


          ! the master processor gathers the test results
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,4
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             test_validated = test_validated.and.test_loc

          end if

        end function test_integrate


        function test_compute_integration_step(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(field_par)                   :: field_used
          real(rkind), dimension(100,110,3) :: nodes1

          real(rkind), dimension(nx,ny,ne)  :: interior_nodes_tmp
          real(rkind), dimension(nx,ny,ne)  :: interior_timedev

          logical :: test_loc
          integer :: ios

          integer :: ierror
          integer :: i
          integer :: rank
          logical, dimension(4) :: test_loc_gathered


          test_validated = .true.


          !output
          !------------------------------------------------------------
          call field_used%ini()

          interior_timedev = field_used%compute_time_dev()

          call field_used%compute_integration_step(
     $         dt,
     $         interior_nodes_tmp,
     $         interior_timedev,
     $         compute_1st_step)


          !validation
          !------------------------------------------------------------
          !validation: read the results from the small domain
          open(unit=2,
     $         file='field_nodes1st.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes1
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          select case(field_used%usr_rank)
      
             !SW tile
             case(0)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(1:50,1:55,:),
     $               nodes1(1:50,1:55,:),
     $               .false.)

             !NW tile
             case(1)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(1:50,3:57,:),
     $               nodes1(1:50,56:110,:),
     $               .false.)
                
             !SE tile
             case(2)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(3:52,1:55,:),
     $               nodes1(51:100,1:55,:),
     $               .false.)
                
             !NE tile
             case(3)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(3:52,3:57,:),
     $               nodes1(51:100,56:110,:),
     $               .false.)

          end select

          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes('',I2,'') failed'')', field_used%usr_rank
          end if


          ! the master processor gathers the test results
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,4
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             test_validated = test_validated.and.test_loc

          end if

        end function test_compute_integration_step


        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)                 :: mpi_process_used
          type(field_par)                   :: field_used
          real(rkind), dimension(100,110,3) :: nodes0

          logical :: test_loc
          integer :: ios

          integer :: ierror
          integer :: rank
          integer :: i

          logical, dimension(4) :: test_loc_gathered


          test_validated = .true.


          call mpi_process_used%ini_mpi()


          !output
          !------------------------------------------------------------
          call field_used%ini()


          !validation
          !------------------------------------------------------------

          !validation: read the results from the small domain
          open(unit=2,
     $         file='field_nodes0.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes0
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          select case(field_used%usr_rank)
      
             !SW tile
             case(0)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(1:52,1:57,:),
     $               detailled)

             !NW tile
             case(1)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(1:52,54:110,:),
     $               detailled)
                
             !SE tile
             case(2)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(49:100,1:57,:),
     $               detailled)
                
             !NE tile
             case(3)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(49:100,54:110,:),
     $               detailled)

          end select

          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes('',I2,'') failed'')', field_used%usr_rank
          end if


          ! the master processor gathers the test results
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,4
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             test_validated = test_validated.and.test_loc
          end if
          

        end function test_ini


        subroutine generate_small_domain_results()

          implicit none

          type(field) :: field_used

          real(rkind), dimension(nx,ny,ne) :: time_dev
          real(rkind), dimension(nx,ny,ne) :: interior_nodes_tmp

          integer :: ios


          call check_inputs_small_domain()


          !> initialization
          call field_used%ini()


          !> write the nodes0 on an output file
          open(unit=2,
     $         file='field_nodes0.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes written in field_nodes0.out'')'


          !> compute the time derivatives
          time_dev = field_used%compute_time_dev()

          !> write the time derivatives on an output file
          open(unit=2,
     $         file='field_timedev.out',
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) time_dev
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''time derivatives written in field_timedev.out'')'


          !> compute the time integration step
          call field_used%compute_integration_step(
     $         dt,
     $         interior_nodes_tmp,
     $         time_dev,
     $         compute_1st_step)

          !> write the nodes1st on an output file
          open(unit=2,
     $         file='field_nodes1st.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes written in field_nodes1st.out'')'


          !> re-initialization
          call field_used%ini()

          !> complete integration cycle
          call field_used%integrate(dt)

          !> write the nodes1st on an output file
          open(unit=2,
     $         file='field_nodesInt.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes after integration written in field_nodesInt.out'')'

        end subroutine generate_small_domain_results


        subroutine check_inputs_small_domain()

          if(.not.(
     $       (npx.eq.1).and.
     $       (npy.eq.1).and.
     $       (ntx.eq.100).and.
     $       (nty.eq.110).and.
     $       (ne.eq.3).and.
     $       is_real_validated(x_min,-10.0d0,.true.).and.
     $       is_real_validated(x_max,-0.50d0,.true.).and.
     $       is_real_validated(y_min,-10.0d0,.true.).and.
     $       is_real_validated(y_max, 11.0d0,.true.).and.
     $       (ic_choice.eq.sincos).and.
     $       (bc_choice.eq.reflection_xy_choice))) then

             print '(''the test requires:'')'
             print '(''   - npx=1'')'
             print '(''   - npy=1'')'
             print '(''   - nx=100'')'
             print '(''   - ny=110'')'
             print '(''   - ne=3'')'
             print '(''   - x_min=-10'')'
             print '(''   - x_max=-0.5'')'
             print '(''   - y_min=-10'')'
             print '(''   - y_max= 11'')'
             print '(''   - wave2d model'')'
             print '(''   - sincos i.c.'')'
             print '(''   - reflection i.c.'')'
             stop ''

          end if

        end subroutine check_inputs_small_domain


        subroutine check_inputs()

          if(.not.(
     $       (npx.eq.2).and.
     $       (npy.eq.2).and.
     $       (ntx.eq.104).and.
     $       (nty.eq.114).and.
     $       (ne.eq.3).and.
     $       is_real_validated(x_min,-10.0d0,.true.).and.
     $       is_real_validated(x_max,-0.50d0,.true.).and.
     $       is_real_validated(y_min,-10.0d0,.true.).and.
     $       is_real_validated(y_max, 11.0d0,.true.).and.
     $       (ic_choice.eq.sincos).and.
     $       (bc_choice.eq.reflection_xy_choice))) then

             print '(''the test requires:'')'
             print '(''   - npx=2'')'
             print '(''   - npy=2'')'
             print '(''   - ntx=104'')'
             print '(''   - nty=114'')'
             print '(''   - ne=3'')'
             print '(''   - x_min=-10.0'')'
             print '(''   - x_max=-0.50'')'
             print '(''   - y_min=-10.0'')'
             print '(''   - y_max= 11.0'')'
             print '(''   - wave2d model'')'
             print '(''   - sincos i.c.'')'
             print '(''   - reflection i.c.'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_field_par
