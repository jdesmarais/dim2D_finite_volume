      program test_verify_y_symmetry

       use cmd_verify_symmetry_class, only :
     $       cmd_verify_symmetry

        use netcdf

        use nf90_operators_module, only :
     $       nf90_close_file

        use nf90_operators_read_module, only :
     $       nf90_open_file_for_reading,
     $       nf90_get_varid,
     $       nf90_read_borders,
     $       nf90_get_var_model,
     $       nf90_get_var_model_nopt

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_constant, only :
     $       vector_y
        
        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none

        logical, parameter :: detailled = .true.        

        type(cmd_verify_symmetry) :: cmd_operators_used
        character(len=1024)       :: filename1
        character(len=1024)       :: filename2
        logical                   :: analyze_grdpts_id
        integer, dimension(4)     :: nb_bf_layers

        !analyse command line arguments
        call cmd_operators_used%analyze_cmd_line_arg()

        if(cmd_operators_used%is_analyze_activated()) then
           filename1 = cmd_operators_used%get_symmetry_filename1()

        else
           print '(''no file to be analyzed'')'
           stop ''

        end if

        analyze_grdpts_id = cmd_operators_used%get_analyze_bf_file()

        
        !if there are two files to analyze
        if(cmd_operators_used%get_analyze_two_files()) then

           filename2 = cmd_operators_used%get_symmetry_filename2()

           call verify_field_y_symmetry_two_files(
     $          filename1,
     $          filename2,
     $          analyze_grdpts_id,
     $          detailled)

        !if there is only one file to analyze
        else

           call verify_field_y_symmetry(
     $          filename1,
     $          analyze_grdpts_id,
     $          detailled)

        end if

        contains


        !check if two doubles are the same
        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: detailled
          logical                 :: test_validated


          if(detailled) then
             print *, var
             print *, cst
          end if

          test_validated = abs(var-cst)<1e-10

          if(.not.test_validated) then

             print *, 'var', var
             print *, 'cst', cst

          end if
          
        end function is_test_validated


        !check the y-symmetry in the field
        subroutine verify_field_y_symmetry(
     $       filename,
     $       analyze_grdpts_id,
     $       detailled)

          implicit none
          
          character(len=1024), intent(in) :: filename
          logical            , intent(in) :: analyze_grdpts_id
          logical            , intent(in) :: detailled

          type(pmodel_eq)                            :: p_model
          real(rkind)                                :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes
          integer    , dimension(:,:)  , allocatable :: grdpts_id


          !extract netcdf data
          call extract_data(
     $         filename1,
     $         analyze_grdpts_id,
     $         time,
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id)


          !check whether symmetry is respected along the y-axis
          if(analyze_grdpts_id) then
             call check_y_symmetry(
     $            p_model,
     $            nodes,
     $            detailled,
     $            grdpts_id=grdpts_id)

          else
             call check_y_symmetry(
     $            p_model,
     $            nodes,detailled)

          end if

          !deallocate tables
          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)
          if(analyze_grdpts_id) then
             deallocate(grdpts_id)
          end if

        end subroutine verify_field_y_symmetry



        !check the y-symmetry in the field
        subroutine verify_field_y_symmetry_two_files(
     $       filename1,
     $       filename2,
     $       analyze_grdpts_id,
     $       detailled)

          implicit none
          
          character(len=1024), intent(in) :: filename1
          character(len=1024), intent(in) :: filename2
          logical            , intent(in) :: analyze_grdpts_id
          logical            , intent(in) :: detailled

          type(pmodel_eq)                            :: p_model

          real(rkind)                                :: time1
          real(rkind), dimension(:)    , allocatable :: x_map1
          real(rkind), dimension(:)    , allocatable :: y_map1
          real(rkind), dimension(:,:,:), allocatable :: nodes1
          integer    , dimension(:,:)  , allocatable :: grdpts_id1

          real(rkind)                                :: time2
          real(rkind), dimension(:)    , allocatable :: x_map2
          real(rkind), dimension(:)    , allocatable :: y_map2
          real(rkind), dimension(:,:,:), allocatable :: nodes2
          integer    , dimension(:,:)  , allocatable :: grdpts_id2

          
          !extract netcdf data
          call extract_data(
     $         filename1,
     $         analyze_grdpts_id,
     $         time1,
     $         x_map1,
     $         y_map1,
     $         nodes1,
     $         grdpts_id1)

          call extract_data(
     $         filename2,
     $         analyze_grdpts_id,
     $         time2,
     $         x_map2,
     $         y_map2,
     $         nodes2,
     $         grdpts_id2)

          !check whether symmetry is respected along the y-axis
          if(analyze_grdpts_id) then
             call check_y_symmetry_two_files(
     $            p_model,
     $            grdpts_id1,grdpts_id2,
     $            nodes1,nodes2,
     $            detailled)

          else
             print '(''analyze_grdpts_id should be true'')'
             print '(''when analyzing two files...'')'
          end if

          !deallocate tables
          deallocate(x_map1)
          deallocate(y_map1)
          deallocate(nodes1)
          if(analyze_grdpts_id) then
             deallocate(grdpts_id1)
          end if

          deallocate(x_map2)
          deallocate(y_map2)
          deallocate(nodes2)
          if(analyze_grdpts_id) then
             deallocate(grdpts_id2)
          end if

        end subroutine verify_field_y_symmetry_two_files


        subroutine check_y_symmetry(
     $       p_model,
     $       nodes,
     $       detailled,
     $       grdpts_id)

          implicit none

          type(pmodel_eq)                        , intent(in) :: p_model
          real(rkind), dimension(:,:,:)          , intent(in) :: nodes
          logical                                , intent(in) :: detailled
          integer    , dimension(:,:)  , optional, intent(in) :: grdpts_id

          integer(ikind) :: i,j
          integer        :: k
          integer        :: size_x
          integer        :: size_y
          integer        :: j_sym

          integer, dimension(ne) :: var_type
          integer                :: sign_y
          logical                :: symmetric_grdpts_id
          logical                :: symmetric_all_field
          logical                :: symmetric_var
          logical                :: symmetric_loc
          
          var_type = p_model%get_var_type()
          size_x   = size(nodes,1)
          size_y   = size(nodes,2)

          symmetric_all_field = .true.

          if(present(grdpts_id)) then
             symmetric_grdpts_id=.true.

             do j=1, int(size_y/2.0)

                j_sym = size_y-j+1
                
                do i=1, size_x

                   symmetric_loc =
     $                  grdpts_id(i,j).eq.grdpts_id(i,j_sym)
                   symmetric_grdpts_id =
     $                  symmetric_grdpts_id.and.symmetric_loc

                   if(detailled.and.(.not.symmetric_loc)) then

                      print '(''['',2I4'']->['',2I4,'']: '',I2,''->'',I2)',
     $                     i,j,
     $                     i,j_sym,
     $                     grdpts_id(i,j),
     $                     grdpts_id(i,j_sym)

                   end if 

                end do

             end do

             print '(''symmetric grdpts_id?: '',L1)',
     $            symmetric_grdpts_id
             print '()'


             do k=1, ne
                
                if(var_type(k).eq.vector_y) then
                   sign_y = -1
                else
                   sign_y = +1
                end if

                symmetric_var = .true.

                do j=1,int(size_y/2.0)
                   
                   j_sym = size_y-j+1

                   do i=1, size(nodes,1)

                      if(grdpts_id(i,j).ne.no_pt) then

                         symmetric_loc = is_test_validated(
     $                        nodes(i,j,k),
     $                        sign_y*nodes(i,j_sym,k),
     $                        .false.)

                      else

                         symmetric_loc=.true.

                      end if

                      symmetric_var =
     $                     symmetric_var.and.symmetric_loc
                      
                      if(detailled.and.(.not.symmetric_loc)) then
                         
                         print '(''['',3I4'']->['',3I4,'']: '',F19.14,''->'',F19.14)',
     $                        i,j,k,
     $                        i,j_sym,k,
     $                        nodes(i,j,k),
     $                        sign_y*nodes(i,j_sym,k)
                         
                      end if
                         
                   end do

                end do

                print '(''var('',I1,'') symmetric ?: '',L1)', k, symmetric_var
                print '()'
                symmetric_all_field = symmetric_all_field.and.symmetric_var

             end do


          else

             do k=1, ne
                
                if(var_type(k).eq.vector_y) then
                   sign_y = -1
                else
                   sign_y = +1
                end if

                symmetric_var = .true.

                do j=1,int(size_y/2.0)
                   
                   j_sym = size_y-j+1

                   do i=1, size(nodes,1)

                      symmetric_loc = is_test_validated(
     $                     nodes(i,j,k),
     $                     sign_y*nodes(i,j_sym,k),
     $                     .false.)
                      symmetric_var =
     $                     symmetric_var.and.symmetric_loc

                      if(detailled.and.(.not.symmetric_loc)) then

                         print '(''['',3I4'']->['',3I4,'']: '',F19.14,''->'',F19.14)',
     $                        i,j,k,
     $                        i,j_sym,k,
     $                        nodes(i,j,k),
     $                        sign_y*nodes(i,j_sym,k)

                      end if                   

                   end do

                end do

                print '(''var('',I1,'') symmetric ?: '',L1)', k, symmetric_var
                print '()'
                symmetric_all_field = symmetric_all_field.and.symmetric_var

             end do

          end if

          print '(''field symmetric ?: '',L1)', symmetric_all_field
          print '()'

        end subroutine check_y_symmetry


        subroutine check_y_symmetry_two_files(
     $     p_model,
     $     grdpts_id1,grdpts_id2,
     $     nodes1,nodes2,
     $     detailled)

          implicit none

          type(pmodel_eq)              , intent(in) :: p_model
          integer    , dimension(:,:)  , intent(in) :: grdpts_id1
          integer    , dimension(:,:)  , intent(in) :: grdpts_id2
          real(rkind), dimension(:,:,:), intent(in) :: nodes1
          real(rkind), dimension(:,:,:), intent(in) :: nodes2
          logical                      , intent(in) :: detailled


          integer(ikind) :: i,j
          integer        :: k
          integer        :: size_x
          integer        :: size_y
          integer        :: j_sym

          integer, dimension(ne) :: var_type
          integer                :: sign_y
          logical                :: symmetric_grdpts_id
          logical                :: symmetric_all_field
          logical                :: symmetric_var
          logical                :: symmetric_loc
          

          var_type = p_model%get_var_type()
          size_x   = size(nodes1,1)
          size_y   = size(nodes1,2)


          if(size_x.ne.size(nodes2,1)) then
             print '(''symmetry pb: '')'
             print '(''both nodes do not have the same nx'')'
             stop ''
          end if

          if(size_y.ne.size(nodes2,2)) then
             print '(''symmetry pb: '')'
             print '(''both nodes do not have the same ny'')'
             stop ''
          end if


          symmetric_all_field = .true.
          symmetric_grdpts_id=.true.


          do j=1, int(size_y)

             j_sym = size_y-j+1
             
             do i=1, size_x
                
                symmetric_loc =
     $               grdpts_id1(i,j).eq.grdpts_id2(i,j_sym)
                symmetric_grdpts_id =
     $               symmetric_grdpts_id.and.symmetric_loc
                
                if(detailled.and.(.not.symmetric_loc)) then
                   
                   print '(''['',2I4'']->['',2I4,'']: '',I2,''->'',I2)',
     $                  i,j,
     $                  i,j_sym,
     $                  grdpts_id1(i,j),
     $                  grdpts_id2(i,j_sym)

                end if 

             end do

          end do

          print '(''symmetric grdpts_id?: '',L1)',
     $         symmetric_grdpts_id
          print '()'


          do k=1, ne
             
             if(var_type(k).eq.vector_y) then
                sign_y = -1
             else
                sign_y = +1
             end if

             symmetric_var = .true.

             do j=1,int(size_y)
                
                j_sym = size_y-j+1

                do i=1, size_x

                   if(grdpts_id1(i,j).ne.no_pt) then

                      symmetric_loc = is_test_validated(
     $                     nodes1(i,j,k),
     $                     sign_y*nodes2(i,j_sym,k),
     $                     .false.)

                   else

                      symmetric_loc=.true.

                   end if

                   symmetric_var =
     $                  symmetric_var.and.symmetric_loc
                   
                   if(detailled.and.(.not.symmetric_loc)) then
                      
                      print '(''['',3I4'']->['',3I4,'']: '',F19.14,''->'',F19.14)',
     $                     i,j,k,
     $                     i,j_sym,k,
     $                     nodes1(i,j,k),
     $                     sign_y*nodes2(i,j_sym,k)
                      
                   end if
                   
                end do

             end do

             print '(''var('',I1,'') symmetric ?: '',L1)', k, symmetric_var
             print '()'
             symmetric_all_field = symmetric_all_field.and.symmetric_var

          end do

          print '(''field symmetric ?: '',L1)', symmetric_all_field
          print '()'

        end subroutine check_y_symmetry_two_files


        subroutine extract_data(
     $     filename,
     $     analyze_grdpts_id,
     $     time,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id)

          implicit none

          character*(*)                             , intent(in)  :: filename
          logical                                   , intent(in)  :: analyze_grdpts_id
          real(rkind)                               , intent(out) :: time
          real(rkind), dimension(:)    , allocatable, intent(out) :: x_map
          real(rkind), dimension(:)    , allocatable, intent(out) :: y_map
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: nodes
          integer    , dimension(:,:)  , allocatable, intent(out) :: grdpts_id


          type(pmodel_eq)        :: p_model
          integer                :: ncid
          integer, dimension(3)  :: coordinates_id
          integer, dimension(ne) :: data_id
          integer                :: grdptsid_id

          real(rkind), dimension(2) :: x_borders
          real(rkind), dimension(2) :: y_borders
          integer    , dimension(2) :: sizes


          !open netcdf file
          call nf90_open_file_for_reading(
     $         trim(filename),
     $         ncid)

          !get the identifiers for the variables inside
          if(analyze_grdpts_id) then
             call nf90_get_varid(
     $            ncid,
     $            p_model,
     $            coordinates_id,
     $            data_id,
     $            grdptsid_id=grdptsid_id)
          else
             call nf90_get_varid(
     $            ncid,
     $            p_model,
     $            coordinates_id,
     $            data_id)
          end if

          !read the x_map,y_map,nodes (and grdpts_id)
          if(analyze_grdpts_id) then

             call nf90_read_borders(
     $            ncid,
     $            coordinates_id, 
     $            x_borders,
     $            y_borders,
     $            sizes)

             allocate(x_map(sizes(1)))
             allocate(y_map(sizes(2)))
             allocate(nodes(sizes(1),sizes(2),ne))
             allocate(grdpts_id(sizes(1),sizes(2)))

             call nf90_get_var_model_nopt(
     $            ncid,
     $            coordinates_id,
     $            data_id,
     $            time,
     $            nodes,
     $            x_map,
     $            y_map,
     $            [1,sizes(1),sizes(2)],
     $            grdptsid_id=grdptsid_id,
     $            grdpts_id=grdpts_id)

          else
             
             allocate(x_map(nx))
             allocate(y_map(ny))
             allocate(nodes(nx,ny,ne))

             call nf90_get_var_model(
     $            ncid,
     $            coordinates_id,
     $            data_id,
     $            time,
     $            nodes,
     $            x_map,
     $            y_map)

          end if

          !close file
          call nf90_close_file(ncid)

        end subroutine extract_data

      end program test_verify_y_symmetry
