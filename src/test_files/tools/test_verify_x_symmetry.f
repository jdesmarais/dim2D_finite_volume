      program test_verify_x_symmetry

       use check_data_module, only :
     $       is_real_validated

       use cmd_operators_class, only :
     $       cmd_operators

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
     $       vector_x
        
        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none

        logical, parameter :: detailled = .true.        

        type(cmd_operators)   :: cmd_operators_used
        character(len=1024)   :: filename
        logical               :: analyse_grdpts_id
        integer, dimension(4) :: nb_bf_layers

        !analyse command line arguments
        call cmd_operators_used%analyse_cmd_line_arg()

        if(cmd_operators_used%is_restart_activated()) then
           filename = cmd_operators_used%get_restart_filename()

        else
           print '(''no file to be analyzed'')'
           stop ''

        end if

        nb_bf_layers = cmd_operators_used%get_nb_bf_layers()
        analyse_grdpts_id = nb_bf_layers(1).ne.0

        
        !verify the y-symmetry of the field
        call verify_field_x_symmetry(
     $       filename,
     $       analyse_grdpts_id,
     $       detailled)

        contains


        !check the y-symmetry in the field
        subroutine verify_field_x_symmetry(
     $       filename,
     $       analyse_grdpts_id,
     $       detailled)

          implicit none
          
          character(len=1024), intent(in) :: filename
          logical            , intent(in) :: analyse_grdpts_id
          logical            , intent(in) :: detailled

          type(pmodel_eq)        :: p_model
          integer                :: ncid
          integer, dimension(3)  :: coordinates_id
          integer, dimension(ne) :: data_id
          integer                :: grdptsid_id

          real(rkind)                                :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes
          integer    , dimension(:,:)  , allocatable :: grdpts_id

          real(rkind), dimension(2) :: x_borders
          real(rkind), dimension(2) :: y_borders
          integer    , dimension(2) :: sizes


          !open netcdf file
          call nf90_open_file_for_reading(
     $         trim(filename),
     $         ncid)

          !get the identifiers for the variables inside
          if(analyse_grdpts_id) then
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
          if(analyse_grdpts_id) then

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

          !check whether symmetry is respected along the y-axis
          if(analyse_grdpts_id) then
             call check_x_symmetry(
     $            p_model,
     $            nodes,
     $            detailled,
     $            grdpts_id=grdpts_id)
          else
             call check_x_symmetry(
     $            p_model,
     $            nodes,detailled)
          end if

          !deallocate tables
          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)
          if(analyse_grdpts_id) then
             deallocate(grdpts_id)
          end if

        end subroutine verify_field_x_symmetry


        subroutine check_x_symmetry(
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
          integer        :: i_sym

          integer, dimension(ne) :: var_type
          integer                :: sign_x
          logical                :: symmetric_grdpts_id
          logical                :: symmetric_all_field
          logical                :: symmetric_var
          logical                :: symmetric_loc

          logical                      :: symmetry_loss
          real(rkind)                  :: diff
          real(rkind)                  :: diff_loc
          integer(ikind), dimension(2) :: diff_indices

          symmetry_loss = .false.
          
          var_type = p_model%get_var_type()
          size_x   = size(nodes,1)
          size_y   = size(nodes,2)

          symmetric_all_field = .true.

          if(present(grdpts_id)) then
             symmetric_grdpts_id=.true.

             do j=1, size_y

                do i=1, int(size_x/2.0)

                   i_sym = size_x-i+1

                   symmetric_loc =
     $                  grdpts_id(i,j).eq.grdpts_id(i_sym,j)
                   symmetric_grdpts_id =
     $                  symmetric_grdpts_id.and.symmetric_loc

                   if(detailled.and.(.not.symmetric_loc)) then

                      print '(''['',2I4'']: '',I2,''->'',I2)',
     $                     i,j,
     $                     grdpts_id(i,j),
     $                     grdpts_id(i_sym,j)

                   end if 

                end do

             end do

             print '(''symmetric grdpts_id?: '',L1)',
     $            symmetric_grdpts_id
             print '()'

          else
             symmetric_grdpts_id = .true.
          end if

          if(symmetric_grdpts_id) then

             do k=1, ne
                
                if(var_type(k).eq.vector_x) then
                   sign_x = -1
                else
                   sign_x = +1
                end if

                symmetric_var = .true.

                do j=1,size_y
                   
                   do i=1, int(size_x/2.0)

                      i_sym = size_x-i+1


                      if(present(grdpts_id)) then
                         if(grdpts_id(i,j).ne.no_pt) then
                         
                            symmetric_loc = is_real_validated(
     $                           nodes(i,j,k),
     $                           sign_x*nodes(i_sym,j,k),
     $                           .true.)

                         else

                            symmetric_loc = .true.

                         end if
                      else

                         symmetric_loc = is_real_validated(
     $                           nodes(i,j,k),
     $                           sign_x*nodes(i_sym,j,k),
     $                           .true.)

                      end if
                      if(.not.symmetric_loc) then
                         diff_loc = abs(nodes(i,j,k) - sign_x*nodes(i_sym,j,k))
                         if(symmetry_loss) then
                            if(diff_loc.gt.diff) then
                               diff = diff_loc
                               diff_indices = [i,j]
                            end if
                         else
                            diff = diff_loc
                            symmetry_loss = .true.
                         end if
                      end if

                      symmetric_var =
     $                     symmetric_var.and.symmetric_loc

                      if(detailled.and.(.not.symmetric_loc)) then
                         
                         print '(''['',3I4'']: '',F16.13,''->'',F16.13)',
     $                        i,j,k,
     $                        nodes(i,j,k),
     $                        sign_x*nodes(i_sym,j,k)
                         
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

          if(symmetry_loss) then
             print *, 'largest difference', diff
             print '(''largest difference at : '',2I4)', diff_indices
          end if

        end subroutine check_x_symmetry

      end program test_verify_x_symmetry
