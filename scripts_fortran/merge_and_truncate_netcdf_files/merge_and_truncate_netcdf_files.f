      program merge_and_truncate_netcdf_files
      
        use cmd_operators_mt_class, only :
     $       cmd_operators_mt

        use file_merge_param_class, only :
     $       file_merge_param

        use nf90_operators_mt_module, only :
     $       nf90_get_coordinate_borders,
     $       nf90_get_lg_domain_param,
     $       get_rank_tr_corner_files,
     $       get_filename,
     $       get_corner_filenames,
     $       nf90_merge_files_at_timestep

        use parameters_kind, only :
     $       rkind

        use progress_mg_module, only :
     $       ini_progress_mg,
     $       update_progress_mg,
     $       finalize_progress_mg

        implicit none

        character(len=1024)         :: pr_folder
        character(len=1024)         :: tr_folder

        real(rkind), dimension(2,2) :: borders_tr_domain
        integer    , dimension(2)   :: nb_tiles
        integer                     :: bc_size
        integer                     :: nb_timesteps

        real(rkind), dimension(2,2) :: borders_lg_domain
        real(rkind), dimension(2)   :: map_spacing
        integer    , dimension(2)   :: nb_pts_per_tile
        integer    , dimension(2,2) :: corner_ranks

        type(file_merge_param), dimension(:,:), allocatable :: fm_params

        real(rkind), dimension(1)                  :: time
        real(rkind), dimension(:)    , allocatable :: x_map
        real(rkind), dimension(:)    , allocatable :: y_map
        real(rkind), dimension(:,:,:), allocatable :: nodes

        integer :: nx,ny
        integer :: t

        integer, parameter :: unit_mg_progress=6
        integer, parameter :: ne=4
        
        ! analyse the command line arguments
        call get_command_line_arguments(
     $       pr_folder,
     $       tr_folder,
     $       borders_tr_domain,
     $       nb_tiles,
     $       bc_size,
     $       nb_timesteps)

        ! determine the characteristic of the large domain
        ! and how it should be processed to create the merge
        ! truncated domain
        call get_lg_domain_param(
     $       pr_folder,
     $       nb_tiles,
     $       borders_tr_domain,
     $       borders_lg_domain,
     $       map_spacing,
     $       nb_pts_per_tile,
     $       corner_ranks)

        ! determine the characteristics of the files that should
        ! be merged
        call get_file_merged_parameters(
     $       nb_tiles,
     $       borders_lg_domain,
     $       borders_tr_domain,
     $       map_spacing,
     $       nb_pts_per_tile,
     $       bc_size,
     $       fm_params)

        ! merge the files from t=0,nb_timesteps
        nx = nint((borders_tr_domain(1,2)-borders_tr_domain(1,1))/map_spacing(1))+1
        ny = nint((borders_tr_domain(2,2)-borders_tr_domain(2,1))/map_spacing(2))+1

        print '(''nx: '',I4)', nx
        print '(''ny: '',I4)', ny

        call ini_progress_mg(unit_mg_progress)

        allocate(x_map(nb_pts_per_tile(1)))
        allocate(y_map(nb_pts_per_tile(2)))
        allocate(nodes(nb_pts_per_tile(1),nb_pts_per_tile(2),ne))

        do t=0,nb_timesteps-1

           call nf90_merge_files_at_timestep(
     $          pr_folder,
     $          tr_folder,
     $          nx,ny,
     $          t,
     $          fm_params,
     $          time,
     $          x_map,
     $          y_map,
     $          nodes)

           call update_progress_mg(unit_mg_progress,t+1,nb_timesteps)

        end do

        deallocate(x_map)
        deallocate(y_map)
        deallocate(nodes)

        call finalize_progress_mg(unit_mg_progress)


        contains

        subroutine get_file_merged_parameters(
     $       nb_tiles,
     $       borders_lg_domain,
     $       borders_tr_domain,
     $       map_spacing,
     $       nb_pts_per_tile,
     $       bc_size,
     $       fm_params)
        
          implicit none

          integer               , dimension(2)               , intent(in)  :: nb_tiles
          real(rkind)           , dimension(2,2)             , intent(in)  :: borders_lg_domain
          real(rkind)           , dimension(2,2)             , intent(in)  :: borders_tr_domain
          real(rkind)           , dimension(2)               , intent(in)  :: map_spacing
          integer               , dimension(2)               , intent(in)  :: nb_pts_per_tile
          integer                                            , intent(in)  :: bc_size
          type(file_merge_param), dimension(:,:), allocatable, intent(out) :: fm_params


          integer :: i,j,i_s,j_s

          print '(''extraction_indices'')'
          print '(''----------------------------------------'')'

          allocate(fm_params(
     $         corner_ranks(1,2)-corner_ranks(1,1)+1,
     $         corner_ranks(2,2)-corner_ranks(2,1)+1))

          do j=corner_ranks(2,1),corner_ranks(2,2)

             j_s = j-corner_ranks(2,1)+1

             do i=corner_ranks(1,1),corner_ranks(1,2)

                i_s = i-corner_ranks(1,1)+1
                
                call fm_params(i_s,j_s)%set_ranks(
     $               i,j,
     $               nb_tiles(2))
                
                call fm_params(i_s,j_s)%determine_indices_for_extracting_and_saving(
     $               borders_lg_domain,
     $               borders_tr_domain,
     $               map_spacing,
     $               nb_pts_per_tile,
     $               bc_size)

                print '(''[i_min,i_max]-x: '',2I4)', fm_params(i_s,j_s)%extracted_from_file(1,:)
                print '(''[j_min,j_max]-y: '',2I4)', fm_params(i_s,j_s)%extracted_from_file(2,:)
                print '(''[i_min,i_max]-x: '',2I4)', fm_params(i_s,j_s)%saved_in_lg_domain(1,:)
                print '(''[j_min,j_max]-y: '',2I4)', fm_params(i_s,j_s)%saved_in_lg_domain(2,:)
                print '()'

             end do
          end do

          print '(''----------------------------------------'')'
          print '()'


        end subroutine get_file_merged_parameters

        subroutine get_lg_domain_param(
     $       pr_folder,
     $       nb_tiles,
     $       borders_tr_domain,
     $       borders_lg_domain,
     $       map_spacing,
     $       nb_pts_per_tile,
     $       corner_ranks)

          implicit none

          character*(*)              , intent(in)  :: pr_folder
          integer    , dimension(2)  , intent(in)  :: nb_tiles
          real(rkind), dimension(2,2), intent(in)  :: borders_tr_domain
          real(rkind), dimension(2,2), intent(out) :: borders_lg_domain
          real(rkind), dimension(2)  , intent(out) :: map_spacing
          integer    , dimension(2)  , intent(out) :: nb_pts_per_tile
          integer    , dimension(2,2), intent(out) :: corner_ranks

          
          character(len=16) :: SW_corner_filename
          character(len=16) :: NE_corner_filename

          character(len=1024) :: SW_path
          character(len=1024) :: NE_path

          integer :: len_pr_folder
          integer :: len_SW_file
          integer :: len_NE_file


          call get_corner_filenames(
     $         nb_tiles,
     $         SW_corner_filename,
     $         NE_corner_filename)

          ! get the path for the parallel files
          len_pr_folder = len(trim(pr_folder))
          len_SW_file   = len(trim(SW_corner_filename))
          len_NE_file   = len(trim(NE_corner_filename))

          SW_path(1:len_pr_folder) = trim(pr_folder)
          SW_path(len_pr_folder+1:len_pr_folder+len_SW_file) = trim(SW_corner_filename)

          NE_path(1:len_pr_folder) = trim(pr_folder)
          NE_path(len_pr_folder+1:len_pr_folder+len_NE_file) = trim(NE_corner_filename)


          call nf90_get_lg_domain_param(
     $         trim(SW_path),
     $         trim(NE_path),
     $         borders_lg_domain,
     $         map_spacing,
     $         nb_pts_per_tile)
          
          call get_rank_tr_corner_files(
     $         borders_lg_domain,
     $         borders_tr_domain,
     $         map_spacing,
     $         nb_tiles,
     $         nb_pts_per_tile,
     $         bc_size,
     $         corner_ranks)

          print '(''large domain parameters'')'
          print '(''--------------------------------------------------'')'
          print '(''[x_min_lg,x_max_lg]: '',2F8.2)', borders_lg_domain(1,:)
          print '(''[y_min_lg,y_max_lg]: '',2F8.2)', borders_lg_domain(2,:)
          print '(''[dx,dy]            : '',2F8.4)', map_spacing
          print '(''corner_ranks_x     : '',2I4)', corner_ranks(1,:)
          print '(''corner_ranks_y     : '',2I4)', corner_ranks(2,:)
          print '(''--------------------------------------------------'')'
          print '()'

        end subroutine get_lg_domain_param


        subroutine get_command_line_arguments(
     $     pr_folder,
     $     tr_folder,
     $     borders_tr_domain,
     $     nb_tiles,
     $     bc_size,
     $     nb_timesteps)

          implicit none

          character(len=1024)        , intent(out) :: pr_folder
          character(len=1024)        , intent(out) :: tr_folder
          real(rkind), dimension(2,2), intent(out) :: borders_tr_domain
          integer    , dimension(2)  , intent(out) :: nb_tiles
          integer                    , intent(out) :: bc_size
          integer                    , intent(out) :: nb_timesteps

          type(cmd_operators_mt)    :: cmd_op

          call cmd_op%analyse_cmd_line_arg()

          pr_folder = cmd_op%get_parallel_folder()
          tr_folder = cmd_op%get_truncation_folder()

          borders_tr_domain  = cmd_op%get_borders()
        
          nb_tiles(1)  = cmd_op%get_nb_tiles_x()
          nb_tiles(2)  = cmd_op%get_nb_tiles_y()
          
          bc_size      = cmd_op%get_bc_size()
          nb_timesteps = cmd_op%get_nb_timesteps()

          print '()'
          print '(''input parameters'')'
          print '(''--------------------------------------------------'')'
          print '(''[x_min,x_max]: '',2F9.4)', borders_tr_domain(1,:)
          print '(''[y_min,y_max]: '',2F9.4)', borders_tr_domain(2,:)
          print '('' nb_tiles    : '',2I2)'  , nb_tiles
          print '('' bc_size     : '',2I2)'  , bc_size
          print '('' nb_timesteps: '',I4)'   , nb_timesteps
          print '('' input folder: '',A20)'  , pr_folder
          print '('' output folder: '',A20)' , tr_folder
          print '(''--------------------------------------------------'')'
          print '()'

        end subroutine get_command_line_arguments

      end program merge_and_truncate_netcdf_files
