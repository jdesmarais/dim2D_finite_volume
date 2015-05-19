      module nf90_operators_mt_module

        use file_merge_param_class, only :
     $     file_merge_param

        use netcdf

        use parameters_kind, only :
     $       rkind

        use parameters_cst, only :
     $       SUCCESS,
     $       NOT_SUCCESS

        implicit none

        private
        public ::
     $       nf90_get_coordinate_borders,
     $       nf90_get_lg_domain_param,
     $       get_rank_tr_corner_files,
     $       get_filename,
     $       get_corner_filenames,
     $       nf90_def_header_merge,
     $       nf90_def_var_merge,
     $       nf90_merge_files_at_timestep

        contains


        subroutine nf90_merge_files_at_timestep(
     $       pr_folder,
     $       tr_folder,
     $       nx,ny,
     $       timestep,
     $       fm_params,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes)

          implicit none

          character*(*)                           , intent(in)    :: pr_folder
          character*(*)                           , intent(in)    :: tr_folder
          integer                                 , intent(in)    :: nx
          integer                                 , intent(in)    :: ny
          integer                                 , intent(in)    :: timestep
          type(file_merge_param), dimension(:,:)  , intent(in)    :: fm_params
          real(rkind)           , dimension(1)    , intent(out)   :: time
          real(rkind)           , dimension(:)    , intent(out)   :: x_map
          real(rkind)           , dimension(:)    , intent(out)   :: y_map
          real(rkind)           , dimension(:,:,:), intent(out)   :: nodes

          integer             :: len_tr_folder
          integer             :: len_file_merge
          character(len=1024) :: path_merge
          character(len=16)   :: filename_merge
          integer             :: ncid_merge
          integer             :: ncid_tile

          integer, dimension(3) :: coordinates_id
          integer, dimension(4) :: var_id          

          integer, dimension(2,2) :: borders_tile
          integer, dimension(2,2) :: borders_merge

          integer :: rank_tile

          integer :: i,j


          !1) open the netcdf file where the different files
          !   are merged
          call get_filename(
     $         filename_merge,
     $         timestep)

          len_tr_folder  = len(trim(tr_folder))
          len_file_merge = len(trim(filename_merge))

          path_merge=''
          path_merge(1:len_tr_folder) = trim(tr_folder)
          path_merge(len_tr_folder+1:len_tr_folder+len_file_merge) = trim(filename_merge)

          call nf90_handle_err(
     $         NF90_CREATE(
     $         trim(path_merge),
     $         NF90_NETCDF4,
     $         ncid_merge))


          !2) open the first file to be merged inside to define
          !   the header and the main variables
          rank_tile = fm_params(1,1)%get_rank()

          call get_filename(
     $         filename_merge,
     $         timestep,
     $         rank=rank_tile)

          len_tr_folder  = len(trim(pr_folder))
          len_file_merge = len(trim(filename_merge))

          path_merge(1:len_tr_folder) = trim(pr_folder)
          path_merge(len_tr_folder+1:len_tr_folder+len_file_merge) = trim(filename_merge)

          call nf90_handle_err(
     $         NF90_OPEN(
     $         trim(path_merge),
     $         NF90_NOWRITE,
     $         ncid_tile))

          call nf90_def_header_merge(
     $         ncid_tile,
     $         ncid_merge)

          call nf90_def_var_merge(
     $         ncid_tile,
     $         ncid_merge,
     $         nx,ny,
     $         coordinates_id,
     $         var_id)

          call nf90_handle_err(
     $         NF90_ENDDEF(ncid_merge))

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_tile))


          !3) merge the files
          do j=1,size(fm_params,2)
             do i=1, size(fm_params,1)

                ! open file
                rank_tile = fm_params(i,j)%get_rank()

                call get_filename(
     $               filename_merge,
     $               timestep,
     $               rank=rank_tile)
                
                len_tr_folder  = len(trim(pr_folder))
                len_file_merge = len(trim(filename_merge))
                
                path_merge(1:len_tr_folder) = trim(pr_folder)
                path_merge(len_tr_folder+1:len_tr_folder+len_file_merge) = trim(filename_merge)
                
                call nf90_handle_err(
     $               NF90_OPEN(
     $               trim(path_merge),
     $               NF90_NOWRITE,
     $               ncid_tile))

                ! put var
                borders_tile  = fm_params(i,j)%get_indices_for_extracting()
                borders_merge = fm_params(i,j)%get_indices_for_saving()

                call nf90_put_var_merge(
     $               ncid_tile,
     $               ncid_merge,
     $               coordinates_id,
     $               var_id,
     $               time,
     $               x_map,
     $               y_map,
     $               nodes,
     $               borders_tile,
     $               borders_merge)

             end do
          end do

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_merge))

        end subroutine nf90_merge_files_at_timestep

        
        !extract the variable from the file and store them in the
        !merge file
        ! - from the rank of the processor which created the output
        !   file, determine where the extracted variables should be
        !   saved
        subroutine nf90_put_var_merge(
     $       ncid_tile,
     $       ncid_merge,
     $       coordinates_id,
     $       var_id,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes,
     $       borders_tile,
     $       borders_merge,
     $       ierror)

          implicit none

          integer                        , intent(in)    :: ncid_tile
          integer                        , intent(in)    :: ncid_merge
          integer      , dimension(3)    , intent(in)    :: coordinates_id
          integer      , dimension(:)    , intent(in)    :: var_id
          real(rkind)  , dimension(1)    , intent(inout) :: time
          real(rkind)  , dimension(:)    , intent(inout) :: x_map
          real(rkind)  , dimension(:)    , intent(inout) :: y_map
          real(rkind)  , dimension(:,:,:), intent(inout) :: nodes
          integer      , dimension(2,2)  , intent(in)    :: borders_tile
          integer      , dimension(2,2)  , intent(in)    :: borders_merge
          integer      , optional        , intent(out)   :: ierror


          integer :: size_x_tile
          integer :: size_y_tile

          integer, dimension(2) :: count
          integer, dimension(2) :: start
          integer, dimension(4) :: limit

          integer :: NF_MYREAL

          integer :: k
          integer :: ne


          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          size_x_tile = size(x_map,1)
          size_y_tile = size(y_map,1)
          ne          = size(nodes,3)


          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)          


          !2) deduce the indices where the variables
          !   extracted from the file will be saved
          !------------------------------------------------------------
          count(1) = borders_tile(1,2)-borders_tile(1,1)+1
          count(2) = borders_tile(2,2)-borders_tile(2,1)+1

          start(1) = borders_merge(1,1)
          start(2) = borders_merge(2,1)

          limit(1) = borders_tile(1,1)
          limit(2) = borders_tile(1,2)
          limit(3) = borders_tile(2,1)
          limit(4) = borders_tile(2,2)


          !3) extract the time, coordinates, and governing variables
          !   from the rank file
          !------------------------------------------------------------
          !3.1) get the time
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_tile,
     $         coordinates_id(1),
     $         time))
          
          !3.2) get the x_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_tile,
     $         coordinates_id(2),
     $         x_map))
          
          !3.3) get the y_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_tile,
     $         coordinates_id(3),
     $         y_map))             
      
          !3.4) get the governing variables
          do k=1,ne
          
             call nf90_handle_err(
     $            NF90_GET_VAR(
     $            ncid_tile,
     $            var_id(k),
     $            nodes(:,:,k),
     $            START=[1,1,1],
     $            COUNT=[1,size_x_tile,size_y_tile]))
          end do

                          
          print '(''extraction'')'
          print '(''----------------'')'
          print '(''[x_min,x_max]: '',2F7.4)', x_map(limit(1)), x_map(limit(2))
          print '(''[y_min,y_max]: '',2F7.4)', y_map(limit(3)), y_map(limit(4))


          !6) save the time, coordinates, and governing variables
          !   in the merge file
          !------------------------------------------------------------
          !6.1) write the time
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_merge,
     $         coordinates_id(1),
     $         time))

          !6.2) write the x-coordinates
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_merge,
     $         coordinates_id(2),
     $         x_map(limit(1):limit(2)),
     $         START=[start(1)],
     $         COUNT=[count(1)]))

          !6.3) write the y-coordinates
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_merge,
     $         coordinates_id(3),
     $         y_map(limit(3):limit(4)),
     $         START=[start(2)],
     $         COUNT=[count(2)]))

          !6.4) write the variables of the governing equations
          do k=1, ne

             call nf90_handle_err(
     $            NF90_PUT_VAR(
     $            ncid_merge,
     $            var_id(k),
     $            nodes(limit(1):limit(2),limit(3):limit(4),k),
     $            START=[1, start(1), start(2)], 
     $            COUNT=[1, count(1), count(2)]))
          end do


          !7) close the processor file
          call nf90_handle_err(
     $         NF90_CLOSE(ncid_tile))

          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end subroutine nf90_put_var_merge


        ! define the variables saved in the file
        ! for the truncated domain knowing the
        ! truncation indices
        subroutine nf90_def_var_merge(
     $       ncid_tile,
     $       ncid_merge,
     $       nx,ny,
     $       coordinates_id,
     $       var_id)

          implicit none

          integer                       , intent(in)  :: ncid_tile
          integer                       , intent(in)  :: ncid_merge
          integer                       , intent(in)  :: nx
          integer                       , intent(in)  :: ny
          integer       , dimension(3)  , intent(out) :: coordinates_id
          integer       , dimension(:)  , intent(out) :: var_id

          integer :: retval

          integer, dimension(3) :: sizes
          integer               :: NF_MYREAL
          integer, dimension(3) :: dims_id
          integer               :: k,l

          character(len=60)   :: nameVar

          integer             :: nAttsVar
          character(len=60)   :: nameAtt
          character(len=1024) :: charVarAtt


          !0) sizes
          sizes = [1,nx,ny]

          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)

          
          !3) extract the name, units and longname
          !   of the coordinates
          do k=1,3
             
             ! extract the name of the coordinate
             ! and save it in the error_file
             retval = NF90_INQUIRE_VARIABLE(
     $            ncid_tile,
     $            k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)
             

             ! define the dimension with the name
             retval = NF90_DEF_DIM(
     $            ncid_merge,
     $            trim(nameVar),
     $            sizes(k),
     $            dims_id(k))
             call nf90_handle_err(retval)


             ! define the corresponding variable
             retval = NF90_DEF_VAR(
     $            ncid_merge,
     $            trim(nameVar),
     $            NF_MYREAL,
     $            dims_id(k),
     $            coordinates_id(k))
             call nf90_handle_err(retval)


             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                retval = NF90_INQ_ATTNAME(
     $               ncid_tile,
     $               k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_tile,
     $               k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                retval = NF90_PUT_ATT(
     $               ncid_merge,
     $               coordinates_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt))
                call nf90_handle_err(retval)
                
             end do

          end do


          !4) extract and define the governing 
          !   variables
          !   [1,3]: coordinates
          !   [4,:]: governing variables
          do k=1, size(var_id,1)

             ! extract the governing variable
             ! from the netcdf file of the
             ! small domain
             retval = NF90_INQUIRE_VARIABLE(
     $            ncid_tile,
     $            3+k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)

             ! define the governing variable to
             ! be saved in the netcdf file 
             ! containing the truncation
             retval = NF90_DEF_VAR(
     $            ncid_merge,
     $            nameVar,
     $            NF_MYREAL,
     $            dims_id,
     $            var_id(k))
             call nf90_handle_err(retval)

             
             ! extract and assign the units to
             ! the governing variables
             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                retval = NF90_INQ_ATTNAME(
     $               ncid_tile,
     $               3+k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_tile,
     $               3+k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                retval = NF90_PUT_ATT(
     $               ncid_merge,
     $               var_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt))
                call nf90_handle_err(retval)
                
             end do

          end do

        end subroutine nf90_def_var_merge


        ! get the NF90_TYPE corresponding to real(rkind)
        function nf90_get_myreal(rkind)
     $     result(NF_MYREAL)

          implicit none

          integer, intent(in) :: rkind
          integer             :: NF_MYREAL

          select case(rkind)
            case(4)
               NF_MYREAL=NF90_FLOAT
            case(8)
               NF_MYREAL=NF90_DOUBLE
            case default
               print '(''nf90_operators_error_module'')'
               print '(''get_nf90_myreal'')'
               stop 'NF_MYREAL'
          end select

        end function nf90_get_myreal


        !define the header of the merge file
        ! - create the file for merging the netcdf files
        ! - read the rank 0 file and copy the header in
        !   the merge file
        ! - define the coordinates and the variables ID
        !   by computing the total size of the domain
        !------------------------------------------------------------
        subroutine nf90_def_header_merge(
     $       ncid_tile,
     $       ncid_merge,
     $       ierror)

          implicit none


          integer                , intent(in)  :: ncid_tile
          integer                , intent(in)  :: ncid_merge
          integer      , optional, intent(out) :: ierror

          integer               :: retval

          integer               :: size_header_tile

          integer               :: k

          character(len=1028)   :: header_attname_tile
          integer               :: header_att_type
          integer               :: header_att_len

          character(len=1028)   :: att_char_var
          integer               :: att_int_var
          real                  :: att_real_var
          real*8                :: att_double_var

        
          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          !1) inquire the number of global
          !   variables (header variables)
          retval = NF90_INQUIRE(
     $         ncid_tile,
     $         nAttributes=size_header_tile)
          call nf90_handle_err(retval)


          !2) inquire the name of the global
          !   variables for the small and large
          !   domains
          do k=1, size_header_tile

             ! inquire the name of the global
             ! variables in the header
             retval = NF90_INQ_ATTNAME(
     $            ncid_tile,
     $            NF90_GLOBAL,
     $            k,
     $            header_attname_tile)
             call nf90_handle_err(retval)


             ! inquire the type and length of
             ! the global variable
             retval = NF90_INQUIRE_ATTRIBUTE(
     $            ncid_tile,
     $            NF90_GLOBAL,
     $            trim(header_attname_tile),
     $            xtype=header_att_type,
     $            len=header_att_len)
             call nf90_handle_err(retval)


             ! inquire the value of the global
             ! variable
             call nf90_get_att_value(
     $            ncid_tile,
     $            NF90_GLOBAL,
     $            header_att_type,
     $            trim(header_attname_tile),
     $            att_char_var,
     $            att_int_var,
     $            att_real_var,
     $            att_double_var)

             
             ! put the value in the merge file
             call nf90_put_att_value(
     $            ncid_merge,
     $            NF90_GLOBAL,
     $            header_att_type,
     $            trim(header_attname_tile),
     $            att_char_var,
     $            att_int_var,
     $            att_real_var,
     $            att_double_var)

          end do


          if(present(ierror)) then
             ierror = SUCCESS
          end if

      end subroutine nf90_def_header_merge


      ! get the value corresponding to the attribute with
      ! the name 'attname'
      subroutine nf90_get_att_value(
     $   ncid,
     $   varid,
     $   att_type,
     $   att_name,
     $   att_char_var,
     $   att_int_var,
     $   att_real_var,
     $   att_double_var)

        implicit none

        integer            , intent(in) :: ncid
        integer            , intent(in) :: varid
        integer            , intent(in) :: att_type
        character*(*)      , intent(in) :: att_name
        character(len=1028), intent(out):: att_char_var
        integer            , intent(out):: att_int_var
        real               , intent(out):: att_real_var
        real*8             , intent(out):: att_double_var

        
        integer :: retval


        select case(att_type)

          case(NF90_CHAR)
             retval = nf90_get_att(
     $            ncid,
     $            varid,
     $            trim(att_name),
     $            att_char_var)

          case(NF90_INT)
             retval = nf90_get_att(
     $            ncid,
     $            varid,
     $            trim(att_name),
     $            att_int_var)

          case(NF90_FLOAT)
             retval = nf90_get_att(
     $            ncid,
     $            varid,
     $            trim(att_name),
     $            att_real_var)

          case(NF90_DOUBLE)
             retval = nf90_get_att(
     $            ncid,
     $            varid,
     $            trim(att_name),
     $            att_double_var)
           
          case default
             print '(''nf90_operators_error_module'')'
             print '(''nf90_get_att_var'')'
             print '(''var type not recognized'')'
             stop ''
             
          end select

          call nf90_handle_err(retval)

      end subroutine nf90_get_att_value


      ! put the value corresponding to the attribute with
      ! the name 'attname'
      subroutine nf90_put_att_value(
     $     ncid,
     $     varid,
     $     att_type,
     $     att_name,
     $     att_char_var,
     $     att_int_var,
     $     att_real_var,
     $     att_double_var)

          implicit none

          integer      , intent(in) :: ncid
          integer      , intent(in) :: varid
          integer      , intent(in) :: att_type
          character*(*), intent(in) :: att_name
          character*(*), intent(in) :: att_char_var
          integer      , intent(in) :: att_int_var
          real         , intent(in) :: att_real_var
          real*8       , intent(in) :: att_double_var


          integer :: retval


          select case(att_type)

            case(NF90_INT)
               retval = nf90_put_att(
     $              ncid,
     $              varid,
     $              trim(att_name),
     $              att_int_var)
               
            case(NF90_CHAR)
               retval = nf90_put_att(
     $              ncid,
     $              varid,
     $              trim(att_name),
     $              att_char_var)

            case(NF90_FLOAT)
               retval = nf90_put_att(
     $              ncid,
     $              varid,
     $              trim(att_name),
     $              att_real_var)

            case(NF90_DOUBLE)
               retval = nf90_put_att(
     $              ncid,
     $              varid,
     $              trim(att_name),
     $              att_double_var)
                       
          end select

          call nf90_handle_err(retval)

        end subroutine nf90_put_att_value


        subroutine get_corner_filenames(
     $       nb_tiles,
     $       SW_corner_filename,
     $       NE_corner_filename)

          implicit none

          integer, dimension(2), intent(in)  :: nb_tiles
          character(len=16)    , intent(out) :: SW_corner_filename
          character(len=16)    , intent(out) :: NE_corner_filename


          call get_filename(SW_corner_filename,0,0)
          call get_filename(NE_corner_filename,0,nb_tiles(1)*nb_tiles(2)-1)


        end subroutine get_corner_filenames


        subroutine get_filename(filename, timestep, rank)

          implicit none

          integer          , intent(in)  :: timestep
          character(len=16), intent(out) :: filename
          integer, optional, intent(in)  :: rank

          integer           :: size_timestep
          integer           :: size_rank
          character(len=16) :: format_string


          ! determine the character size to write the timestep
          if (timestep.eq.0) then
             size_timestep  = 1
          else
             size_timestep  = floor(log(float(timestep))/log(10.))+1
          end if

          if(present(rank)) then

             ! determine the character size to write the timestep
             if (rank.eq.0) then
                size_rank  = 1
             else
                size_rank  = floor(log(float(rank))/log(10.))+1
             end if

             !create the format string
             write (format_string, "(A5,I1,A5,I1,A4)")
     $            '(A4,I', size_timestep,
     $            ',A1,I', size_rank,
     $            ',A3)'

             !create the name of the file
             write (filename, trim(format_string)) 'data', timestep, '_', rank, '.nc'

          else

             !create the format string
             write (format_string, "(A5,I1,A4)")
     $            '(A4,I', size_timestep, ',A3)'

             !create the name of the file
             write (filename, trim(format_string)) 'data', timestep, '.nc'

          end if

        end subroutine get_filename


        subroutine get_rank_tr_corner_files(
     $       borders_lg_domain,
     $       borders_tr_domain,
     $       map_spacing,
     $       nb_tiles,
     $       nb_pts_per_tile,
     $       bc_size,
     $       corner_ranks)

          implicit none

          real(rkind), dimension(2,2), intent(in)  :: borders_lg_domain
          real(rkind), dimension(2,2), intent(in)  :: borders_tr_domain
          real(rkind), dimension(2)  , intent(in)  :: map_spacing
          integer    , dimension(2)  , intent(in)  :: nb_tiles
          integer    , dimension(2)  , intent(in)  :: nb_pts_per_tile
          integer                    , intent(in)  :: bc_size
          integer    , dimension(2,2), intent(out) :: corner_ranks
          
          integer               :: i_corner
          integer               :: i_tile
          integer               :: j,k
          integer, dimension(2) :: n
          

          do j=1,2
             n(j) = nint((borders_lg_domain(j,2)-borders_lg_domain(j,1))/
     $            map_spacing(j)) + 1
          end do

          do k=1,2
             do j=1,2

                i_corner = nint((borders_tr_domain(j,k)-borders_lg_domain(j,1))/
     $               map_spacing(j)) + 1

                

                if(i_corner.le.bc_size) then
                   corner_ranks(j,k) = 0

                else

                   if(i_corner.ge.(n(j)-bc_size+1)) then
                      corner_ranks(j,k) = nb_tiles(j)-1

                   else
                      i_corner = i_corner-bc_size-1
                      i_tile   = mod(i_corner,nb_pts_per_tile(j)-2*bc_size)

                      !print *, 'i_corner: ', j,k, i_corner
                      !print *, 'i_tile:   ', j,k, i_tile
                      !print *, ''

                      corner_ranks(j,k) = (i_corner-i_tile)/(nb_pts_per_tile(j)-2*bc_size)

                   end if

                end if
                
             end do
          end do

        end subroutine get_rank_tr_corner_files


        subroutine nf90_get_lg_domain_param(
     $       SW_corner_filename,
     $       NE_corner_filename,
     $       borders_lg_domain,
     $       map_spacing,
     $       nb_pts_per_tile)

          implicit none

          character*(*)              , intent(in)  :: SW_corner_filename
          character*(*)              , intent(in)  :: NE_corner_filename
          real(rkind), dimension(2,2), intent(out) :: borders_lg_domain
          real(rkind), dimension(2)  , intent(out) :: map_spacing
          integer    , dimension(2)  , intent(out) :: nb_pts_per_tile

          
          real(rkind), dimension(2,2) :: SW_map_borders
          real(rkind), dimension(2)   :: SW_map_spacing
          real(rkind), dimension(2,2) :: NE_map_borders
          real(rkind), dimension(2)   :: NE_map_spacing

          integer :: j


          call nf90_get_coordinate_borders(
     $         SW_corner_filename,
     $         SW_map_borders,
     $         SW_map_spacing)

          call nf90_get_coordinate_borders(
     $         NE_corner_filename,
     $         NE_map_borders,
     $         NE_map_spacing)


          map_spacing = SW_map_spacing

          do j=1,2

             borders_lg_domain(j,1) = SW_map_borders(j,1)
             borders_lg_domain(j,2) = NE_map_borders(j,2)

             nb_pts_per_tile(j) = nint(
     $            (SW_map_borders(j,2)-SW_map_borders(j,1))/
     $            map_spacing(j)) + 1

          end do

        end subroutine nf90_get_lg_domain_param


        subroutine nf90_get_coordinate_borders(
     $       filename,
     $       map_borders,map_spacing)

          implicit none
          
          character*(*)              , intent(in)  :: filename
          real(rkind), dimension(2,2), intent(out) :: map_borders
          real(rkind), dimension(2)  , intent(out) :: map_spacing
          
          integer :: ncid          
          
          ! open the netcdf file for reading
          call nf90_handle_err(
     $         NF90_OPEN(
     $         trim(filename),
     $         NF90_NOWRITE,
     $         ncid))
          
          ! extract the map properties from the netcdf file
          call nf90_get_map_borders(
     $         ncid,
     $         'x',
     $         map_borders(1,:),
     $         map_spacing(1))

          call nf90_get_map_borders(
     $         ncid,
     $         'y',
     $         map_borders(2,:),
     $         map_spacing(2))

          !close netcdf file
          call nf90_handle_err(
     $         NF90_CLOSE(ncid))

        end subroutine nf90_get_coordinate_borders


        subroutine nf90_get_map_borders(
     $     ncid,
     $     map_name,
     $     map_borders,
     $     map_spacing)

          implicit none

          integer                  , intent(in)  :: ncid
          character*(*)            , intent(in)  :: map_name
          real(rkind), dimension(2), intent(out) :: map_borders
          real(rkind)              , intent(out) :: map_spacing

          integer :: varid
          integer :: size_map
          
          real(rkind), dimension(:), allocatable :: map_domain 

          call nf90_handle_err(
     $         NF90_INQ_VARID(
     $         ncid,
     $         map_name,
     $         varid))
          
          call nf90_handle_err(
     $         NF90_INQUIRE_DIMENSION(
     $         ncid,
     $         varid,
     $         len=size_map))
          
          allocate(map_domain(size_map))
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid,
     $         varid,
     $         map_domain))
          
          !extract the borders
          map_borders(1) = map_domain(1)
          map_borders(2) = map_domain(size_map)
          map_spacing    = map_domain(2)-map_domain(1)
          
          deallocate(map_domain)

        end subroutine nf90_get_map_borders


        subroutine nf90_handle_err(nf90_err)

          implicit none

          integer, intent(in) ::  nf90_err
        
          if (nf90_err.ne.NF90_NOERR) then
             print *, 'Error: ', NF90_STRERROR(nf90_err)
             stop 2
          end if

        end subroutine nf90_handle_err

      end module nf90_operators_mt_module
