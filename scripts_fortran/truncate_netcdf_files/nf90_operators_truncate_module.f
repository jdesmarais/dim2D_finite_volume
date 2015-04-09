      module nf90_operators_truncate_module

        use check_data_module, only : 
     $     is_real_validated

        use netcdf

        use parameters_cst, only :
     $       NOT_SUCCESS,
     $       SUCCESS

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: 
     $       nf90_truncate_files,
     $       nf90_truncate_file_at_timestep,
     $       nf90_put_var_truncate,
     $       nf90_def_var_truncate,
     $       nf90_get_extraction_indices,
     $       nf90_def_header_truncate,
     $       nf90_handle_err,
     $       get_filename

        contains


        ! merge all files corresponding to a simulation
        ! usign processors working in parallel
        subroutine nf90_truncate_files(
     $       ne,
     $       nb_timesteps,
     $       borders_tr,
     $       timestep_start,
     $       timestep_increment)

          implicit none

          integer                    , intent(in) :: ne
          integer                    , intent(in) :: nb_timesteps  
          real(rkind), dimension(2,2), intent(in) :: borders_tr
          integer    , optional      , intent(in) :: timestep_start
          integer    , optional      , intent(in) :: timestep_increment
          

          character(len=16) :: filename_lg
          integer           :: ncid_lg
          integer           :: ierror

          integer(ikind), dimension(2,2) :: borders_indices

          real(rkind), dimension(1)                  :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes

          integer :: size_x
          integer :: size_y

          integer :: timestep
          integer :: timestep_start_op
          integer :: timestep_increment_op


          if(present(timestep_start)) then
             timestep_start_op = timestep_start
          else
             timestep_start_op = 0
          end if


          if(present(timestep_increment)) then
             timestep_increment_op = timestep_increment
          else
             timestep_increment_op = 1
          end if


          ! this is just to get the size_x and size_y
          ! such that the allocation of x_map, y_map,
          ! and nodes happens only once in this
          ! subroutine instead of at each timestep
          call get_filename(
     $         filename_lg,
     $         timestep_start_op)

          call nf90_handle_err(
     $         NF90_OPEN(
     $            trim(filename_lg),
     $            NF90_NOWRITE,
     $            ncid_lg))

          call nf90_get_extraction_indices(
     $         ncid_lg,
     $         borders_tr,
     $         borders_indices,
     $         ierror)

          size_x = borders_indices(1,2)-borders_indices(1,1)+1
          size_y = borders_indices(2,2)-borders_indices(2,1)+1
          
          call nf90_handle_err(NF90_CLOSE(ncid_lg))

          allocate(x_map(size_x))
          allocate(y_map(size_y))
          allocate(nodes(size_x,size_y,ne))


          ! truncate the files for all timesteps
          do timestep=timestep_start_op,nb_timesteps-1,timestep_increment_op

             call nf90_truncate_file_at_timestep(
     $            timestep,
     $            ne,
     $            borders_indices,
     $            time,
     $            x_map,
     $            y_map,
     $            nodes)

          end do

          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)

        end subroutine nf90_truncate_files


        !extract the domain for a timestep
        subroutine nf90_truncate_file_at_timestep(
     $       timestep,
     $       ne,
     $       borders_indices,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes)
          
          implicit none

          integer                         , intent(in)  :: timestep
          integer                         , intent(in)  :: ne
          integer(ikind), dimension(2,2)  , intent(in)  :: borders_indices
          real(rkind)   , dimension(1)    , intent(out) :: time
          real(rkind)   , dimension(:)    , intent(out) :: x_map
          real(rkind)   , dimension(:)    , intent(out) :: y_map
          real(rkind)   , dimension(:,:,:), intent(out) :: nodes
          

          character(len=16) :: filename_lg
          character(len=27) :: filename_tr
          
          integer :: ncid_lg
          integer :: ncid_tr

          integer :: ierror

          integer, dimension(3)              :: coordinates_id
          integer, dimension(:), allocatable :: var_id


          allocate(var_id(ne))

          call get_filename(
     $         filename_lg,
     $         timestep)

          filename_tr(1:11) = './tr_files/'
          filename_tr(12:27)= filename_lg

          ! define the header of the truncated file
          call nf90_def_header_truncate(
     $         filename_lg,
     $         filename_tr,
     $         ncid_lg,
     $         ncid_tr,
     $         ierror)

          ! define the variables in the truncated file
          call nf90_def_var_truncate(
     $         ncid_lg,
     $         ncid_tr,
     $         borders_indices,
     $         coordinates_id,
     $         var_id)

          ! put the variables in the truncated file
          call nf90_put_var_truncate(
     $         ncid_lg,
     $         ncid_tr,
     $         borders_indices,
     $         coordinates_id,
     $         var_id,
     $         time,
     $         x_map,
     $         y_map,
     $         nodes,
     $         ierror)

          ! close the files
          call nf90_handle_err(NF90_CLOSE(ncid_lg))
          call nf90_handle_err(NF90_CLOSE(ncid_tr))

        end subroutine nf90_truncate_file_at_timestep


        !extract the variable from the large domain file and
        !store them in the file of the truncated domain
        subroutine nf90_put_var_truncate(
     $       ncid_lg,
     $       ncid_tr,
     $       borders_indices,
     $       coordinates_id,
     $       var_id,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes,
     $       ierror)

          implicit none

          integer                        , intent(in)    :: ncid_lg
          integer                        , intent(in)    :: ncid_tr
          integer      , dimension(2,2)  , intent(in)    :: borders_indices
          integer      , dimension(3)    , intent(in)    :: coordinates_id
          integer      , dimension(:)    , intent(in)    :: var_id
          real(rkind)  , dimension(1)    , intent(inout) :: time
          real(rkind)  , dimension(:)    , intent(inout) :: x_map
          real(rkind)  , dimension(:)    , intent(inout) :: y_map
          real(rkind)  , dimension(:,:,:), intent(inout) :: nodes
          integer      , optional        , intent(out)   :: ierror


          integer :: ne

          integer, dimension(2) :: count
          integer, dimension(2) :: start

          integer :: NF_MYREAL

          integer :: k


          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          ne = size(nodes,3)


          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)


          !4) deduce the indices where the variables
          !   extracted from the file will be saved
          !------------------------------------------------------------
          count(1) = borders_indices(1,2)-borders_indices(1,1)+1
          count(2) = borders_indices(2,2)-borders_indices(2,1)+1

          start(1) = borders_indices(1,1)
          start(2) = borders_indices(2,1)


          !5) extract the time, coordinates, and governing variables
          !   from the rank file
          !------------------------------------------------------------
          !5.1) get the time
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_lg,
     $         coordinates_id(1),
     $         time))
          
          !5.2) get the x_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_lg,
     $         coordinates_id(2),
     $         x_map,
     $         START=[start(1)],
     $         COUNT=[count(1)]))

          !5.3) get the y_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_lg,
     $         coordinates_id(3),
     $         y_map,
     $         START=[start(2)],
     $         COUNT=[count(2)]))

      
          !5.4) get the governing variables
          do k=1,ne
          
             call nf90_handle_err(
     $            NF90_GET_VAR(
     $            ncid_lg,
     $            var_id(k),
     $            nodes(:,:,k),
     $            START=[1,start(1),start(2)],
     $            COUNT=[1,count(1),count(2)]))

          end do


          !6) save the time, coordinates, and governing variables
          !   in the merge file
          !------------------------------------------------------------
          !6.1) write the time
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_tr,
     $         coordinates_id(1),
     $         time))

          !6.2) write the x-coordinates
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_tr,
     $         coordinates_id(2),
     $         x_map))

          !6.3) write the y-coordinates
          call nf90_handle_err(
     $         NF90_PUT_VAR(
     $         ncid_tr,
     $         coordinates_id(3),
     $         y_map))


          !6.4) write the variables of the governing equations
          do k=1, ne

             call nf90_handle_err(
     $            NF90_PUT_VAR(
     $            ncid_tr,
     $            var_id(k),
     $            nodes(:,:,k),
     $            START=[1,1,1],
     $            COUNT=[1,count(1),count(2)]))
          end do

          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end subroutine nf90_put_var_truncate


        ! define the variables saved in the file
        ! for the truncated domain knowing the
        ! truncation indices
        subroutine nf90_def_var_truncate(
     $       ncid_lg_domain,
     $       ncid_truncate,
     $       borders_indices,
     $       coordinates_id,
     $       var_id)

          implicit none

          integer                       , intent(in)  :: ncid_lg_domain
          integer                       , intent(in)  :: ncid_truncate
          integer(ikind), dimension(2,2), intent(in)  :: borders_indices
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
          sizes = [
     $         1,
     $         borders_indices(1,2)-borders_indices(1,1)+1,
     $         borders_indices(2,2)-borders_indices(2,1)+1]



          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)

          
          !3) extract the name, units and longname
          !   of the coordinates
          do k=1,3
             
             ! extract the name of the coordinate
             ! and save it in the error_file
             retval = NF90_INQUIRE_VARIABLE(
     $            ncid_lg_domain,
     $            k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)
             

             ! define the dimension with the name
             retval = NF90_DEF_DIM(
     $            ncid_truncate,
     $            trim(nameVar),
     $            sizes(k),
     $            dims_id(k))
             call nf90_handle_err(retval)


             ! define the corresponding variable
             retval = NF90_DEF_VAR(
     $            ncid_truncate,
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
     $               ncid_lg_domain,
     $               k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_lg_domain,
     $               k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                retval = NF90_PUT_ATT(
     $               ncid_truncate,
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
     $            ncid_lg_domain,
     $            3+k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)

             ! define the governing variable to
             ! be saved in the netcdf file 
             ! containing the truncation
             retval = NF90_DEF_VAR(
     $            ncid_truncate,
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
     $               ncid_lg_domain,
     $               3+k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_lg_domain,
     $               3+k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                retval = NF90_PUT_ATT(
     $               ncid_truncate,
     $               var_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt))
                call nf90_handle_err(retval)
                
             end do

          end do

        end subroutine nf90_def_var_truncate


        ! get the indices matching the truncated
        ! domain in the large domain file
        subroutine nf90_get_extraction_indices(
     $       ncid,
     $       borders_tr,
     $       borders_indices,
     $       ierror)

          implicit none

          integer                       , intent(in)  :: ncid
          real(rkind)   , dimension(2,2), intent(in)  :: borders_tr
          integer(ikind), dimension(2,2), intent(out) :: borders_indices
          integer                       , intent(out) :: ierror

          
          integer(ikind), dimension(2) :: x_borders_indices
          integer(ikind), dimension(2) :: y_borders_indices


          ! get the extraction indices for the x_map
          call nf90_get_map_extraction_indices(
     $         ncid,
     $         'x',
     $         [borders_tr(1,1),borders_tr(1,2)],
     $         x_borders_indices,
     $         ierror)
          if(ierror.ne.SUCCESS) then
             print '(''nf90_operators_truncate_module'')'
             print '(''nf90_get_extraction_indices'')'
             print '(''error when extracting x-borders'')'
          end if


          ! get the extraction indices for the y_map
          call nf90_get_map_extraction_indices(
     $         ncid,
     $         'y',
     $         [borders_tr(2,1),borders_tr(2,2)],
     $         y_borders_indices,
     $         ierror)
          if(ierror.ne.SUCCESS) then
             print '(''nf90_operators_truncate_module'')'
             print '(''nf90_get_extraction_indices'')'
             print '(''error when extracting y-borders'')'
          end if


          ! combine the extraction indices
          borders_indices(1,1) = x_borders_indices(1)
          borders_indices(1,2) = x_borders_indices(2)
          borders_indices(2,1) = y_borders_indices(1)
          borders_indices(2,2) = y_borders_indices(2)

        end subroutine nf90_get_extraction_indices


        !get the indices matching the truncated domain in the
        !large domain file
        subroutine nf90_get_map_extraction_indices(
     $     ncid,
     $     map_name,
     $     borders_tr,
     $     borders_indices,
     $     ierror)

          implicit none

          integer                     , intent(in)  :: ncid
          character*(*)               , intent(in)  :: map_name
          real(rkind)   , dimension(2), intent(in)  :: borders_tr
          integer(ikind), dimension(2), intent(out) :: borders_indices
          integer                     , intent(out) :: ierror

          integer :: varid
          integer :: size_map_lg_domain


          real(rkind), dimension(:), allocatable :: map_lg_domain


          ! extract the map from the netcdf file
          call nf90_handle_err(
     $         NF90_INQ_VARID(
     $            ncid,
     $            map_name,
     $            varid))

          call nf90_handle_err(
     $         NF90_INQUIRE_DIMENSION(
     $            ncid,
     $            varid,
     $            len=size_map_lg_domain))

          allocate(map_lg_domain(size_map_lg_domain))
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $            ncid,
     $            varid,
     $            map_lg_domain))


          ! compute the indices matching the truncated domain
          borders_indices(1) = nf90_get_map_extraction_index(
     $         borders_tr(1),
     $         map_lg_domain,
     $         side=.false.,
     $         ierror=ierror)
          if(ierror.ne.SUCCESS) then
             print '(''nf90_operators_truncate_module'')'
             print '(''***error when matching min in large domain***'')'
             stop ''
          end if

          borders_indices(2) = nf90_get_map_extraction_index(
     $         borders_tr(2),
     $         map_lg_domain,
     $         side=.true.,
     $         ierror=ierror)
          if(ierror.ne.SUCCESS) then
             print '(''nf90_operators_truncate_module'')'
             print '(''***error when matching max in large domain***'')'
             stop ''
          end if


          ! deallocate the temporary map
          deallocate(map_lg_domain)

        end subroutine nf90_get_map_extraction_indices


        ! get the indices matching the truncated domain in the
        ! large domain coordinate-map
        function nf90_get_map_extraction_index(
     $     coord,
     $     coord_map,
     $     side,
     $     detailled,
     $     ierror)
     $     result(index)

          implicit none

          real(rkind)              , intent(in)  :: coord
          real(rkind), dimension(:), intent(in)  :: coord_map
          logical    , optional    , intent(in)  :: side
          logical    , optional    , intent(in)  :: detailled
          integer    , optional    , intent(out) :: ierror
          integer                                :: index

          logical :: side_op
          logical :: detailled_op
          integer :: ierror_op

          integer :: i

          ierror_op = NOT_SUCCESS


          if(present(side)) then
             side_op = side
          else
             side_op = .true.
          end if

          
          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if


          !start from the right side
          if(side_op) then

             do i=size(coord_map,1), 1, -1

                if(is_real_validated(
     $               coord,
     $               coord_map(i),
     $               detailled=detailled_op)
     $            ) then

                   index = i
                   ierror_op = SUCCESS
                   exit
                end if

             end do
             

          !start from the left side
          else

             do i=1, size(coord_map,1), +1

                if(is_real_validated(coord,coord_map(i),detailled=detailled_op)) then
                   index = i
                   ierror_op = SUCCESS
                   exit
                end if

             end do

          end if


          !set ierror
          if(present(ierror)) then
             ierror = ierror_op
          end if

        end function nf90_get_map_extraction_index


        !define the header of the merge file
        ! - create the file for merging the netcdf files
        ! - read the rank 0 file and copy the header in
        !   the merge file
        ! - define the coordinates and the variables ID
        !   by computing the total size of the domain
        !------------------------------------------------------------
        subroutine nf90_def_header_truncate(
     $       filename_rank0,
     $       filename_merge,
     $       ncid_rank0,
     $       ncid_merge,
     $       ierror)

          implicit none


          character*(*)          , intent(in)  :: filename_rank0
          character*(*)          , intent(in)  :: filename_merge
          integer                , intent(out) :: ncid_rank0
          integer                , intent(out) :: ncid_merge
          integer      , optional, intent(out) :: ierror

          integer               :: retval

          integer               :: size_header_rank0

          integer               :: k

          character(len=1028)   :: header_attname_rank0
          integer               :: header_att_type
          integer               :: header_att_len

          character(len=1028)   :: att_char_var
          integer               :: att_int_var
          real                  :: att_real_var
          real*8                :: att_double_var

        
          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          !1) open netcdf file for saving the header
          !   as a combination of the small and large
          !   domains headers
          retval = NF90_CREATE(
     $         trim(filename_merge),
     $         NF90_NETCDF4,
     $         ncid_merge)
          call nf90_handle_err(retval)


          !2) open second netcdf file for reading
          !   small domain header
          retval = NF90_OPEN(
     $         trim(filename_rank0),
     $         NF90_NOWRITE,
     $         ncid_rank0)
          call nf90_handle_err(retval)

          
          !3) inquire the number of global
          !   variables (header variables)
          retval = NF90_INQUIRE(
     $         ncid_rank0,
     $         nAttributes=size_header_rank0)
          call nf90_handle_err(retval)


          !4) inquire the name of the global
          !   variables for the small and large
          !   domains
          do k=1, size_header_rank0

             ! inquire the name of the global
             ! variables in the header
             retval = NF90_INQ_ATTNAME(
     $            ncid_rank0,
     $            NF90_GLOBAL,
     $            k,
     $            header_attname_rank0)
             call nf90_handle_err(retval)


             ! inquire the type and length of
             ! the global variable
             retval = NF90_INQUIRE_ATTRIBUTE(
     $            ncid_rank0,
     $            NF90_GLOBAL,
     $            trim(header_attname_rank0),
     $            xtype=header_att_type,
     $            len=header_att_len)
             call nf90_handle_err(retval)


             ! inquire the value of the global
             ! variable
             call nf90_get_att_value(
     $            ncid_rank0,
     $            NF90_GLOBAL,
     $            header_att_type,
     $            trim(header_attname_rank0),
     $            att_char_var,
     $            att_int_var,
     $            att_real_var,
     $            att_double_var)

             
             ! put the value in the merge file
             call nf90_put_att_value(
     $            ncid_merge,
     $            NF90_GLOBAL,
     $            header_att_type,
     $            trim(header_attname_rank0),
     $            att_char_var,
     $            att_int_var,
     $            att_real_var,
     $            att_double_var)

          end do


          if(present(ierror)) then
             ierror = SUCCESS
          end if

      end subroutine nf90_def_header_truncate


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


        ! print the netcdf error received
        subroutine nf90_handle_err(nf90_err)

          implicit none

          integer, intent(in) ::  nf90_err
        
          if (nf90_err.ne.NF90_NOERR) then
             print *, 'Error: ', NF90_STRERROR(nf90_err)
             stop 2
          end if

        end subroutine nf90_handle_err


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


        !get the filename corresponding to a timestep
        !and a processor rank
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

      end module nf90_operators_truncate_module
