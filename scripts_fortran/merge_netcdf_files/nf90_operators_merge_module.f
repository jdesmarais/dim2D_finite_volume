      module nf90_operators_merge_module

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
     $       nf90_merge_files,
     $       nf90_merge_files_at_timestep,
     $       nf90_put_var_merge,
     $       nf90_def_var_merge,
     $       nf90_def_header_merge,
     $       nf90_handle_err,
     $       get_filename

        contains


        ! merge all files corresponding to a simulation
        subroutine nf90_merge_files(
     $       nb_tiles,
     $       ne,bc_size,
     $       nb_timesteps)

          implicit none

          integer, dimension(2), intent(in) :: nb_tiles
          integer              , intent(in) :: ne
          integer              , intent(in) :: bc_size
          integer              , intent(in) :: nb_timesteps          
          

          character(len=16) :: filename_rank0
          character(len=16) :: filename_merge
          integer           :: ncid_merge
          integer           :: ierror

          integer, dimension(3) :: coordinates_id
          integer, dimension(ne):: var_id
                    
          integer :: timestep

          real(rkind), dimension(1)                  :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes

          integer :: size_x
          integer :: size_y


          ! this is just to get the size_x and size_y
          ! such that the allocation of x_map, y_map,
          ! and nodes happens only once in this
          ! subroutine instead of at each timestep
          call get_filename(
     $         filename_rank0,
     $         0,
     $         rank=0)

          call get_filename(
     $         filename_merge,
     $         0)

          call nf90_def_header_merge(
     $         filename_rank0,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)

          call nf90_def_var_merge(
     $         filename_rank0,
     $         ncid_merge,
     $         bc_size,nb_tiles,
     $         coordinates_id,
     $         var_id,
     $         size_x=size_x,
     $         size_y=size_y,
     $         ierror=ierror)

          call nf90_handle_err(NF90_CLOSE(ncid_merge))

          allocate(x_map(size_x))
          allocate(y_map(size_y))
          allocate(nodes(size_x,size_y,ne))


          ! merge the files for all timesteps
          do timestep=0,nb_timesteps-1

             call nf90_merge_files_at_timestep(
     $            nb_tiles,
     $            ne,
     $            bc_size,
     $            timestep,
     $            time,
     $            x_map,
     $            y_map,
     $            nodes)

             print '(''timestep: '', I3)', timestep

          end do

          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)

        end subroutine nf90_merge_files


        !merge the files corresponding to the same timestep
        subroutine nf90_merge_files_at_timestep(
     $       nb_tiles,
     $       ne,
     $       bc_size,
     $       timestep,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes)

          implicit none

          integer    , dimension(2)    , intent(in)    :: nb_tiles
          integer                      , intent(in)    :: ne
          integer                      , intent(in)    :: bc_size
          integer                      , intent(in)    :: timestep
          real(rkind), dimension(1)    , intent(inout) :: time
          real(rkind), dimension(:)    , intent(inout) :: x_map
          real(rkind), dimension(:)    , intent(inout) :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes


          integer :: ncid_merge
          integer :: ierror

          integer, dimension(3) :: coordinates_id
          integer, dimension(ne):: var_id
                    
          character(len=16) :: filename_rank
          character(len=16) :: filename_merge

          integer :: k

          
          ! define the name of the files
          call get_filename(
     $         filename_rank,
     $         timestep,
     $         rank=0)

          call get_filename(
     $         filename_merge,
     $         timestep)


          ! define the header of the merged file
          call nf90_def_header_merge(
     $         filename_rank,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)

          ! define the variables stored in the merged file
          call nf90_def_var_merge(
     $         filename_rank,
     $         ncid_merge,
     $         bc_size,nb_tiles,
     $         coordinates_id,
     $         var_id,
     $         ierror=ierror)


          ! combine the data from the different processor files 
          ! into one merged file
          do k=0, nb_tiles(1)*nb_tiles(2)-1

             call get_filename(
     $            filename_rank,
     $            timestep,
     $            rank=k)

             call nf90_put_var_merge(
     $            trim(filename_rank),
     $            k,
     $            ncid_merge,
     $            bc_size,
     $            nb_tiles,
     $            coordinates_id,
     $            var_id,
     $            time,
     $            x_map,
     $            y_map,
     $            nodes,
     $            ierror=ierror)

          end do
          
          call nf90_handle_err(NF90_CLOSE(ncid_merge))

        end subroutine nf90_merge_files_at_timestep


        !extract the variable from the file and store them in the
        !merge file
        ! - from the rank of the processor which created the output
        !   file, determine where the extracted variables should be
        !   saved
        subroutine nf90_put_var_merge(
     $       filename_rank,
     $       rank,
     $       ncid_merge,
     $       bc_size,
     $       nb_tiles,
     $       coordinates_id,
     $       var_id,
     $       time,
     $       x_map,
     $       y_map,
     $       nodes,
     $       ierror)

          implicit none

          character*(*)                  , intent(in)    :: filename_rank
          integer                        , intent(in)    :: rank
          integer                        , intent(in)    :: ncid_merge
          integer                        , intent(in)    :: bc_size
          integer      , dimension(2)    , intent(in)    :: nb_tiles
          integer      , dimension(3)    , intent(in)    :: coordinates_id
          integer      , dimension(:)    , intent(in)    :: var_id
          real(rkind)  , dimension(1)    , intent(inout) :: time
          real(rkind)  , dimension(:)    , intent(inout) :: x_map
          real(rkind)  , dimension(:)    , intent(inout) :: y_map
          real(rkind)  , dimension(:,:,:), intent(inout) :: nodes
          integer      , optional        , intent(out)   :: ierror


          integer :: ncid_rank

          integer :: size_x_rank
          integer :: size_y_rank
          integer :: ne

          integer, dimension(2) :: cart_coords
          integer, dimension(2) :: count
          integer, dimension(2) :: start
          integer, dimension(4) :: limit

          integer :: NF_MYREAL

          integer :: k


          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          size_x_rank = size(x_map,1)
          size_y_rank = size(y_map,1)
          ne          = size(nodes,3)


          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)

          
          !2) open the small domain netcdf file
          call nf90_handle_err(
     $         NF90_OPEN(
     $         trim(filename_rank),
     $         NF90_NOWRITE,
     $         ncid_rank))


          !3) extract the cartesian coordinates
          !   from the rank of the processor which
          !   computed it
          !------------------------------------------------------------
          cart_coords(2) = mod(rank,nb_tiles(2))
          cart_coords(1) = (rank-cart_coords(2))/nb_tiles(2)


          !4) deduce the indices where the variables
          !   extracted from the file will be saved
          !------------------------------------------------------------
          count(1) = size_x_rank-2*bc_size
          count(2) = size_y_rank-2*bc_size

          start(1) = 1+bc_size+cart_coords(1)*(size_x_rank-2*bc_size)
          start(2) = 1+bc_size+cart_coords(2)*(size_y_rank-2*bc_size)

          limit(1) = 1+bc_size
          limit(2) = size_x_rank-bc_size
          limit(3) = 1+bc_size
          limit(4) = size_y_rank-bc_size

          if(cart_coords(1).eq.0) then
             limit(1)=limit(1)-bc_size
             start(1)=start(1)-bc_size
             count(1)=count(1)+bc_size
          end if

          if(cart_coords(1).eq.(nb_tiles(1)-1)) then
             limit(2)=limit(2)+bc_size
             count(1)=count(1)+bc_size
          end if

          if(cart_coords(2).eq.0) then
             limit(3)=limit(3)-bc_size
             start(2)=start(2)-bc_size
             count(2)=count(2)+bc_size
          end if

          if(cart_coords(2).eq.(nb_tiles(2)-1)) then
             limit(4)=limit(4)+bc_size
             count(2)=count(2)+bc_size
          end if


          !5) extract the time, coordinates, and governing variables
          !   from the rank file
          !------------------------------------------------------------
          !5.1) get the time
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_rank,
     $         coordinates_id(1),
     $         time))
          
          !5.2) get the x_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_rank,
     $         coordinates_id(2),
     $         x_map))
          
          !5.3) get the y_map
          call nf90_handle_err(
     $         NF90_GET_VAR(
     $         ncid_rank,
     $         coordinates_id(3),
     $         y_map))             
      
          !5.4) get the governing variables
          do k=1,ne
          
             call nf90_handle_err(
     $            NF90_GET_VAR(
     $            ncid_rank,
     $            var_id(k),
     $            nodes(:,:,k),
     $            START=[1,1,1],
     $            COUNT=[1,size_x_rank,size_y_rank]))
          end do


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
     $         NF90_CLOSE(ncid_rank))

          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end subroutine nf90_put_var_merge


        !define the coordinates and the variables saved in the merge
        !file
        ! - extract the name of the variables + units + size from
        !   the rank0 file
        ! - from these sizes and the number of tiles in the x- and y-
        !   directions, deduce the total size of the merged domain
        ! - define the coordinates and the variables in the merge
        !   file
        !------------------------------------------------------------
        subroutine nf90_def_var_merge(
     $       filename_rank0,
     $       ncid_merge,
     $       bc_size,
     $       nb_tiles,
     $       coordinates_id,
     $       var_id,
     $       size_x,
     $       size_y,
     $       ierror)

          implicit none

          character*(*)              , intent(in)    :: filename_rank0
          integer                    , intent(in)    :: ncid_merge
          integer                    , intent(in)    :: bc_size
          integer      , dimension(2), intent(in)    :: nb_tiles
          integer      , dimension(3), intent(out)   :: coordinates_id
          integer      , dimension(:), intent(inout) :: var_id
          integer      , optional    , intent(out)   :: size_x
          integer      , optional    , intent(out)   :: size_y
          integer      , optional    , intent(out)   :: ierror


          integer               :: size_x_rank0
          integer               :: size_y_rank0
          integer, dimension(3) :: sizes

          integer               :: NF_MYREAL
          integer               :: ncid_rank0
          integer, dimension(3) :: dims_id
          integer               :: k,l

          character(len=60)   :: nameVar

          integer             :: nAttsVar
          character(len=60)   :: nameAtt
          character(len=1024) :: charVarAtt

          integer :: ne

          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)

          
          !2) open the small domain netcdf file
          call nf90_handle_err(
     $         NF90_OPEN(
     $         trim(filename_rank0),
     $         NF90_NOWRITE,
     $         ncid_rank0))


          !3) extract the name, units and longname
          !   of the coordinates
          do k=1,3
             
             ! extract the name of the coordinate
             ! and save it in the error_file
             call nf90_handle_err(
     $            NF90_INQUIRE_VARIABLE(
     $            ncid_rank0,
     $            k,
     $            name=nameVar,
     $            nAtts=nAttsVar))
             
             
             ! define the dimension with the name
             ! and compute the size for the coordinate
             select case(nameVar)
             
               case('time')

                  sizes(1) = 1

               case('x')

                  !extract the size of the 'x'-coordinates
                  !in the rank0 file
                  call nf90_handle_err(
     $                 NF90_INQUIRE_DIMENSION(
     $                 ncid_rank0,
     $                 k,
     $                 len=size_x_rank0))

                  sizes(2) = (size_x_rank0 - 2*bc_size)*nb_tiles(1)
     $                       + 2*bc_size

                  if(present(size_x)) then
                     size_x = size_x_rank0
                  end if

               case('y')

                  !extract the size of the 'y'-coordinates
                  !in the rank0 file
                  call nf90_handle_err(
     $                 NF90_INQUIRE_DIMENSION(
     $                 ncid_rank0,
     $                 k,
     $                 len=size_y_rank0))

                  sizes(3) = (size_y_rank0 - 2*bc_size)*nb_tiles(2)
     $                       + 2*bc_size                  

                  if(present(size_y)) then
                     size_y = size_y_rank0
                  end if

             end select             


             ! define the size of the coordinate
             call nf90_handle_err(
     $            NF90_DEF_DIM(
     $            ncid_merge,
     $            trim(nameVar),
     $            sizes(k),
     $            dims_id(k)))


             ! define the corresponding variable
             call nf90_handle_err(
     $            NF90_DEF_VAR(
     $            ncid_merge,
     $            trim(nameVar),
     $            NF_MYREAL,
     $            dims_id(k),
     $            coordinates_id(k)))


             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                call nf90_handle_err(
     $               NF90_INQ_ATTNAME(
     $               ncid_rank0,
     $               k,
     $               l,
     $               nameAtt))

                call nf90_handle_err(
     $               NF90_GET_ATT(
     $               ncid_rank0,
     $               k,
     $               trim(nameAtt),
     $               charVarAtt))
                
                call nf90_handle_err(
     $               NF90_PUT_ATT(
     $               ncid_merge,
     $               coordinates_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt)))
                
             end do

          end do



          !4) extract and define the governing 
          !   variables
          !   [1,3]: coordinates
          !   [4,:]: governing variables
          ne = size(var_id,1)

          do k=1, ne

             ! extract the governing variable
             ! from the netcdf file of the
             ! small domain
             call nf90_handle_err(
     $            NF90_INQUIRE_VARIABLE(
     $            ncid_rank0,
     $            3+k,
     $            name=nameVar,
     $            nAtts=nAttsVar))


             ! define the governing variable to
             ! be saved in the netcdf file 
             ! containing the error
             call nf90_handle_err(
     $            NF90_DEF_VAR(
     $            ncid_merge,
     $            nameVar,
     $            NF_MYREAL,
     $            dims_id,
     $            var_id(k)))

             
             ! extract and assign the units to
             ! the governing variables
             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                call nf90_handle_err(
     $               NF90_INQ_ATTNAME(
     $               ncid_rank0,
     $               3+k,
     $               l,
     $               nameAtt))

                call nf90_handle_err(
     $               NF90_GET_ATT(
     $               ncid_rank0,
     $               3+k,
     $               trim(nameAtt),
     $               charVarAtt))

                call nf90_handle_err(
     $               NF90_PUT_ATT(
     $               ncid_merge,
     $               var_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt)))
                
             end do

          end do

          call nf90_handle_err(NF90_CLOSE(ncid_rank0))

          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end subroutine nf90_def_var_merge


        !define the header of the merge file
        ! - create the file for merging the netcdf files
        ! - read the rank 0 file and copy the header in
        !   the merge file
        ! - define the coordinates and the variables ID
        !   by computing the total size of the domain
        !------------------------------------------------------------
        subroutine nf90_def_header_merge(
     $       filename_rank0,
     $       filename_merge,
     $       ncid_merge,
     $       ierror)

          implicit none


          character*(*)          , intent(in)  :: filename_rank0
          character*(*)          , intent(in)  :: filename_merge
          integer                , intent(out) :: ncid_merge
          integer      , optional, intent(out) :: ierror

          integer               :: retval
          integer               :: ncid_rank0

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


          !5) close netcdf files
          retval = NF90_CLOSE(ncid_rank0)
          call nf90_handle_err(retval)
          
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


      end module nf90_operators_merge_module
