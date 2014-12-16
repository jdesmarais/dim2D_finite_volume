      program test_write_netcdf_detectors

        use netcdf

        use nf90_operators_module, only :
     $       nf90_handle_err,
     $       nf90_open_file,
     $       nf90_close_file

        use parameters_bf_layer, only :
     $       N,S,E,W

        use parameters_kind, only :
     $       rkind

        
        real(rkind), dimension(2,4) :: N_detectors
        real(rkind), dimension(2,4) :: S_detectors
        real(rkind), dimension(2,4) :: E_detectors
        real(rkind), dimension(2,4) :: W_detectors

        character*(*), parameter       :: filename = 'detectors0.nc'
        integer                        :: ncid
        integer                        :: retval

        integer                        :: k

        character(8) , dimension(4)    :: dct_size_name
        integer      , dimension(4)    :: dct_size
        character(8) , dimension(4)    :: dct_name
        integer      , dimension(4)    :: dct_id

        character(9)                   :: xy_size_name
        integer                        :: xy_size
        integer                        :: xy_dimid

        integer      , dimension(4)    :: dct_dimid
        integer      , dimension(2)    :: dimids
        integer                        :: NF_MYREAL

        integer      , dimension(2)    :: start_op
        integer      , dimension(2)    :: count_op


        !detector definitions
        !------------------------------------------------------------
        dct_size_name = [
     $       'size_N_dct',
     $       'size_S_dct',
     $       'size_E_dct',
     $       'size_W_dct']

        dct_size = [4,4,4,4]

        dct_name = [
     $       'N_detectors',
     $       'S_detectors',
     $       'E_detectors',
     $       'W_detectors']

        xy_size_name = 'xy_coords'
        xy_size      = 2



        N_detectors = reshape((/
     $       -1.0, -1.0,
     $        1.0, -1.0,
     $        1.0,  1.0,
     $       -1.0,  1.0/),
     $       (/2,4/))

        S_detectors = reshape((/
     $       -2.0, -2.0,
     $        2.0, -2.0,
     $        2.0,  2.0,
     $       -2.0,  2.0/),
     $       (/2,4/))

        E_detectors = reshape((/
     $       -3.0, -3.0,
     $        3.0, -3.0,
     $        3.0,  3.0,
     $       -3.0,  3.0/),
     $       (/2,4/))

        W_detectors = reshape((/
     $       -4.0, -4.0,
     $        4.0, -4.0,
     $        4.0,  4.0,
     $       -4.0,  4.0/),
     $       (/2,4/))


        !variable type defintion
        !------------------------------------------------------------
        select case(RKIND)
           case(4)
              NF_MYREAL=NF90_FLOAT
           case(8)
              NF_MYREAL=NF90_DOUBLE
           case default
              print '(''nf90_operators_wr_class :'')'
              print '(''nf90_def_var_model'')'
              stop 'NF_MYREAL'
        end select
           

        !open file for writing
        !------------------------------------------------------------
        call nf90_open_file(filename,ncid)


        !define the variables
        !------------------------------------------------------------
        !dimension of the detector tables
        !............................................................
        do k=1,4

           retval = nf90_def_dim(
     $          ncid,
     $          dct_size_name(k),
     $          dct_size(k),
     $          dct_dimid(k))
           !dec$ forceinline recursive
           call nf90_handle_err(retval)

        end do

        retval = nf90_def_dim(
     $          ncid,
     $          xy_size_name,
     $          xy_size,
     $          xy_dimid)
        !dec$ forceinline recursive
        call nf90_handle_err(retval)

        
        !define the variables
        !............................................................
        dimids(1) = xy_dimid

        do k=1,4

           dimids(2) = dct_dimid(k)

           retval = NF90_DEF_VAR(
     $          ncid,
     $          dct_name(k),
     $          NF_MYREAL,
     $          dimids,
     $          dct_id(k))
           !DEC$ FORCEINLINE RECURSIVE
           call nf90_handle_err(retval)

        end do

        !end the defintion of variables
        !............................................................
        retval = NF90_ENDDEF(ncid)
        !DEC$ FORCEINLINE RECURSIVE
        call nf90_handle_err(retval)
        

        !put the variables
        !------------------------------------------------------------
        start_op = [1,1]

        !N detectors
        count_op = [2,dct_size(N)]
        
        retval = NF90_PUT_VAR(
     $       ncid,
     $       dct_id(N),
     $       N_detectors(:,:),
     $       START=start_op,
     $       COUNT=count_op)
        !DEC$ FORCEINLINE RECURSIVE
        call nf90_handle_err(retval)


        !S detectors
        count_op = [2,dct_size(S)]
        
        retval = NF90_PUT_VAR(
     $       ncid,
     $       dct_id(S),
     $       S_detectors(:,:),
     $       START=start_op,
     $       COUNT=count_op)
        !DEC$ FORCEINLINE RECURSIVE
        call nf90_handle_err(retval)


        !E detectors
        count_op = [2,dct_size(E)]
        
        retval = NF90_PUT_VAR(
     $       ncid,
     $       dct_id(E),
     $       E_detectors(:,:),
     $       START=start_op,
     $       COUNT=count_op)
        !DEC$ FORCEINLINE RECURSIVE
        call nf90_handle_err(retval)


        !W detectors
        count_op = [2,dct_size(W)]
        
        retval = NF90_PUT_VAR(
     $       ncid,
     $       dct_id(W),
     $       W_detectors(:,:),
     $       START=start_op,
     $       COUNT=count_op)
        !DEC$ FORCEINLINE RECURSIVE
        call nf90_handle_err(retval)


        !close file
        !------------------------------------------------------------
        call nf90_close_file(ncid)

      end program test_write_netcdf_detectors
