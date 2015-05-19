      module nf90_operators_mt_module

        use netcdf

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: nf90_get_map_borders
        

        contains


        function nf90_get_map_borders(
     $       filename, map_name)
     $       result(map_borders)

          implicit none
          
          character*(*), intent(in) :: filename
          character*(*), intent(in) :: map_name
          real(rkind), dimension(2) :: map_borders
          
          
          integer :: ncid
          integer :: varid
          integer :: size_map
          
          real(rkind), dimension(:), allocatable :: map_domain          
          
          
          ! open the netcdf file for reading
          call nf90_handle_err(
     $         NF90_OPEN(
     $         trim(filename),
     $         NF90_NOWRITE,
     $         ncid))
          
          ! extract the map from the netcdf file
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
          
          deallocate(map_domain)

        end function nf90_get_map_borders

      end module nf90_operators_mt_modulex
