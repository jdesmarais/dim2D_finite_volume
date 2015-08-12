      program extract_energyTr

        use cmd_operators_extract_class, only :
     $       cmd_operators_extract

        use energyTr_computation_module, only :
     $       compute_energyTr

        use nf90_operators_error_module, only :
     $       nf90_get_maps,
     $       nf90_get_gov_var

        use parameters_cst, only :
     $       SUCCESS

        use parameters_kind, only :
     $       rkind


        implicit none

        type(cmd_operators_extract)                :: cmd_op

        character(len=100)                         :: filename

        real(rkind), dimension(:)    , allocatable :: x_map
        real(rkind), dimension(:)    , allocatable :: y_map
        real(rkind), dimension(:,:,:), allocatable :: nodes

        real(rkind) :: energyTr

        integer :: ierror
        

        ! analyse the command line arguments
        call cmd_op%analyse_cmd_line_arg()
        filename = cmd_op%get_filename()


        ! extract the data from the netcdf file
        ! (mass, momentum_x, momentum_y, energy)
        call nf90_get_maps(
     $       filename,
     $       x_map,
     $       y_map)

        call nf90_get_gov_var(
     $       filename,
     $       ['mass','momentum_x','momentum_y','energy'],
     $       nodes,
     $       ierror=ierror)

        if(ierror.ne.SUCCESS) then
           print '(''data extraction failed'')'
        end if


        ! extract the energy balance through the domain borders
        energyTr = compute_energyTr(
     $       x_map,
     $       y_map,
     $       nodes,
     $       borders=[.true.,.false.,.true.,.true.])


        ! output the energy balance
        print *, 'energyTr: ', energyTr

      end program extract_energyTr
