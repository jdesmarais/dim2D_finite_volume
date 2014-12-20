      program test_dct_verify_x_symmetry

        use bf_restart_module, only :
     $       read_detectors_from_file

        use cmd_operators_class, only :
     $       cmd_operators

        use parameters_bf_layer, only :
     $       N,S,E,W

        use parameters_kind, only :
     $       ikind,
     $       rkind
        

        type(cmd_operators)                      :: cmd_operators_used
        integer    , dimension(4)                :: nb_bf_layers
        integer                                  :: timestep
        real(rkind), dimension(:,:), allocatable :: N_dct_rcoords
        real(rkind), dimension(:,:), allocatable :: S_dct_rcoords
        real(rkind), dimension(:,:), allocatable :: E_dct_rcoords
        real(rkind), dimension(:,:), allocatable :: W_dct_rcoords
        logical, parameter                       :: detailled = .true.

        integer :: k

        !get the timestep to be analyzed
        call cmd_operators_used%analyse_cmd_line_arg()

        nb_bf_layers = cmd_operators_used%get_nb_bf_layers()

        timestep = nb_bf_layers(1)
        do k=2,4
           if(nb_bf_layers(k).gt.timestep) then
              timestep = nb_bf_layers(k)
           end if
        end do        

        !get detectors
        call get_detectors(
     $       timestep,
     $       N_dct_rcoords,
     $       S_dct_rcoords,
     $       E_dct_rcoords,
     $       W_dct_rcoords)

        
        !analyze the symmetry for the detectors
        if(nb_bf_layers(N).ne.0) then
           call analyse_x_symmetry(N_dct_rcoords,detailled)
        end if

        if(nb_bf_layers(S).ne.0) then
           call analyse_x_symmetry(S_dct_rcoords,detailled)
        end if

        if(nb_bf_layers(E).ne.0) then
           call analyse_x_symmetry(E_dct_rcoords,detailled)
        end if

        if(nb_bf_layers(W).ne.0) then
           call analyse_x_symmetry(W_dct_rcoords,detailled)
        end if


        contains

        !analyse the symmetry in the coords
        subroutine analyse_x_symmetry(dct_rcoords,detailled)

          implicit none

          real(rkind), dimension(:,:), intent(in) :: dct_rcoords
          logical                    , intent(in) :: detailled

          integer :: k
          integer :: k_sym
          integer :: nb_dct
          logical :: symmetric_dct
          logical :: symmetric_loc


          symmetric_dct = .true.
          nb_dct        = size(dct_rcoords,2)

          do k=1, int(nb_dct/2)

             k_sym = nb_dct-k+1

             symmetric_loc = is_test_validated(
     $             dct_rcoords(1,k),
     $            -dct_rcoords(1,k_sym),
     $            .false.)
             symmetric_dct = symmetric_dct.and.symmetric_loc

             if(detailled.and.(.not.symmetric_loc)) then

                print '(''['',I4,'']: '',F8.4,''->'',F8.4)',
     $               k,
     $               dct_rcoords(1,k),
     $               dct_rcoords(1,k_sym)

             end if

          end do

          print '(''symmetric detectors?: '',L1)', symmetric_dct
          print '()'

        end subroutine analyse_x_symmetry


        !get the detectors
        subroutine get_detectors(
     $       timestep,
     $       N_dct_rcoords,
     $       S_dct_rcoords,
     $       E_dct_rcoords,
     $       W_dct_rcoords)

          implicit none

          integer                                 , intent(in)  :: timestep
          real(rkind), dimension(:,:), allocatable, intent(out) :: N_dct_rcoords
          real(rkind), dimension(:,:), allocatable, intent(out) :: S_dct_rcoords
          real(rkind), dimension(:,:), allocatable, intent(out) :: E_dct_rcoords
          real(rkind), dimension(:,:), allocatable, intent(out) :: W_dct_rcoords

          character(len=20)     :: filename


          !generate the filename for the file
          !containing the detectors
          filename = generate_dct_filename(timestep)
          
          
          !extract the detectors
          call read_detectors_from_file(
     $         filename,
     $         N_dct_rcoords,
     $         S_dct_rcoords,
     $         E_dct_rcoords,
     $         W_dct_rcoords)

        end subroutine get_detectors


        !generate the filename for the file containing
        !the detectors from the timestep
        function generate_dct_filename(index) result(filename)

          implicit none

          integer, intent(in) :: index
          character(len=20)   :: filename

          character(len=11) :: filename_format
          integer           :: index_format
   

          !determine the number of character to print the
          !timestep in the filename
          if(index.eq.0) then
             index_format = 1
          else
             index_format = floor(log10(real(index)))+1
          end if

          !determine the format to write the fielname
          write(filename_format, '(''(A9,I'',I1,'',A6)'')') index_format

          !determine the filename
          write(filename, filename_format)
     $         'detectors',
     $         index,
     $         '.curve'

        end function generate_dct_filename


        !check if two doubles are the same
        function is_test_validated(var,cst,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, nint(var*1e4)
             print *, nint(cst*1e4)
          end if
          
          test_validated=abs(
     $         nint(var*1e4)-
     $         nint(cst*1e4)).le.1
          
        end function is_test_validated

      end program test_dct_verify_x_symmetry
