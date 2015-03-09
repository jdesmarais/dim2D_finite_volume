      module bf_interface_print_class

        use bf_interface_basic_class, only :
     $     bf_interface_basic

        type, extends(bf_interface_basic) :: bf_interface_print

          contains

          procedure, pass :: ini
          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

        end type bf_interface_print


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer/interior domain interface
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        ! the data between them
        !--------------------------------------------------------------
        subroutine ini(this,interior_x_map,interior_y_map)

          implicit none

          class(bf_interface_print) , intent(inout) :: this
          real(rkind), dimension(nx), intent(in)    :: interior_x_map
          real(rkind), dimension(ny), intent(in)    :: interior_y_map


          !initialize the references to the mainlayer pointers
          !and the mainlayer interfaces
          call this%bf_interface_basic%ini(
     $         interior_x_map,
     $         interior_y_map)


          !create a netcdf file with the grdpts_id for the interior
          !computational domain
          call print_interior_grdpts_id_on_netcdf(
     $         'interior_grdpts_id.nc',
     $         interior_x_map,
     $         interior_y_map)


        end subroutine ini


               !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the content of the interface on external binary files
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param suffix_nodes
       !> suffix for the name of the files storing the nodes content
       !
       !>@param suffix_grdid
       !> suffix for the name of the files storing the grdpts_id content
       !
       !>@param suffix_sizes
       !> suffix for the name of the files storing the profile of the 
       !> nodes and grdpts_id arrays
       !
       !>@param suffix_nb_sublayers_max
       !> suffix for the name of the files storing the maximum number of
       !> sublayers per main buffer layer
       !--------------------------------------------------------------
       subroutine print_binary(
     $     this,
     $     suffix_x_map,
     $     suffix_y_map,
     $     suffix_nodes,
     $     suffix_grdid,
     $     suffix_sizes,
     $     suffix_nb_sublayers_max)

         implicit none

         class(bf_interface_print), intent(in) :: this
         character(*)             , intent(in) :: suffix_x_map
         character(*)             , intent(in) :: suffix_y_map
         character(*)             , intent(in) :: suffix_nodes
         character(*)             , intent(in) :: suffix_grdid
         character(*)             , intent(in) :: suffix_sizes
         character(*)             , intent(in) :: suffix_nb_sublayers_max

         integer           :: i
         integer           :: nb_sublayers_max
         character(len=18) :: filename_format
         character(len=28) :: nb_sublayers_filename

                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files and 
         !determine the maximum number of sublayers
         nb_sublayers_max = 0
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_binary(
     $              suffix_x_map,
     $              suffix_y_map,
     $              suffix_nodes,
     $              suffix_grdid,
     $              suffix_sizes)

               nb_sublayers_max = max(
     $              nb_sublayers_max,
     $              this%mainlayer_pointers(i)%get_nb_sublayers())
            end if            
         end do


         !print the maximum number of sublayers in an output
         !binary file
         write(filename_format,
     $        '(''(A12,A'',I2,'')'')')
     $        len(suffix_nb_sublayers_max)

         write(nb_sublayers_filename, filename_format)
     $        'sublayers_nb',
     $        suffix_nb_sublayers_max

         call print_nb_sublayers_max(
     $        nb_sublayers_filename, nb_sublayers_max)

        end subroutine print_binary


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the maxmimum number of sublayers per main layer
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param filename
       !> name of the file for storing the maximum number of sublayer
       !> per mainlayer
       !
       !>@param nb_sublayers
       !> maximum number of sublayers per main layer
       !--------------------------------------------------------------
       subroutine print_nb_sublayers_max(filename, nb_sublayers)

          implicit none

          character(*), intent(in) :: filename
          integer     , intent(in) :: nb_sublayers

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers_max


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the content of the interface on external netcdf files
       !
       !> @date
       !> 11_07_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param timestep_written
       !> integer identifying the timestep written
       !
       !>@param name_var
       !> table with the short name for the governing variables saved
       !> in the netcdf file
       !
       !>@param longname_var
       !> table with the long name for the governing variables saved
       !> in the netcdf file
       !
       !>@param unit_var
       !> table with the units of the governing variables saved
       !> in the netcdf file
       !
       !>@param time
       !> time corresponding to the data for the grdpts and the nodes
       !--------------------------------------------------------------
       subroutine print_netcdf(
     $     this,
     $     timestep_written,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     time)

         implicit none

         class(bf_interface_print)  , intent(in) :: this
         integer                    , intent(in) :: timestep_written
         character(*), dimension(ne), intent(in) :: name_var
         character(*), dimension(ne), intent(in) :: longname_var
         character(*), dimension(ne), intent(in) :: unit_var
         real(rkind)                , intent(in) :: time

         integer :: i
                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_netcdf(
     $              timestep_written,
     $              name_var,
     $              longname_var,
     $              unit_var,
     $              time)

            end if
         end do

        end subroutine print_netcdf

      end module bf_interface_print_class
