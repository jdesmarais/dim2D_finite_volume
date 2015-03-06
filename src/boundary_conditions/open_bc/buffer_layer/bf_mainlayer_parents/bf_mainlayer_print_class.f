      !> @file
      !> bf_mainlayer_basic augmented with print procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_mainlayer_basic augmented with print procedures
      !
      !> @date
      ! 06_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_print_class

        use bf_mainlayer_basic_class, only :
     $       bf_mainlayer_basic
      
        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: bf_mainlayer_print
        
        
        !> @class bf_mainlayer_print
        !> bf_mainlayer_basic augmented with print procedures
        !
        !> @param print_binary
        !> print the content of the mainlayer on binary files, one
        !> file per buffer layer
        !
        !> @param print_netcdf
        !> print the content of the mainlayer on netcdf files, one
        !> file per buffer layer
        !---------------------------------------------------------------
        type, extends(bf_mainlayer_basic) :: bf_mainlayer_print

          contains

          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

        end type bf_mainlayer_print


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param suffix_nodes
        !> suffix for the name of the output files storing the nodes
        !> of the bf_sublayers
        !
        !> @param suffix_grdid
        !> suffix for the name of the output files storing the grdpts_id
        !> of the bf_sublayers
        !
        !> @param suffix_sizes
        !> suffix for the name of the output files storing the sizes
        !> of the bf_sublayers        
        !--------------------------------------------------------------
        subroutine print_binary(
     $     this,
     $     suffix_x_map,
     $     suffix_y_map,
     $     suffix_nodes,
     $     suffix_grdid,
     $     suffix_sizes)

          implicit none

          class(bf_mainlayer_print), intent(in) :: this
          character(*)             , intent(in) :: suffix_x_map
          character(*)             , intent(in) :: suffix_y_map
          character(*)             , intent(in) :: suffix_nodes
          character(*)             , intent(in) :: suffix_grdid
          character(*)             , intent(in) :: suffix_sizes


          character(2), dimension(4) :: bf_layer_char
          integer                    :: i
          character(len=14)          :: filename_format
          character(len=30)          :: sizes_filename
          character(len=30)          :: x_map_filename
          character(len=30)          :: y_map_filename
          character(len=30)          :: nodes_filename
          character(len=30)          :: grdid_filename
          type(bf_sublayer), pointer :: current_sublayer


          bf_layer_char = ['N_','S_','E_','W_']


          current_sublayer => this%head_sublayer


          !go through the chained list and write the content
          !of each element
          do i=1, this%nb_sublayers

             !determine the name of the output files
             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_nodes)

             write(x_map_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_x_map

             write(y_map_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_y_map

             write(nodes_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_nodes


             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_grdid)
             write(grdid_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_grdid

             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_sizes)
             write(sizes_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_sizes


             !write the content of the sublayer corresponding
             !to the index i
             call current_sublayer%print_binary(
     $            x_map_filename,
     $            y_map_filename,
     $            nodes_filename,
     $            grdid_filename,
     $            sizes_filename)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_binary


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on a netcdf output file
        !
        !> @date
        !> 11_07_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param timestep_written
        !> integer identifying the timestep written
        !
        !> @param name_var
        !> table with the short name for the governing variables saved
        !> in the netcdf file
        !
        !> @param longname_var
        !> table with the long name for the governing variables saved
        !> in the netcdf file
        !
        !> @param unit_var
        !> table with the units of the governing variables saved
        !> in the netcdf file
        !
        !> @param x_min_interior
        !> x-coordinate of interior grid point next to the left
        !> boundary layer
        !
        !> @param y_min_interior
        !> y-coordinate of interior grid point next to the lower
        !> boundary layer
        !
        !> @param dx
        !> grid size along the x-coordinate
        !
        !> @param dy
        !> grid size along the y-coordinate
        !
        !> @param time
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

          class(bf_mainlayer_print)  , intent(in) :: this
          integer                    , intent(in) :: timestep_written
          character(*), dimension(ne), intent(in) :: name_var
          character(*), dimension(ne), intent(in) :: longname_var
          character(*), dimension(ne), intent(in) :: unit_var
          real(rkind)                , intent(in) :: time


          character(2), dimension(4) :: bf_layer_char
          integer                    :: size_timestep
          integer                    :: size_bf_order
          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i
          character(len=16)          :: filename_format
          character(len=30)          :: filename


          bf_layer_char = ['N_','S_','E_','W_']


          !determine the size needed to write the timestep
          if (timestep_written.eq.0) then
             size_timestep  = 1
          else
             size_timestep  = floor(
     $            log(float(timestep_written))/log(10.))+1
          end if

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and write the content
          !of each element
          do i=1, this%nb_sublayers

             !determine the size needed to write the bf_order
             size_bf_order = floor(log(float(i))/log(10.))+1

             !determine the name of the netcdf file
             write(filename_format,
     $            '(''(A2,I'',I1,'',A1,I'',I1,'',A3)'')')
     $            size_bf_order, size_timestep

             write(filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_',
     $            timestep_written, '.nc'

             !write the content of the sublayer corresponding
             !to the index i
             call current_sublayer%print_netcdf(
     $            trim(filename),
     $            name_var,
     $            longname_var,
     $            unit_var,
     $            i,
     $            time)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_netcdf       

      end module bf_mainlayer_print_class
