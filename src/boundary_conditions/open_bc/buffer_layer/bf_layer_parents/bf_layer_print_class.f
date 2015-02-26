      !> @file
      !> module encapsulating the bf_layer_print object.
      !> bf_layer_basic enhanced with i/o intrinsic procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_basic enhanced with i/o intrinsic procedures
      !
      !> @date
      ! 23_02_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_print_class

        use bf_layer_basic_class, only :
     $       bf_layer_basic

        use bf_layer_nf90_operators_module, only :
     $       print_bf_layer_on_netcdf

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       rkind

        private
        public :: bf_layer_print


        !> @class bf_layer_print
        !> class encapsulating the i/o intrinsic procedures
        !> for the buffer layer object
        !
        !> @param print_binary
        !> print the nodes and the grdpts_id attributes
        !> as well as the size of the previous tables in
        !> output binary files
        !
        !> @param print_netcdf
        !> print the nodes and the grdpts_id attributes
        !> on a netcdf file
        !-------------------------------------------------------------
        type, extends(bf_layer_basic) :: bf_layer_print

          contains

          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

        end type bf_layer_print

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the nodes and the grdpts_id attributes
        !> as well as the size of the previous tables in
        !> output binary files
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param filename_nodes
        !> name of the output file for the data of the nodes attribute
        !
        !>@param filename_grdpts_id
        !> name of the output file for the data of the grdpts_id attribute
        !
        !>@param filename_sizes
        !> name of the output file for the sizes of the nodes and
        !> grpdts_id attributes
        !--------------------------------------------------------------
        subroutine print_binary(
     $     this,
     $     filename_x_map,
     $     filename_y_map,
     $     filename_nodes,
     $     filename_grdpts_id,
     $     filename_sizes)

          implicit none

          class(bf_layer_print), intent(in) :: this
          character(*)         , intent(in) :: filename_x_map
          character(*)         , intent(in) :: filename_y_map
          character(*)         , intent(in) :: filename_nodes
          character(*)         , intent(in) :: filename_grdpts_id
          character(*)         , intent(in) :: filename_sizes

          integer :: ios
          
          !x_map
          open(unit=2,
     $         file=filename_x_map,
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%x_map
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          !y_map
          open(unit=2,
     $         file=filename_y_map,
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%y_map
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          !nodes
          open(unit=3,
     $          file=filename_nodes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=3, iostat=ios) this%nodes
             close(unit=3)
          else
             stop 'file opening pb'
          end if

          !grdpts_id
          open(unit=2,
     $          file=filename_grdpts_id,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%grdpts_id
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          !sizes
          open(unit=2,
     $          file=filename_sizes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios)
     $            size(this%nodes,1),
     $            size(this%nodes,2),
     $            size(this%nodes,3),
     $            this%alignment(1,1),
     $            this%alignment(1,2),
     $            this%alignment(2,1),
     $            this%alignment(2,2)
             close(unit=2)
          else
             stop 'file opening pb'
          end if

        end subroutine print_binary


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the nodes and the grdpts_id attributes
        !> on a netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param filename
        !> name of the output file for the grdpts_id of the nodes
        !> attributes
        !
        !>@param name_var
        !> physical model used to know the name of the governing
        !> variables when writing the netcdf file
        !
        !>@param bf_order
        !> identification of the buffer layer in the main layer to
        !> determine the title of the netcdf file
        !
        !>@param x_min_interior
        !> x-coordinate corresponding to the grid point next to the
        !> left boundary layer in the interior domain
        !
        !>@param y_min_interior
        !> y-coordinate corresponding to the grid point next to the
        !> lower boundary layer in the interior domain
        !
        !>@param dx
        !> buffer layer grid size along the x-direction
        !
        !>@param dy
        !> buffer layer grid size along the y-direction
        !
        !>@param time
        !> time corresponding to the buffer layer data
        !--------------------------------------------------------------
        subroutine print_netcdf(
     $     this,
     $     filename,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     bf_order,
     $     time)

          implicit none

          class(bf_layer_print)      , intent(inout) :: this
          character(*)               , intent(in)    :: filename
          character(*), dimension(ne), intent(in)    :: name_var
          character(*), dimension(ne), intent(in)    :: longname_var
          character(*), dimension(ne), intent(in)    :: unit_var
          integer                    , intent(in)    :: bf_order
          real(rkind)                , intent(in)    :: time

          call print_bf_layer_on_netcdf(
     $         filename,
     $         name_var, longname_var, unit_var,
     $         this%localization, bf_order,
     $         this%grdpts_id,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         time)

        end subroutine print_netcdf

      end module bf_layer_print_class
