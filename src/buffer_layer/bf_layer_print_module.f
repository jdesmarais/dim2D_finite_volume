      !> @file
      !> module encapsulating the subroutines related to printing
      !> the main attributes of the buffer layer object
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to printing the main attributes of the
      !> buff erlayer object
      !
      !> @date
      ! 07_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_print_module

        use bf_layer_abstract_class, only : bf_layer_abstract

        implicit none

        private
        public :: print_nodes,
     $            print_grdpts_id,
     $            print_sizes


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the nodes table in a binary
        !> file
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param filename
        !> name of the binary file where the nodes are
        !> written
        !--------------------------------------------------------------
        subroutine print_nodes(this,filename)

          implicit none

          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

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
              write(unit=1, iostat=ios) this%nodes
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the grdpt_id table in a binary
        !> file
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param filename
        !> name of the binary file where the grdpt_id are
        !> written
        !--------------------------------------------------------------
        subroutine print_grdpts_id(this,filename)

          implicit none
          
          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

          integer :: ios

          open(unit=2,
     $          file=filename,
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

        end subroutine print_grdpts_id
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the sizes of the main table in
        !> a binary file
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param filename
        !> name of the binary file where the sizes are
        !> written
        !--------------------------------------------------------------
        subroutine print_sizes(this, filename)

          implicit none
          
          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

          integer :: ios

          open(unit=2,
     $          file=filename,
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

        end subroutine print_sizes

      end module bf_layer_print_module
