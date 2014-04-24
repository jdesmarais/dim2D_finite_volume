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

        use parameters_kind, only : rkind

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
        subroutine print_nodes(nodes,filename)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          character(*)                 , intent(in) :: filename

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
              write(unit=1, iostat=ios) nodes
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
        subroutine print_grdpts_id(grdpts_id,filename)

          implicit none
          
          integer, dimension(:,:) , intent(in) :: grdpts_id
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
             write(unit=2, iostat=ios) grdpts_id
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
        subroutine print_sizes(nodes, alignment, filename)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer    , dimension(2,2)  , intent(in) :: alignment
          character(*)                 , intent(in) :: filename

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
     $            size(nodes,1),
     $            size(nodes,2),
     $            size(nodes,3),
     $            alignment(1,1),
     $            alignment(1,2),
     $            alignment(2,1),
     $            alignment(2,2)
             close(unit=2)
          else
             stop 'file opening pb'
          end if

        end subroutine print_sizes

      end module bf_layer_print_module
