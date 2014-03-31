      module bf_layer_abstract_class

        use parameters_bf_layer, only : no_pt, interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: bf_layer_abstract


        !localization   : N,S,E,W,N_E,N_W,S_E,S_W
        !alignment      : position of the buffer layer compared
        !                 to the interior nodes table
        !                   - alignment(1,1) : i_min
        !                   - alignment(1,2) : i_max
        !                   - alignment(2,1) : j_min
        !                   - alignment(2,2) : j_max
        !
        !initialize     : initialize the localization of the buffer
        !                 layer compared to the interior table
        !allocate_nodes : allocate the nodes for the buffer layer
        type :: bf_layer_abstract

          integer                        :: localization
          integer(ikind), dimension(2,2) :: alignment
          
          real(rkind)   , dimension(:,:,:), allocatable :: nodes
          integer       , dimension(:,:)  , allocatable :: grdpts_id

          contains

          procedure, pass :: ini
          procedure, pass :: allocate_bf_layer
          procedure, pass :: print_nodes
          procedure, pass :: print_grdpts_id
          procedure, pass :: print_sizes

        end type bf_layer_abstract


        contains


        subroutine ini(this,localization)

          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer(ikind)          , intent(in)    :: localization
          
          this%localization = localization

        end subroutine ini
        

        subroutine allocate_bf_layer(this, alignment, nodes)

          implicit none
          
          class(bf_layer_abstract)        , intent(inout) :: this
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes

          integer(ikind), dimension(:,:)  , allocatable   :: list_new_grdpts

          !set the alignment between the interior table
          !and the buffer layer
          this%alignment = alignment

          !allocate the number of gridpoints needed for
          !the buffer layer,
          !copy gridpoints from the interior table,
          !and create the list of new gridpoints
          !that should be computed
          call initial_grdpt_copy_from_interior(
     $         this,
     $         nodes,
     $         list_new_grdpts)

          !compute new gridpoints
          call compute_new_gridpoints(this, list_new_grdpts)

          !update interior points list

          !update bc points list

        end subroutine allocate_bf_layer


        subroutine initial_grdpt_copy_from_interior(
     $     this,
     $     nodes,
     $     list_new_grdpts)

          implicit none

          class(bf_layer_abstract)                     , intent(inout) :: this
          real(rkind)   , dimension(:,:,:)             , intent(in)    :: nodes
          integer(ikind), dimension(:,:)  , allocatable, intent(out)   :: list_new_grdpts
          
          integer(ikind), dimension(3) :: size_nodes
          integer(ikind) :: i,j,k
          integer(ikind) :: i_match, j_match
          integer(ikind) :: nb_new_grdpts


          !allocate nodes for buffer layer
          select case(this%localization)
          
            case(N,S)
               size_nodes(1) = this%alignment(1,2) - this%alignment(1,1) + 2*bc_size
               size_nodes(2) = 2*bc_size+1

            case(E,W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = this%alignment(2,2) - this%alignment(2,1) + 2*bc_size

            case(N_E,N_W,S_E,S_W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = 2*bc_size+1

            case default
               print *, 'bf_layer_abstract_class'
               print *, 'initial_grdpt_copy_from_interior'
               stop 'localization not recognized'

          end select

          size_nodes(3) = ne
          allocate(this%nodes(size_nodes(1),size_nodes(2),size_nodes(3)))

          select case(this%localization)

            case(N)

               !copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 2*bc_size+1
               end do

            case(S)               

               !copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j+1,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               
            case(E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 2*bc_size+1
                  list_new_grdpts(j,2) = j
               end do

            case(W)

               !copy of grid points from the interior
               i_match = 0
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i+1,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do
               
               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 1
                  list_new_grdpts(j,2) = j
               end do

            case(N_E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size + 1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, 2*bc_size
                  list_new_grdpts(j,1) = 2*bc_size+1
                  list_new_grdpts(j,2) = j
               end do
               do i=1, 2*bc_size+1
                  list_new_grdpts(2*bc_size+i,1) = i
                  list_new_grdpts(2*bc_size+i,2) = 2*bc_size+1
               end do
               

            case(N_W)

               !copy of grid points from the interior
               i_match = 0
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size + 1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, 2*bc_size
                  list_new_grdpts(j,1) = 1
                  list_new_grdpts(j,2) = j
               end do
               do i=1, 2*bc_size+1
                  list_new_grdpts(2*bc_size+i,1) = i
                  list_new_grdpts(2*bc_size+i,2) = 2*bc_size+1
               end do

            case(S_E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size+1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, 2*bc_size+1
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               do j=2, 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,1) = 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,2) = j
               end do
               

            case(S_W)

               !copy of grid points from the interior
               i_match = 0
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size+1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, 2*bc_size+1
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               do j=2, 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,1) = 1
                  list_new_grdpts(2*bc_size+j,2) = j
               end do

          end select


          !identification of the buffer layer gridpoints
          allocate(this%grdpts_id(size_nodes(1),size_nodes(2)))
          do j=1, size_nodes(2)
             do i=1, size_nodes(1)
                this%grdpts_id(i,j) = bc_pt
             end do
          end do

          select case(this%localization)

            case(N,S)
               j = bc_size+1
               do i=bc_size+1, size_nodes(1)-bc_size
                  this%grdpts_id(i,j) = interior_pt
               end do

            case(E,W)
               i = bc_size+1
               do j=bc_size+1, size_nodes(2)-bc_size
                  this%grdpts_id(i,j) = interior_pt
               end do

            case(N_E,N_W,S_E,S_W)
               i = bc_size+1
               j = bc_size+1
               this%grdpts_id(i,j) = interior_pt

          end select

        end subroutine initial_grdpt_copy_from_interior

      
        subroutine compute_new_gridpoints(this, list_new_gridpts)

          implicit none

          class(bf_layer_abstract)      , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: list_new_gridpts
          
          integer(ikind) :: i,j,k,l

          do l=1, size(list_new_gridpts,1)
             
             i = list_new_gridpts(l,1)
             j = list_new_gridpts(l,2)
          
             do k=1, ne
                this%nodes(i,j,k) = 3
             end do

          end do

        end subroutine compute_new_gridpoints

      
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

      end module bf_layer_abstract_class
