      !> @file
      !> module encapsulating the subroutines related to the
      !> exchange of data between the interior and the buffer
      !> layers and between the different buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to allocation and reallocation of the
      !> buffer layer object
      !
      !> @date
      ! 04_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_exchange_module

        use bf_layer_abstract_class, only : bf_layer_abstract

        use parameters_bf_layer, only : exchange_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind


        implicit none


        private
        public :: first_exchange_with_interior,
     $            copy_interior_data_after_reallocation

        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the data contained in the
        !> main tables of the buffer layer (grdptid, nodes)
        !> using the data of the interior domain
        !> the grid points that can not be initialized by the
        !> exchange are identified and their position is saved
        !> in an additional table
        !
        !> @date
        !> 04_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param nodes
        !> table encapsulating the data of the internal
        !> grid points
        !
        !>@param list_new_grdpts
        !> table encapsulating the coordinates of the new grid points
        !--------------------------------------------------------------
        subroutine first_exchange_with_interior(
     $       this,
     $       nodes,
     $       list_new_grdpts)

          implicit none

          class(bf_layer_abstract)                     , intent(inout) :: this
          real(rkind)   , dimension(:,:,:)             , intent(in)    :: nodes
          integer(ikind), dimension(:,:)  , allocatable, intent(out)   :: list_new_grdpts
          

          integer(ikind) :: i,j,k
          integer(ikind) :: i_match, j_match
          integer(ikind) :: nb_new_grdpts


          !< copy of the interior grid points
          !> and identification of the new grid
          !> points
          select case(this%localization)

            case(N)

               !< copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 2*bc_size+1
               end do

            case(S)               

               !< copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j+1,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               
            case(E)

               !< copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 2*bc_size+1
                  list_new_grdpts(j,2) = j
               end do

            case(W)

               !< copy of grid points from the interior
               i_match = 0
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i+1,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do
               
               !< determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 1
                  list_new_grdpts(j,2) = j
               end do

            case(N_E)

               !< copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
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

               !< copy of grid points from the interior
               i_match = 0
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
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

               !< copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
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

               !< copy of grid points from the interior
               i_match = 0
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< determination of the list of new gridpoints
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

        end subroutine first_exchange_with_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine copying the data from the interior domain
        !> to the reallocated table
        !
        !> @date
        !> 22_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !
        !>@param nodes
        !> table encapsulating the data of the internal
        !> grid points
        !--------------------------------------------------------------
        subroutine copy_interior_data_after_reallocation(
     $       this,
     $       nodes,
     $       border_changes)

          implicit none

          class(bf_layer_abstract)          , intent(inout) :: this
          real(rkind)   , dimension(:,:,:)  , intent(in)    :: nodes
          integer(ikind), dimension(2,2)    , intent(in)    :: border_changes

          integer(ikind) :: i,j, i_match, j_match
          integer        :: k

          !< copy of the interior grid points
          !> and identification of the new grid
          !> points
          select case(this%localization)

            case(N)

               !< copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, -border_changes(1,1)
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               do k=1, ne
                  do j=1, 2*bc_size
                     do i=size(this%nodes,1)-border_changes(1,2)+1, size(this%nodes,1)
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< update the exchanged gridpoints
               do j=1, bc_size
                  do i=1, -border_changes(1,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

               do j=1, bc_size
                  do i=size(this%grdpts_id,1)-border_changes(1,2)+1, size(this%grdpts_id,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

            case(S)               

               !< copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = size(this%nodes,2)-(2*bc_size)
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, -border_changes(1,1)
                        this%nodes(i,j_match+j,k) = nodes(i_match+i,j,k)
                     end do
                  end do
               end do

               do k=1, ne
                  do j=1, 2*bc_size
                     do i=size(this%nodes,1)-border_changes(1,2)+1, size(this%nodes,1)
                        this%nodes(i,j_match+j,k) = nodes(i_match+i,j,k)
                     end do
                  end do
               end do

               !< update the exchanged gridpoints
               do j=size(this%grdpts_id,2)-bc_size+1, size(this%grdpts_id,2)
                  do i=1, -border_changes(1,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

               do j=size(this%grdpts_id,2)-bc_size+1, size(this%grdpts_id,2)
                  do i=size(this%grdpts_id,1)-border_changes(1,2)+1, size(this%grdpts_id,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do
               
            case(E)

               !< copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, -border_changes(2,1)
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               do k=1, ne
                  do j=size(this%nodes,2)-border_changes(2,2)+1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !< update the exchanged gridpoints
               do j=1, -border_changes(2,1)
                  do i=1, bc_size
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

               do j=size(this%grdpts_id,2)-border_changes(2,2)+1, size(this%grdpts_id,2)
                  do i=1, bc_size
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

            case(W)

               !< copy of grid points from the interior
               i_match = size(this%nodes,1)-(2*bc_size)
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, -border_changes(2,1)
                     do i=1, 2*bc_size
                        this%nodes(i_match+i,j,k) = nodes(i,j_match+j,k)
                     end do
                  end do
               end do
               
               do k=1, ne
                  do j=size(this%nodes,2)-border_changes(2,2)+1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i_match+i,j,k) = nodes(i,j_match+j,k)
                     end do
                  end do
               end do

               !< update the exchanged gridpoints
               do j=1, -border_changes(2,1)
                  do i=size(this%grdpts_id,1)-bc_size+1, size(this%grdpts_id,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

               do j=size(this%grdpts_id,2)-border_changes(2,2)+1, size(this%grdpts_id,2)
                  do i=size(this%grdpts_id,1)-bc_size+1, size(this%grdpts_id,1)
                     this%grdpts_id(i,j) = exchange_pt
                  end do
               end do

          end select

        end subroutine copy_interior_data_after_reallocation

      end module bf_layer_exchange_module
