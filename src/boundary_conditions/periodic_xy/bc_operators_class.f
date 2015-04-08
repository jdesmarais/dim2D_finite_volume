      !> @file
      !> class encapsulating subroutines to apply periodic boundary
      !> conditions at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the 
      !> gridpoints at the edge of the computational domain
      !
      !> @date
      !> 23_09_2013 - initial version               - J.L. Desmarais
      !> 14_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_constant, only :
     $       bc_nodes_choice

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,rkind

        use pmodel_eq_class, only :
     $       pmodel_eq
        
        use sd_operators_class, only :
     $       sd_operators

        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> periodic boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !
        !> @param period_x
        !> period along the x-direction
        !
        !> @param period_y
        !> period along the y-direction
        ! 
        !> @param initialize
        !> initialize the period_x and period_y
        !> attributes of the boundary conditions
        !
        !> @param apply_bc_on_nodes
        !> apply the periodic boundary conditions along the x and y
        !> directions at the edge of the computational domain
        !> for the field
        !>
        !> @param apply_bc_on_nodes_nopt
        !> apply the reflection boundary conditions along
        !> the x and y directions at the edge of the
        !> computational domain for the field but only 
        !> on the bc_sections
        !---------------------------------------------------------------
        type, extends(bc_operators_default) :: bc_operators

          integer(ikind) :: period_x
          integer(ikind) :: period_y

          contains

          procedure,   pass :: ini
          procedure,   pass :: apply_bc_on_nodes
          procedure,   pass :: apply_bc_on_nodes_nopt

        end type bc_operators


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param s
        !> spatial discretisation operators
        !
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine ini(this, p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model
          
          integer :: neq

          neq     = p_model%get_eq_nb()

          this%period_x = nx-2*bc_size
          this%period_y = ny-2*bc_size

          this%bc_type = [
     $         bc_nodes_choice,
     $         bc_nodes_choice,
     $         bc_nodes_choice,
     $         bc_nodes_choice]

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,nodes)

          implicit none

          class(bc_operators)             , intent(in)    :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          integer(ikind) :: i,j
          integer        :: k

          !<compute the east and west boundary layers
          !>without the north and south corners
          do k=1, ne
             !DEC$ IVDEP
             do j=1+bc_size, ny-bc_size
                !DEC$ UNROLL(2)
                do i=1, bc_size

                   nodes(i,j,k)=
     $                  nodes(i+this%period_x,j,k)
                   nodes(i+this%period_x+bc_size,j,k)= 
     $                  nodes(i+bc_size,j,k)
                   
                end do
             end do
          end do


          !<compute the south and north layers
          !>with the east and west corners
          do k=1, ne
             !DEC$ UNROLL(2)
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx

                   nodes(i,j,k)=
     $                  nodes(i,j+this%period_y,k)
                   nodes(i,j+this%period_y+bc_size,k)=
     $                  nodes(i,j+bc_size,k)

                end do
             end do
          end do

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 08_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
        !
        !>@param nodes
        !> grid-points
        !
        !>@param bc_sections
        !> boundary sections computed
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes_nopt(this,nodes,bc_sections)

          implicit none

          class(bc_operators)                        , intent(in)    :: this
          real(rkind)   , dimension(nx,ny,ne)        , intent(inout) :: nodes
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections


          integer(ikind) :: i,j
          integer(ikind) :: i_p,j_p
          integer        :: k
          integer        :: m


          integer(ikind) :: i_min,j_min
          integer(ikind) :: i_max,j_max


          if(allocated(bc_sections)) then
             
             do m=1, size(bc_sections,2)

                i_min = bc_sections(2,m)
                j_min = bc_sections(3,m)

                select case(bc_sections(1,m))

                  case(NE_corner_type,NW_corner_type)
                     do k=1,ne
                        do j=j_min,j_min+1

                           j_p = j-this%period_y

                           do i=i_min,i_min+1
                              
                              nodes(i,j,k) = nodes(i,j_p,k)
                              
                           end do
                        end do
                     end do

                  case(SE_corner_type,SW_corner_type)
                     do k=1,ne
                        do j=j_min,j_min+1

                           j_p = j+this%period_y

                           do i=i_min,i_min+1
                              
                              nodes(i,j,k) = nodes(i,j_p,k)
                              
                           end do
                        end do
                     end do

                  case(N_edge_type)

                     i_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_min+1

                           j_p = j-this%period_y

                           do i=i_min,i_max
                              
                              nodes(i,j,k) = nodes(i,j_p,k)
                              
                           end do
                        end do
                     end do

                  case(S_edge_type)

                     i_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_min+1

                           j_p = j+this%period_y

                           do i=i_min,i_max
                              
                              nodes(i,j,k) = nodes(i,j_p,k)
                              
                           end do
                        end do
                     end do

                  case(E_edge_type)

                     j_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_max

                           do i=i_min,i_min+1

                              i_p = i-this%period_x
                              
                              nodes(i,j,k) = nodes(i_p,j,k)
                              
                           end do
                        end do
                     end do

                  case(W_edge_type)

                     j_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_max

                           do i=i_min,i_min+1

                              i_p = i+this%period_x
                              
                              nodes(i,j,k) = nodes(i_p,j,k)
                              
                           end do
                        end do
                     end do

                end select

             end do

          end if

        end subroutine apply_bc_on_nodes_nopt

      end module bc_operators_class
