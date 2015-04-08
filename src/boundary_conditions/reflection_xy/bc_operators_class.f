      !> @file
      !> class encapsulating subroutines to apply reflection boundary
      !> conditions at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply reflection boundary
      !> conditions at the edge of the computational domain
      !
      !> @date
      !> 07_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use pmodel_eq_class, only :
     $       pmodel_eq

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
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use reflection_xy_module, only :
     $       reflection_x_prefactor,
     $       reflection_y_prefactor

        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> reflection boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !>
        !> @param prefactor_x
        !> prefactor for the computation of the reflection
        !> boundary conditions along the x-direction
        !>
        !> @param prefactor_y
        !> prefactor for the computation of the reflection
        !> boundary conditions along the y-direction
        !>
        !> @param initialize
        !> initialize the prefactor_x and prefactor_y
        !> attributes of the boundary conditions
        !>
        !> @param apply_bc_on_nodes
        !> apply the reflection boundary conditions along
        !> the x and y directions at the edge of the
        !> computational domain for the field
        !>
        !> @param apply_bc_on_nodes_nopt
        !> apply the reflection boundary conditions along
        !> the x and y directions at the edge of the
        !> computational domain for the field but only 
        !> on the bc_sections
        !---------------------------------------------------------------
        type, extends(bc_operators_default) :: bc_operators

          integer, dimension(ne) :: prefactor_x
          integer, dimension(ne) :: prefactor_y

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
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine ini(this, p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model

          this%prefactor_x = reflection_x_prefactor(p_model)
          this%prefactor_y = reflection_y_prefactor(p_model)

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
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
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


          integer(ikind)         :: i,j
          integer                :: k
          

          !< compute the reflection b.c. in E and W boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                  this%prefactor_x(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                  this%prefactor_x(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do


          !< compute the reflection b.c. in N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  this%prefactor_y(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  this%prefactor_y(k)*nodes(i,ny-bc_size-j+1,k)
                   
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
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes_nopt(this,nodes,bc_sections)

          implicit none

          class(bc_operators)                        , intent(in)    :: this
          real(rkind)   , dimension(nx,ny,ne)        , intent(inout) :: nodes
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections


          integer(ikind) :: i,j
          integer(ikind) :: i_r,j_r
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

                           j_r = 2*j_min-(j+1)

                           do i=i_min,i_min+1
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_y(k)*nodes(i,j_r,k)
                              
                           end do
                        end do
                     end do

                  case(SE_corner_type,SW_corner_type)
                     do k=1,ne
                        do j=j_min,j_min+1

                           j_r = 2*j_min+3-j

                           do i=i_min,i_min+1
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_y(k)*nodes(i,j_r,k)
                              
                           end do
                        end do
                     end do

                  case(N_edge_type)

                     i_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_min+1

                           j_r = 2*j_min-(j+1)

                           do i=i_min,i_max
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_y(k)*nodes(i,j_r,k)
                              
                           end do
                        end do
                     end do

                  case(S_edge_type)

                     i_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_min+1

                           j_r = 2*j_min+3-j

                           do i=i_min,i_max
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_y(k)*nodes(i,j_r,k)
                              
                           end do
                        end do
                     end do

                  case(E_edge_type)

                     j_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_max

                           do i=i_min,i_min+1

                              i_r = 2*i_min-(i+1)
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_x(k)*nodes(i_r,j,k)
                              
                           end do
                        end do
                     end do

                  case(W_edge_type)

                     j_max = bc_sections(4,m)

                     do k=1,ne
                        do j=j_min,j_max

                           do i=i_min,i_min+1

                              i_r = 2*i_min+3-i
                              
                              nodes(i,j,k) = 
     $                             this%prefactor_x(k)*nodes(i_r,j,k)
                              
                           end do
                        end do
                     end do

                end select

             end do

          end if

        end subroutine apply_bc_on_nodes_nopt

      end module bc_operators_class
