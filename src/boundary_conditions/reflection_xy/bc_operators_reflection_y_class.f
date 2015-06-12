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
      module bc_operators_reflection_y_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use check_data_module, only :
     $       is_real_validated

        use errors_module, only :
     $       error_bc_section_type

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       bc_nodes_choice,
     $       N_edge_type,
     $       S_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       debug_initialize_bc_nodes,
     $       debug_real

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use reflection_xy_module, only :
     $       reflection_y_prefactor

        implicit none


        private
        public :: bc_operators_reflection_y


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> reflection boundary conditions in the y
        !> direction at the edges of the computational
        !> domain
        !
        !> @param apply_bc_on_nodes
        !> apply the reflection boundary conditions along
        !> the y direction at the edge of the
        !> computational domain for the field
        !
        !> @param apply_bc_on_nodes_nopt
        !> apply the reflection boundary conditions along
        !> the y direction at the edge of the
        !> computational domain for the field but only 
        !> on the bc_sections
        !---------------------------------------------------------------
        type, extends(bc_operators_default) :: bc_operators_reflection_y

          contains

          procedure, nopass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_nodes_nopt

        end type bc_operators_reflection_y


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes for a specific boundary section
        !> on the computational domain
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param p_model
        !> physical model
        !
        !>@param nodes
        !> governing variables at t
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(
     $       bc_section,
     $       t,x_map,y_map,nodes_tmp,
     $       p_model,
     $       nodes)

          implicit none

          integer    , dimension(4)       , intent(in)    :: bc_section
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          
          call apply_bc_on_nodes_nopt(
     $         [bc_section,0],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes for a specific boundary section
        !> on the sub-domain
        !
        !> @date
        !> 10_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed on sub-domain
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction on sub-domain
        !
        !>@param y_map
        !> coordinate map along the y-direction on sub-domain
        !
        !>@param nodes_tmp
        !> governing variables at t-dt on sub-domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param nodes
        !> governing variables at t on sub-domain
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes_nopt(
     $       bc_section,
     $       t,x_map,y_map,nodes_tmp,
     $       p_model,
     $       nodes)

          implicit none

          integer    , dimension(5)    , intent(in)    :: bc_section
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(in)    :: nodes_tmp
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          real(rkind)            :: s
          integer, dimension(ne) :: prefactor_y
          integer(ikind)         :: reflection_N
          integer(ikind)         :: reflection_S
          integer(ikind)         :: i,j
          integer                :: k

          s = t + x_map(1) + y_map(1) + nodes_tmp(1,1,1)

          prefactor_y  = reflection_y_prefactor(p_model)

          reflection_N = 2*bc_section(3)-1
          reflection_S = 2*bc_section(3)+2*bc_size-1


          select case(bc_section(1))

            case(N_edge_type)
               do k=1, ne
                  do j=bc_section(3), bc_section(3)+bc_size-1
                     do i=bc_section(2),bc_section(4)
                        
                        nodes(i,j,k) = prefactor_y(k)*
     $                       nodes(i,reflection_N-j,k)

                     end do
                  end do
               end do

            case(S_edge_type)
               do k=1, ne
                  do j=bc_section(3), bc_section(3)+bc_size-1
                     do i=bc_section(2),bc_section(4)
                        
                        nodes(i,j,k) = prefactor_y(k)*
     $                       nodes(i,reflection_S-j,k)

                     end do
                  end do
               end do

            case(SW_corner_type,SE_corner_type)
               do k=1, ne
                  do j=bc_section(3), bc_section(3)+bc_size-1
                     do i=bc_section(2),bc_section(2)+bc_size-1
                        
                        if(debug_initialize_bc_nodes) then

                           if(is_real_validated(nodes(i,reflection_S-j,k),debug_real)) then
                              nodes(i,j,k) = debug_real

                           else
                              nodes(i,j,k) = prefactor_y(k)*
     $                             nodes(i,reflection_S-j,k)
                           end if
                        else
                           nodes(i,j,k) = prefactor_y(k)*
     $                          nodes(i,reflection_S-j,k)
                        end if                        

                     end do
                  end do
               end do

            case(NW_corner_type,NE_corner_type)
               do k=1, ne
                  do j=bc_section(3), bc_section(3)+bc_size-1
                     do i=bc_section(2),bc_section(2)+bc_size-1
                        
                        if(debug_initialize_bc_nodes) then

                           if(is_real_validated(nodes(i,reflection_N-j,k),debug_real)) then
                              nodes(i,j,k) = debug_real

                           else
                              nodes(i,j,k) = prefactor_y(k)*
     $                             nodes(i,reflection_N-j,k)
                           end if
                        else
                           nodes(i,j,k) = prefactor_y(k)*
     $                          nodes(i,reflection_N-j,k)
                        end if

                     end do
                  end do
               end do

            case default
               call error_bc_section_type(
     $              'bc_operators_reflection_y_class',
     $              'apply_bc_on_nodes',
     $              bc_section(1))

          end select

        end subroutine apply_bc_on_nodes_nopt

      end module bc_operators_reflection_y_class
