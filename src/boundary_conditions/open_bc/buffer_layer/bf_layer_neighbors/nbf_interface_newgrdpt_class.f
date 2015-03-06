      !> @file
      !> nbf_interface augmented with functions computing new grid points
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object inheriting from nbf_interface
      !> to encasulate the computation of new grid points
      !
      !> @date
      ! 21_11_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_interface_newgrdpt_class

        use bf_layer_extract_module, only :
     $       get_grdpts_id_from_interior

        use bf_newgrdpt_extract_module, only :
     $       get_interior_data_for_newgrdpt,
     $       are_intermediate_newgrdpt_data_needed,
     $       get_x_map_for_newgrdpt,
     $       get_y_map_for_newgrdpt

        use bf_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure

        use bf_suspicious_bc_interior_pt_module, only :
     $       verify_if_all_grdpts_exist

        use bf_bc_crenel_module, only :
     $       is_temp_array_needed_for_bc_crenel,
     $       detect_and_curb_bc_crenels

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_interface_sync_class, only :
     $       nbf_interface_sync

        use parameters_bf_layer, only :
     $       BF_SUCCESS, bc_pt

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        private
        public ::
     $       nbf_interface_newgrdpt


        !>@class nbf_interface-newgrdpt
        !> object augmenting the nbf_interface with functions 
        !> computing the new gridpoints
        !
        !> @param update_grdpts_after_increase
        !> turn the grdpts_id identified by general coordinates
        !> from bc_interior_pt to interior_pt and reallocate the
        !> buffer layer such that the neighboring points around it
        !> are allocated. Then compute these new grid points
        !
        !> @param compute_newgrdpt
        !> compute the new grid point in the buffer layer resulting
        !> from the computational domain extension
        !--------------------------------------------------------------
        type, extends(nbf_interface_sync) :: nbf_interface_newgrdpt

          contains

          procedure, pass :: get_data_for_newgrdpt
          procedure, pass :: get_grdpts_id_part
          procedure, pass :: set_grdpts_id_part

          procedure, pass :: update_bf_grdpts_after_increase
          procedure, pass :: update_grdpts_around_new_interior_pt
          procedure, pass :: compute_newgrdpt

        end type nbf_interface_newgrdpt


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> complete the temporary arrays containing the grid point ID
        !> and the nodes needed for the computation of new grid points
        !> using the grid points of the neighboring buffer layers
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the buffer layer position
        !
        !>@param bf_neighbor_id
        !> integer identifying the type of neighbor
        !
        !>@param tmp_grdpts_id0
        !> temporary array with the grid point ID at t=t-dt
        !
        !>@param tmp_nodes0
        !> temporary array with the nodes at t=t-dt
        !
        !>@param tmp_nodes1
        !> temporary array with the nodes at t=t
        !
        !>@param gen_borders
        !> array with integers identifying the extent of the data
        !> extracted using general coordinates
        !--------------------------------------------------------------
        subroutine get_data_for_newgrdpt(
     $     this,
     $     bf_localization,
     $     bf_neighbor_id,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_borders)

          implicit none

          class(nbf_interface_newgrdpt)                             , intent(in)    :: this
          integer                                                   , intent(in)    :: bf_localization
          integer                                                   , intent(in)    :: bf_neighbor_id
          integer    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1)   , intent(inout) :: tmp_grdpts_id0
          real(rkind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes0
          real(rkind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes1
          integer(ikind), dimension(2,2)                            , intent(in)    :: gen_borders


          call this%nbf_links(bf_localization,bf_neighbor_id)%get_data_for_newgrdpt(
     $         tmp_grdpts_id0,
     $         tmp_nodes0,
     $         tmp_nodes1,
     $         gen_borders)

        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> substitute the content of the the temporary array
        !> containing the grid point ID at t by the content
        !> of grdpts_id matching the gen_borders from the
        !> neighboring buffer layer identified by
        !> (bf_localization,bf_neighbor_id)
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the buffer layer position
        !
        !>@param bf_neighbor_id
        !> integer identifying the type of neighbor
        !
        !>@param tmp_grdpts_id1
        !> temporary array with the grid point ID at t=t
        !
        !>@param gen_borders
        !> array with integers identifying the extent of the data
        !> extracted using general coordinates
        !--------------------------------------------------------------
        subroutine get_grdpts_id_part(
     $     this,
     $     bf_localization,
     $     bf_neighbor_id,
     $     tmp_grdpts_id1,
     $     gen_borders)

          implicit none

          class(nbf_interface_newgrdpt) , intent(in)    :: this
          integer                       , intent(in)    :: bf_localization
          integer                       , intent(in)    :: bf_neighbor_id
          integer       , dimension(:,:), intent(inout) :: tmp_grdpts_id1
          integer(ikind), dimension(2,2), intent(in)    :: gen_borders


          call this%nbf_links(bf_localization,bf_neighbor_id)%get_grdpts_id_part(
     $         tmp_grdpts_id1,
     $         gen_borders)

        end subroutine get_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> replace the content of the grdpts_id of the neighboring
        !> buffer layer identified by (bf_localization,bf_neighbor_id)
        !> matching the gen_borders by the the temporary array
        !> containing the grid point id
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the buffer layer position
        !
        !>@param bf_neighbor_id
        !> integer identifying the type of neighbor
        !
        !>@param tmp_grdpts_id1
        !> temporary array with the grid point ID at t=t
        !
        !>@param gen_borders
        !> array with integers identifying the extent of the data
        !> extracted using general coordinates
        !--------------------------------------------------------------
        subroutine set_grdpts_id_part(
     $     this,
     $     bf_localization,
     $     bf_neighbor_id,
     $     tmp_grdpts_id1,
     $     gen_borders)

          implicit none

          class(nbf_interface_newgrdpt) , intent(in) :: this
          integer                       , intent(in) :: bf_localization
          integer                       , intent(in) :: bf_neighbor_id
          integer       , dimension(:,:), intent(in) :: tmp_grdpts_id1
          integer(ikind), dimension(2,2), intent(in) :: gen_borders


          call this%nbf_links(bf_localization,bf_neighbor_id)%set_grdpts_id_part(
     $         tmp_grdpts_id1,
     $         gen_borders)

        end subroutine set_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> turn the grdpts_id identified by general coordinates
        !> from bc_interior_pt to interior_pt and reallocate the buffer
        !> layer such that the neighboring points around it are allocated.
        !> Then compute these new grid points.
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time 
        !
        !> @param dt
        !> time step
        !
        !> @param interior_x_map
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !
        !>@param selected_grdpts
        !> list containing the general coordinates of the grid points to be
        !> turned from bc_interior_pt to interior_pt
        !--------------------------------------------------------------
        subroutine update_bf_grdpts_after_increase(
     $       this,
     $       bf_sublayer_updated,
     $       p_model,
     $       t,dt,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       selected_grdpts)

          implicit none

          class(nbf_interface_newgrdpt)   , intent(in)    :: this
          type(bf_sublayer)               , intent(inout) :: bf_sublayer_updated
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          integer(ikind), dimension(:,:)  , intent(in)    :: selected_grdpts

          integer(ikind), dimension(2) :: match_table
          integer(ikind)               :: k
          integer(ikind)               :: i_prev
          integer(ikind)               :: j_prev
          integer(ikind)               :: i_center
          integer(ikind)               :: j_center


          match_table = bf_sublayer_updated%get_general_to_local_coord_tab()
          
          
          !we have a list of gridpoints that should be turned
          !from bc_interior_pt to interior_pt. For a point to be
          !considered an interior point, we need to make sure
          !that the grid points it needs to compute its time
          !derivatives are available
          !
          !previously, the nodes table of the buffer layer was
          !increased to take into account the space needed for
          !new gridpoints at the boundary
          !
          !in this function, we go through the list of gridpoint
          !whose neighbours need to be tested. As it is possible
          !to reduce the number of tests by considering the
          !previous gridpoint tested, there is a distinction
          !between k=1 and k>1 in the list of gridpoints
          !----------------------------------------------------

          !for the first grid point, there is no simplification
          !possible in checking the neighbours so the centers
          !are initialized to do as if the previous central point
          !checked was far away
          !----------------------------------------------------
          i_center = -match_table(1)+nx/2
          j_center = -match_table(2)+ny/2

          do k=1, size(selected_grdpts,2)

             !update the position of the gridpoint previously
             !tested
             i_prev   =   i_center
             j_prev   =   j_center
             i_center = - match_table(1) + selected_grdpts(1,k)
             j_center = - match_table(2) + selected_grdpts(2,k)

             !check its neighbours
             call update_grdpts_around_new_interior_pt(
     $            this,
     $            bf_sublayer_updated,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            i_prev, j_prev,
     $            i_center, j_center)

             !update the status of the central gridpoint
             !to interior_pt
             call bf_sublayer_updated%set_new_interior_pt(
     $            i_center,
     $            j_center)

             !check whether the neighboring bc_interior_pt
             !should be updated to interior_pt
             call finalize_grdpts_for_suspicious_bc_interior_pt(
     $            this,
     $            bf_sublayer_updated,
     $            [i_center,j_center],
     $            match_table)

             !check whether the update the grdpts_id lead to
             !a crenel of bc_pt
             call finalize_grdpts_for_bc_pt_crenel(
     $            this,
     $            bf_sublayer_updated,
     $            [i_center,j_center],
     $            match_table)

          end do

        end subroutine update_bf_grdpts_after_increase


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> update the grid points around the new interior_pt
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time 
        !
        !> @param dt
        !> time step
        !
        !> @param interior_x_map
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !
        !> @param i_prev
        !> x-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param j_prev
        !> y-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param i_center
        !> x-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param j_center
        !> y-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !---------------------------------------------------------------
        subroutine update_grdpts_around_new_interior_pt(
     $     this,
     $     bf_sublayer_updated,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     i_prev, j_prev,
     $     i_center, j_center)

          implicit none

          class(nbf_interface_newgrdpt)   , intent(in)    :: this
          type(bf_sublayer)               , intent(inout) :: bf_sublayer_updated
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          integer(ikind)                  , intent(in)    :: i_prev
          integer(ikind)                  , intent(in)    :: j_prev
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center

          integer(ikind) :: min_j, max_j
          integer(ikind) :: i,j

          logical :: compute_grdpt


          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)


          !check whether the grid points around the new interior_pt
          !exist: if they do not exist, they are computed
          !in any case, the status of the neighboring points may be
          !updated: from bc_pt -> bc_interior_pt or no_pt -> bc_pt
          do j=j_center-bc_size, j_prev - bc_size + min(j_center-j_prev+2*bc_size,-1)
             do i=i_center-bc_size,i_center+bc_size
                
                compute_grdpt = bf_sublayer_updated%update_neighboring_grdpt(
     $               i,j,i_center,j_center)

                if(compute_grdpt) then
                   call compute_newgrdpt(
     $                  this,
     $                  bf_sublayer_updated,
     $                  p_model,t,dt,
     $                  i,j,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes0,
     $                  interior_nodes1)
                end if

             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_center-bc_size, i_prev-bc_size+min(i_center-i_prev+2*bc_size,-1)

                compute_grdpt = bf_sublayer_updated%update_neighboring_grdpt(
     $               i,j,i_center,j_center)

                if(compute_grdpt) then
                   call compute_newgrdpt(
     $                  this,
     $                  bf_sublayer_updated,
     $                  p_model,t,dt,
     $                  i,j,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes0,
     $                  interior_nodes1)
                end if

             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_prev+bc_size+max(i_center-i_prev-2*bc_size,1),i_center+bc_size

                compute_grdpt = bf_sublayer_updated%update_neighboring_grdpt(
     $               i,j,i_center,j_center)

                if(compute_grdpt) then
                   call compute_newgrdpt(
     $                  this,
     $                  bf_sublayer_updated,
     $                  p_model,t,dt,
     $                  i,j,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes0,
     $                  interior_nodes1)
                end if

             end do
          end do

          do j=j_prev+bc_size+max(j_center-j_prev-2*bc_size,1), j_center+bc_size
             do i=i_center-bc_size,i_center+bc_size

                compute_grdpt = bf_sublayer_updated%update_neighboring_grdpt(
     $               i,j,i_center,j_center)

                if(compute_grdpt) then
                   call compute_newgrdpt(
     $                  this,
     $                  bf_sublayer_updated,
     $                  p_model,t,dt,
     $                  i,j,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes0,
     $                  interior_nodes1)
                end if

             end do
          end do


          !finalize the update of the neighboring gridpoints
          !grid points that do not need to be computed but whose
          !status is modified from bc_pt to bc_interior_pt
          call bf_sublayer_updated%finalize_neighboring_grdpts_update(
     $         i_prev, j_prev,
     $         i_center, j_center)

        end subroutine update_grdpts_around_new_interior_pt

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point of the buffer layer
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !> @param bf_sublayer_ptr
        !> buffer layer whose new grid point is computed
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time 
        !
        !> @param dt
        !> time step
        !
        !> @param i1
        !> x-index of the new grid point
        !
        !> @param j1
        !> y-index of the new grid point
        !
        !> @param interior_x_map
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !--------------------------------------------------------------
        subroutine compute_newgrdpt(
     $     this,
     $     bf_sublayer_updated,
     $     p_model,t,dt,
     $     i1,j1,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(nbf_interface_newgrdpt)   , intent(in)    :: this
          type(bf_sublayer)               , intent(inout) :: bf_sublayer_updated
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          integer(ikind)                  , intent(in)    :: i1
          integer(ikind)                  , intent(in)    :: j1
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2)   :: gen_coords
          logical                        :: tmp_needed
          integer(ikind), dimension(2,2) :: gen_borders

          integer    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1)    :: tmp_grdpts_id0
          real(rkind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne) :: tmp_nodes0
          real(rkind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne) :: tmp_nodes1

          integer                                    :: localization
          type(bf_newgrdpt)                          :: bf_newgrdpt_used
          integer                                    :: procedure_type
          integer                                    :: gradient_type
          integer(ikind), dimension(2,2)             :: bf_align
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_x_map
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_y_map
          real(rkind)   , dimension(ne)              :: new_grdpt
          integer                                    :: k


          !localization of the buffer layer
          localization = bf_sublayer_updated%get_localization()

          !compute the general coordinates of the new grid point
          match_table   = bf_sublayer_updated%get_general_to_local_coord_tab()
          gen_coords(1) = match_table(1) + i1
          gen_coords(2) = match_table(2) + j1
          
          !determine the extent of the data needed
          gen_borders(1,1) = gen_coords(1)-(bc_size+1)
          gen_borders(1,2) = gen_coords(1)+(bc_size+1)
          gen_borders(2,1) = gen_coords(2)-(bc_size+1)
          gen_borders(2,2) = gen_coords(2)+(bc_size+1)

          !determine whether the new grid point is at the interface
          !between buffer layers or at the interface with the interior
          !and so requires to create intermediate arrays gathering data
          !from the interior and the neighboring buffer layers
          tmp_needed = are_intermediate_newgrdpt_data_needed(
     $         localization,
     $         gen_borders)

          if(.not.tmp_needed) then
             tmp_needed = .not.(bf_sublayer_updated%does_previous_timestep_exist())
          end if


          !if intermediate data are required, data are gathered from
          !the interior, the buffer layer and its neighbors
          if(tmp_needed) then

             !gather the data from the interior domain
             call get_interior_data_for_newgrdpt(
     $            interior_nodes0,
     $            interior_nodes1,
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            gen_borders)

             !gather the data from the current buffer
             !layer
             call bf_sublayer_updated%get_data_for_newgrdpt(
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            gen_borders)

             !gather the data from the potential neighbors
             !of the buffer layer
             !- gather data from neighbor of type 1
             if(gen_borders(1,1).le.0) then
                
                call this%get_data_for_newgrdpt(
     $               localization,1,
     $               tmp_grdpts_id0,
     $               tmp_nodes0,
     $               tmp_nodes1,
     $               gen_borders)

             end if

             !- gather data from neighbor of type 2
             if(gen_borders(1,2).ge.(nx+1)) then

                call this%get_data_for_newgrdpt(
     $               localization,2,
     $               tmp_grdpts_id0,
     $               tmp_nodes0,
     $               tmp_nodes1,
     $               gen_borders)

             end if

             
             !determine the procedure and gradient type
             call get_newgrdpt_procedure(
     $            bc_size+2,bc_size+2,
     $            tmp_grdpts_id0,
     $            procedure_type,
     $            gradient_type)


             !compute the new grdpt
             ! - construct the bf_align array to localize
             !   the buffer layer
             bf_align(1,1) = gen_borders(1,1)+bc_size
             bf_align(1,2) = gen_borders(1,2)-bc_size
             bf_align(2,1) = gen_borders(2,1)+bc_size
             bf_align(2,2) = gen_borders(2,2)-bc_size

             ! - construct the bf_x_map array
             tmp_x_map = get_x_map_for_newgrdpt(interior_x_map, gen_borders)

             ! - construct the bf_y_map array
             tmp_y_map = get_y_map_for_newgrdpt(interior_y_map, gen_borders)

             ! - compute the new grdpt
             new_grdpt = bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            bf_align, tmp_x_map, tmp_y_map, tmp_nodes0,
     $            bf_align, tmp_x_map, tmp_y_map, tmp_nodes1,
     $            bc_size+2, bc_size+2,
     $            procedure_type, gradient_type)


             !set the new grdpt in the buffer layer
             do k=1,ne
                call bf_sublayer_updated%set_nodes_pt(i1,j1,k,new_grdpt(k))
             end do


          !otherwise, the grid point can be directly computed
          !from the buffer layer
          else

              call bf_sublayer_updated%compute_newgrdpt(
     $            p_model,
     $            t,dt,
     $            i1,j1)

          end if

        end subroutine compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the bc_interior_pt around the new
        !> interior_pt should be updated to interior_pt as well
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param newgrdpt_local_coords
        !> coordinates of the bc_interior_pt turned into interior_pt
        !> expressed as indices of the buffer layer bf_sublayer_updated
        !
        !> @param match_table
        !> array used to convert coordinates expressed in the general
        !> reference frame of the interior domain into local coordinates
        !> of the buffer layer bf_sublayer_updated
        !---------------------------------------------------------------
        subroutine finalize_grdpts_for_suspicious_bc_interior_pt(
     $     this,
     $     bf_sublayer_updated,
     $     newgrdpt_local_coords,
     $     match_table)

          implicit none

          class(nbf_interface_newgrdpt), intent(in)    :: this
          type(bf_sublayer)            , intent(inout) :: bf_sublayer_updated
          integer(ikind), dimension(2) , intent(in)    :: newgrdpt_local_coords
          integer(ikind), dimension(2) , intent(in)    :: match_table

          integer(ikind)                                             :: i,j
          logical                                                    :: has_a_bc_pt_neighbor
          logical                                                    :: ierror
          integer(ikind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1) :: tmp_grdpts_id


          do j=newgrdpt_local_coords(2)-bc_size, newgrdpt_local_coords(2)+bc_size
             do i=newgrdpt_local_coords(1)-bc_size, newgrdpt_local_coords(1)+bc_size

                !check whether there is a bc_interior_pt in the grid points
                !next to the grid points updated by the interior_pt
                if(bf_sublayer_updated%is_bc_interior_pt(i,j)) then
                
                   !check whether the neighboring bc_interior_pt is suspicious
                   !i.e. the grid point has no bc_pt as neighbor
                   has_a_bc_pt_neighbor = bf_sublayer_updated%has_a_bc_pt_neighbor(
     $                  [i,j],
     $                  ierror)

                   !if there were enough neighbors to check whether the
                   !neighboring bc_interior_pt is suspicious or not
                   if(ierror.eqv.BF_SUCCESS) then

                      if(.not.has_a_bc_pt_neighbor) then
                      
                         call update_bc_interior_pt_to_interior_pt(
     $                        this,
     $                        bf_sublayer_updated,
     $                        [i,j],
     $                        match_table)

                      end if

                   !if there were not enough neighbors to check whether
                   !the neighboring bc_interior_pt is suspicious or not,
                   !a temporary array is created
                   else
                      
                      !extract the grdpts_id around the bc_interior_pt
                      !that could not be checked
                      tmp_grdpts_id = get_grdpts_id_tmp_to_check_neighbors(
     $                     this,
     $                     bf_sublayer_updated,
     $                     [i,j],
     $                     match_table)

                      !then the neighboring points can be tested to see
                      !if the bc_interior_pt is a suspicious grid point
                      has_a_bc_pt_neighbor = verify_if_has_a_bc_pt_neighbor(
     $                     tmp_grdpts_id,
     $                     [bc_size+2,bc_size+2])

                      !if the bc_interior_pt has no bc_pt neighbor, then
                      !it is a suspicious bc_interior_pt that may be
                      !turned into an interior_pt
                      if(.not.has_a_bc_pt_neighbor) then

                         call update_bc_interior_pt_to_interior_pt(
     $                        this,
     $                        bf_sublayer_updated,
     $                        [i,j],
     $                        match_table,
     $                        tmp_grdpts_id=tmp_grdpts_id)

                      end if

                   end if

                end if

             end do
          end do
          
        end subroutine finalize_grdpts_for_suspicious_bc_interior_pt


        !update the bc_interior_pt to interior_pt if all the
        !grid points around the modified grid point exist
        subroutine update_bc_interior_pt_to_interior_pt(
     $     this,
     $     bf_sublayer_updated,
     $     local_coords,
     $     match_table,
     $     tmp_grdpts_id)

          implicit none

          class(nbf_interface_newgrdpt)                                , intent(in)    :: this
          type(bf_sublayer)                                            , intent(inout) :: bf_sublayer_updated
          integer(ikind), dimension(2)                                 , intent(in)    :: local_coords
          integer(ikind), dimension(2)                                 , intent(in)    :: match_table
          integer, dimension(2*(bc_size+1)+1,2*(bc_size+1)+1), optional, intent(in)    :: tmp_grdpts_id
          
          logical                                             :: all_grdpts_exist
          logical                                             :: ierror
          integer(ikind), dimension(2)                        :: gen_coords
          integer(ikind), dimension(2,2)                      :: gen_borders
          integer, dimension(2*(bc_size+1)+1,2*(bc_size+1)+1) :: tmp_grdpts_id1

          integer :: j


          if(present(tmp_grdpts_id)) then

             all_grdpts_exist = verify_if_all_grdpts_exist(
     $            bc_size+2,
     $            bc_size+2,
     $            tmp_grdpts_id)

          else

             !ask the buffer layer to check whether all the grid points
             !around the local grid point (local_coord) exist such that
             !it can be turned into an interior_pt
             all_grdpts_exist =
     $            bf_sublayer_updated%verify_if_all_bf_grdpts_exist(
     $            local_coords(1),local_coords(2),ierror)
             
             
             !if the buffer layer does not contain all the grid points
             !needed to determine whether all the grid points exist,
             !a temporary array with the grid points collected from
             !the interior, the neighbors and then the buffer layer
             !is initialized
             !the order is important (interior,neighbors,buffer layer)
             !as the status of the grid points in the buffer layer
             !should erase the other ones as the currently updated
             !buffer layer is the one with the last updated grid points
             !only grid points not appearing in the buffer layer should
             !be gathered
             if(ierror.neqv.BF_SUCCESS) then
             
                !1) extract the grdpts_id around the suspicious
                !   bc_interior_pt
                tmp_grdpts_id1 = get_grdpts_id_tmp_to_check_neighbors(
     $               this,
     $               bf_sublayer_updated,
     $               local_coords,
     $               match_table)
                
                !2) verify whether all grid points are present
                !   to turn the bc_interior_pt into interior_pt
                all_grdpts_exist = verify_if_all_grdpts_exist(
     $               bc_size+2,
     $               bc_size+2,
     $               tmp_grdpts_id1)
     $               
             end if

          end if


          !if all grid points exist to set as new interior_pt
          !the bc_interior_pt, the buffer layer is updated
          if(all_grdpts_exist) then

             call bf_sublayer_updated%set_new_interior_pt(
     $            local_coords(1),
     $            local_coords(2))

          !otherwise, this is a case which is not handled by
          !the program and further investigation is needed
          else
             
             print '(''nbf_interface_newgrdpt_class'')'
             print '(''update_bc_interior_pt_to_interior_pt'')'
             print '(''not all grdpts exist to turn the'')'
             print '(''bc_interior_pt into an interior_pt'')'
             print '()'
             
             print '(''this problem arises in the buffer layer:'')'
             print '(''localization: '',I2)', bf_sublayer_updated%get_localization()
             print '(''bf_align: '',2I2)', bf_sublayer_updated%get_alignment_tab()
             print '(''grid point coords: '',2I2)', local_coords
             print '()'

             print '(''neighboring grid points:'')'
             if(ierror.eqv.BF_SUCCESS) then

                gen_borders(1,1) = gen_coords(1)-(bc_size+1)
                gen_borders(1,2) = gen_coords(1)+(bc_size+1)
                gen_borders(2,1) = gen_coords(2)-(bc_size+1)
                gen_borders(2,2) = gen_coords(2)+(bc_size+1)

                call bf_sublayer_updated%get_grdpts_id_part(
     $               tmp_grdpts_id1,
     $               gen_borders)

             end if

             do j=1,2*(bc_size+1)+1
                print '(5I2)', tmp_grdpts_id1(:,2*bc_size+4-j)
             end do
             print '()'
             stop ''

          end if

        end subroutine update_bc_interior_pt_to_interior_pt


        !extract the grdpts_id
        function get_grdpts_id_tmp_to_check_neighbors(
     $     this,
     $     bf_sublayer_updated,
     $     bf_local_coords,
     $     bf_match_table)
     $     result(tmp_grdpts_id)

          implicit none
          
          class(nbf_interface_newgrdpt)                 , intent(in) :: this
          type(bf_sublayer)                             , intent(in) :: bf_sublayer_updated
          integer(ikind), dimension(2)                  , intent(in) :: bf_local_coords
          integer(ikind), dimension(2)                  , intent(in) :: bf_match_table
          integer(ikind), dimension(2*(bc_size+1)+1,2*(bc_size+1)+1) :: tmp_grdpts_id
          

          integer                        :: bf_localization
          integer(ikind), dimension(2)   :: gen_coords
          integer(ikind), dimension(2,2) :: gen_borders


          !get the cardinal coordinate of the buffer layer
          bf_localization = bf_sublayer_updated%get_localization()

          !get the general coordinates of the grid point
          gen_coords(1) = bf_local_coords(1) + bf_match_table(1)
          gen_coords(2) = bf_local_coords(2) + bf_match_table(2)

          !i_min,i_max,j_min,j_max
          gen_borders(1,1) = gen_coords(1)-(bc_size+1)
          gen_borders(1,2) = gen_coords(1)+(bc_size+1)
          gen_borders(2,1) = gen_coords(2)-(bc_size+1)
          gen_borders(2,2) = gen_coords(2)+(bc_size+1)


          !1) collect data from interior
          call get_grdpts_id_from_interior(
     $         tmp_grdpts_id,
     $         gen_borders)


          !2) collect data from potential neighbors
          !- gather data from neighbor of type 1
          !if(gen_borders(1,1).le.0) then
          if(bf_sublayer_updated%can_exchange_with_neighbor1()) then
             
             call this%get_grdpts_id_part(
     $            bf_localization,1,
     $            tmp_grdpts_id,
     $            gen_borders)

          end if

          !- gather data from neighbor of type 2
          !if(gen_borders(1,2).ge.(nx+1)) then
          if(bf_sublayer_updated%can_exchange_with_neighbor2()) then

             call this%get_grdpts_id_part(
     $            bf_localization,2,
     $            tmp_grdpts_id,
     $            gen_borders)

          end if


          !3) collect data from buffer layer
          call bf_sublayer_updated%get_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_borders)

        end function get_grdpts_id_tmp_to_check_neighbors


        function verify_if_has_a_bc_pt_neighbor(
     $     tmp_grdpts_id,
     $     local_coords)
     $     result(has_a_bc_pt_neighbor)

          implicit none

          integer       , dimension(:,:), intent(in) :: tmp_grdpts_id
          integer(ikind), dimension(2)  , intent(in) :: local_coords
          logical                                    :: has_a_bc_pt_neighbor

          integer(ikind) :: i,j

          has_a_bc_pt_neighbor = .false.

          do j=local_coords(2)-1,local_coords(2)+1
             do i=local_coords(1)-1, local_coords(1)+1
                if(tmp_grdpts_id(i,j).eq.bc_pt) then
                   has_a_bc_pt_neighbor = .true.
                   exit
                end if
             end do
          end do

        end function verify_if_has_a_bc_pt_neighbor


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> control whether the update the bc_interior_pt into
        !> an interior_pt lead to a crenel of bc_pt. If so, the
        !> crenel should be curb
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param bc_grdpt_local_coords
        !> coordinates of the bc_pt which is checked for bc_pt_crenel
        !
        !> @param match_table
        !> array used to convert coordinates expressed in the general
        !> reference frame of the interior domain into local coordinates
        !> of the buffer layer bf_sublayer_updated
        !---------------------------------------------------------------
        subroutine finalize_grdpts_for_bc_pt_crenel(
     $     this,
     $     bf_sublayer_updated,
     $     cpt_local_coords,
     $     match_table)

          implicit none

          class(nbf_interface_newgrdpt), intent(in)    :: this
          type(bf_sublayer)            , intent(inout) :: bf_sublayer_updated
          integer(ikind), dimension(2) , intent(in)    :: cpt_local_coords
          integer(ikind), dimension(2) , intent(in)    :: match_table

          
          integer(ikind) :: i,j


          !check whether the potential bc_pt around the cpt_local_coords
          !lead to a bc_pt_crenel
          j = cpt_local_coords(2)-bc_size
          do i=cpt_local_coords(1)-bc_size,cpt_local_coords(1)+bc_size

             if(bf_sublayer_updated%is_bc_pt(i,j)) then

                call control_bc_pt_crenel(
     $               this,
     $               bf_sublayer_updated,
     $               [i,j],
     $               match_table)

             end if
             
          end do

          do j=cpt_local_coords(2)-bc_size+1,cpt_local_coords(2)+bc_size-1

             i=cpt_local_coords(1)-bc_size
             if(bf_sublayer_updated%is_bc_pt(i,j)) then

                call control_bc_pt_crenel(
     $               this,
     $               bf_sublayer_updated,
     $               [i,j],
     $               match_table)

             end if

             i=cpt_local_coords(1)+bc_size
             if(bf_sublayer_updated%is_bc_pt(i,j)) then

                call control_bc_pt_crenel(
     $               this,
     $               bf_sublayer_updated,
     $               [i,j],
     $               match_table)

             end if

          end do

          j = cpt_local_coords(2)+bc_size
          do i=cpt_local_coords(1)-bc_size,cpt_local_coords(1)+bc_size

             if(bf_sublayer_updated%is_bc_pt(i,j)) then

                call control_bc_pt_crenel(
     $               this,
     $               bf_sublayer_updated,
     $               [i,j],
     $               match_table)

             end if
             
          end do

        end subroutine finalize_grdpts_for_bc_pt_crenel


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> control whether the update the bc_interior_pt into
        !> an interior_pt lead to a crenel of bc_pt. If so, the
        !> crenel should be curb
        !
        !> @date
        !> 21_01_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param bc_grdpt_local_coords
        !> coordinates of the bc_pt which is checked for bc_pt_crenel
        !
        !> @param match_table
        !> array used to convert coordinates expressed in the general
        !> reference frame of the interior domain into local coordinates
        !> of the buffer layer bf_sublayer_updated
        !---------------------------------------------------------------
        subroutine control_bc_pt_crenel(
     $     this,
     $     bf_sublayer_updated,
     $     bc_grdpt_local_coords,
     $     match_table)

          implicit none

          class(nbf_interface_newgrdpt), intent(in)    :: this
          type(bf_sublayer)            , intent(inout) :: bf_sublayer_updated
          integer(ikind), dimension(2) , intent(in)    :: bc_grdpt_local_coords
          integer(ikind), dimension(2) , intent(in)    :: match_table


          integer                                            :: bf_localization
          integer(ikind), dimension(2)                       :: bf_sizes
          integer       , dimension(2*bc_size+1,2*bc_size+1) :: tmp_grdpts_id
          logical                                            :: tmp_grdpts_id_needed
          integer(ikind), dimension(2,2)                     :: gen_borders
          logical                                            :: bc_pt_crenel_exists

          bf_localization = bf_sublayer_updated%get_localization()
          bf_sizes        = bf_sublayer_updated%get_sizes()


          !1) determine whether the control of the bc_pt_crenel
          !   requires grdpts_id from neighboring buffer layers
          tmp_grdpts_id_needed = is_temp_array_needed_for_bc_crenel(
     $         bf_localization,
     $         bf_sizes,
     $         bc_grdpt_local_coords)


          !2) if a temporary array is needed, the grdpts_id are
          !   extracted from the neighboring buffer layers, then
          !   the bc_pt_crenel is checked. If there is a
          !   bc_pt_crenel, it is curb and the temporary array
          !   is substituted in the buffer layers matching the
          !   gen_borders
          if(tmp_grdpts_id_needed) then

             !2.1) extract the grpts_id

             !2.1.1) determine the general coordinates identifying
             !       the borders of the grid points ID extracted
             gen_borders(1,1) = match_table(1) + bc_grdpt_local_coords(1) - bc_size
             gen_borders(1,2) = match_table(1) + bc_grdpt_local_coords(1) + bc_size
             gen_borders(2,1) = match_table(2) + bc_grdpt_local_coords(2) - bc_size
             gen_borders(2,2) = match_table(2) + bc_grdpt_local_coords(2) + bc_size


             !2.1.2) collect data from interior
             call get_grdpts_id_from_interior(
     $            tmp_grdpts_id,
     $            gen_borders)
             

             !2.1.3) collect data from potential neighbors
             !- gather data from neighbor of type 1
             !if(gen_borders(1,1).le.0) then
             if(bf_sublayer_updated%can_exchange_with_neighbor1()) then
                
                call this%get_grdpts_id_part(
     $               bf_localization,1,
     $               tmp_grdpts_id,
     $               gen_borders)
                
             end if
   
             !- gather data from neighbor of type 2
             !if(gen_borders(1,2).ge.(nx+1)) then
             if(bf_sublayer_updated%can_exchange_with_neighbor2()) then
   
                call this%get_grdpts_id_part(
     $               bf_localization,2,
     $               tmp_grdpts_id,
     $               gen_borders)
   
             end if
   
   
             !2.1.4) collect data from buffer layer
             call bf_sublayer_updated%get_grdpts_id_part(
     $            tmp_grdpts_id,
     $            gen_borders)



             !2.2) control and curb the bc_pt_crenel if any
             bc_pt_crenel_exists = detect_and_curb_bc_crenels(
     $            [bc_size+1,bc_size+1],
     $            [2*bc_size+1,2*bc_size+1],
     $            tmp_grdpts_id)


             !2.3) if the bc_pt_crenel exists, substitute the 
             !     grdpts_id in the buffer layers sharing
             !     gridpoints with this temporary array
             if(bc_pt_crenel_exists) then

                !- set data in neighbors of type 1
                if(bf_sublayer_updated%can_exchange_with_neighbor1()) then
                
                   call this%set_grdpts_id_part(
     $                  bf_localization,1,
     $                  tmp_grdpts_id,
     $                  gen_borders)
                   
                end if
   
                !- set data in neighbors of type 2
                if(bf_sublayer_updated%can_exchange_with_neighbor2()) then
                   
                   call this%set_grdpts_id_part(
     $                  bf_localization,2,
     $                  tmp_grdpts_id,
     $                  gen_borders)
                   
                end if   
   
                !- set data in current buffer layer
                call bf_sublayer_updated%set_grdpts_id_part(
     $               tmp_grdpts_id,
     $               gen_borders)

             end if

          !3) if no temporary array is needed, the grdpts_id are
          !   directly checked in the current buffer layer
          else

             bc_pt_crenel_exists =
     $            bf_sublayer_updated%detect_and_curb_bc_pt_crenels(
     $            bc_grdpt_local_coords)

          end if

        end subroutine control_bc_pt_crenel

      end module nbf_interface_newgrdpt_class
