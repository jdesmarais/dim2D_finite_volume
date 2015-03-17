      !> @file
      !> mainlayer_interface_newgrdpt augmented with procedures
      !> updating the configuration of the grdpts_id
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_newgrdpt augmented with procedures
      !> updating the configuration of the grdpts_id
      !
      !> @date
      ! 17_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_grdpts_id_update_class

        use bf_sublayer_class, only :
     $       bf_sublayer

        use mainlayer_interface_newgrdpt_class, only :
     $       mainlayer_interface_newgrdpt

        use parameters_bf_layer, only :
     $       interior_pt

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       debug_geometry_update

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        private
        public :: mainlayer_interface_grdpts_id_update


        !>@class mainlayer_interface_grdpts_id_update
        !> mainlayer_interface_newgrdpt augmented with procedures
        !> updating the configuration of the grdpts_id
        !
        !> @param update_grdpts_id_in_bf_layer
        !> turn the bc_interior_pt grid-points identified by their
        !> general coordinates into interior_pt and update the
        !> configuration of the grid-points around it to make sure
        !> they are computed during the time integration
        !--------------------------------------------------------------
        type, extends(mainlayer_interface_newgrdpt) :: mainlayer_interface_grdpts_id_update

          contains

          procedure, pass :: update_grdpts_id_in_bf_layer

          ! procedures for the update the grid points ID
          procedure, pass :: update_grdpts_id_around_new_interior_pt
          procedure, pass :: update_grdpts_id_around_new_interior_pt_local

          ! procedures for the finalization of crenels
          

        end type mainlayer_interface_grdpts_id_update


        contains


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
        subroutine update_grdpts_id_in_bf_layer(
     $       this,
     $       p_model,
     $       t,
     $       dt,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       bf_sublayer_ptr,
     $       selected_grdpts)

          implicit none

          class(mainlayer_interface_grdpts_id_update), intent(inout) :: this
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          real(rkind)   , dimension(nx)              , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes1
          type(bf_sublayer)                          , intent(inout) :: bf_sublayer_ptr
          integer(ikind), dimension(:,:)             , intent(in)    :: selected_grdpts


          integer(ikind), dimension(2) :: match_table
          integer(ikind)               :: k
          integer(ikind)               :: gen_i_prev
          integer(ikind)               :: gen_j_prev
          integer(ikind)               :: gen_i_center
          integer(ikind)               :: gen_j_center


          match_table = bf_sublayer_ptr%get_general_to_local_coord_tab()
          
          
          !we have a list of gridpoints that should be turned
          !from bc_interior_pt to interior_pt. For a point to be
          !considered an interior point, we need to make sure
          !that the grid points needed to compute its time
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
          gen_i_center = nx/2
          gen_j_center = ny/2

          do k=1, size(selected_grdpts,2)

             !update the position of the gridpoint
             !previously tested
             gen_i_prev   = gen_i_center
             gen_j_prev   = gen_j_center
             gen_i_center = selected_grdpts(1,k)
             gen_j_center = selected_grdpts(2,k)

             !update the configuration of the grid-points
             !around the new interior point
             call update_grdpts_id_around_new_interior_pt(
     $            this,
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            bf_sublayer_ptr,
     $            gen_i_prev,
     $            gen_j_prev,
     $            gen_i_center,
     $            gen_j_center,
     $            match_table)

             

c$$$             !check whether the neighboring bc_interior_pt
c$$$             !should be updated to interior_pt
c$$$             call finalize_grdpts_id_for_suspicious_bc_interior_pt(
c$$$     $            this,
c$$$     $            bf_sublayer_ptr,
c$$$     $            [gen_i_center,gen_j_center],
c$$$     $            match_table)
c$$$
c$$$             !check whether the update the grdpts_id lead to
c$$$             !a crenel of bc_pt
c$$$             call finalize_grdpts_id_for_bc_pt_crenel(
c$$$     $            this,
c$$$     $            bf_sublayer_ptr,
c$$$     $            [gen_i_center,gen_j_center],
c$$$     $            match_table)

          end do

        end subroutine update_grdpts_id_in_bf_layer


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> update the grid points around the new interior_pt
        !
        !> @date
        !> 17_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
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
        !> @param bf_sublayer_ptr
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !        
        !> @param gen_i_prev
        !> x-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_j_prev
        !> y-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_i_center
        !> x-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_j_center
        !> y-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param match_table
        !> correspondance table between the local and
        !> the general coordinates
        !---------------------------------------------------------------
        subroutine update_grdpts_id_around_new_interior_pt(
     $     this,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_ptr,
     $     gen_i_prev,
     $     gen_j_prev,
     $     gen_i_center,
     $     gen_j_center,
     $     match_table)

          implicit none

          class(mainlayer_interface_grdpts_id_update), intent(inout) :: this
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          real(rkind)   , dimension(nx)              , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes1
          type(bf_sublayer)                          , intent(inout) :: bf_sublayer_ptr
          integer(ikind)                             , intent(in)    :: gen_i_prev
          integer(ikind)                             , intent(in)    :: gen_j_prev
          integer(ikind)                             , intent(in)    :: gen_i_center
          integer(ikind)                             , intent(in)    :: gen_j_center
          integer(ikind), dimension(2)               , intent(in)    :: match_table

          integer(ikind) :: min_j, max_j
          integer(ikind) :: i,j   
          integer(ikind) :: loc_i_prev
          integer(ikind) :: loc_j_prev
          integer(ikind) :: loc_i_center
          integer(ikind) :: loc_j_center                


          min_j = min(gen_j_center-gen_j_prev,0)
          max_j = max(gen_j_center-gen_j_prev,0)


          !check whether the grid points around the new interior_pt
          !exist: if they do not exist, they are computed
          !in any case, the status of the neighboring points may be
          !updated: from bc_pt -> bc_interior_pt or no_pt -> bc_pt
          do j=gen_j_center-bc_size, gen_j_prev - bc_size + min(gen_j_center-gen_j_prev+2*bc_size,-1)
             do i=gen_i_center-bc_size,gen_i_center+bc_size
                
                call this%update_grdpts_id_around_new_interior_pt_local(
     $               p_model,t,dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               i,j,
     $               gen_i_center,
     $               gen_j_center,
     $               match_table)

             end do
          end do

          do j=gen_j_center-bc_size-min_j, gen_j_center+bc_size-max_j
             do i=gen_i_center-bc_size, gen_i_prev-bc_size+min(gen_i_center-gen_i_prev+2*bc_size,-1)

                call this%update_grdpts_id_around_new_interior_pt_local(
     $               p_model,t,dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               i,j,
     $               gen_i_center,
     $               gen_j_center,
     $               match_table)

             end do
          end do

          do j=gen_j_center-bc_size-min_j, gen_j_center+bc_size-max_j
             do i=gen_i_prev+bc_size+max(gen_i_center-gen_i_prev-2*bc_size,1),gen_i_center+bc_size

                call this%update_grdpts_id_around_new_interior_pt_local(
     $               p_model,t,dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               i,j,
     $               gen_i_center,
     $               gen_j_center,
     $               match_table)

             end do
          end do

          do j=gen_j_prev+bc_size+max(gen_j_center-gen_j_prev-2*bc_size,1), gen_j_center+bc_size
             do i=gen_i_center-bc_size,gen_i_center+bc_size

                call this%update_grdpts_id_around_new_interior_pt_local(
     $               p_model,t,dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               i,j,
     $               gen_i_center,
     $               gen_j_center,
     $               match_table)

             end do
          end do


          !finalize the update of the neighboring gridpoints
          !grid points that do not need to be computed but whose
          !status is modified from bc_pt to bc_interior_pt
          loc_i_prev   = gen_i_prev   - match_table(1)
          loc_j_prev   = gen_j_prev   - match_table(2)
          loc_i_center = gen_i_center - match_table(1)
          loc_j_center = gen_j_center - match_table(2)

          call bf_sublayer_ptr%finalize_update_grdpts_id_around_new_interior_pt(
     $         loc_i_prev,
     $         loc_j_prev,
     $         loc_i_center,
     $         loc_j_center)


          !update the status of the central gridpoint
          !to interior_pt
          call bf_sublayer_ptr%set_grdpts_id_pt(
     $         loc_i_center,
     $         loc_j_center,
     $         interior_pt)

        end subroutine update_grdpts_id_around_new_interior_pt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> update the grid points around the new interior_pt
        !
        !> @date
        !> 17_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> mainlayer_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
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
        !> @param bf_sublayer_ptr
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !        
        !> @param gen_i
        !> x-index of the bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_j
        !> y-index of the bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_i_center
        !> x-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param gen_j_center
        !> y-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !
        !> @param match_table
        !> correspondance table between the local and
        !> the general coordinates
        !---------------------------------------------------------------
        subroutine update_grdpts_id_around_new_interior_pt_local(
     $     this,
     $     p_model,
     $     t,
     $     dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_ptr,
     $     gen_i,
     $     gen_j,
     $     gen_i_center,
     $     gen_j_center,
     $     match_table)

          implicit none

          class(mainlayer_interface_grdpts_id_update), intent(inout) :: this
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          real(rkind)   , dimension(nx)              , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes1
          type(bf_sublayer)                          , intent(inout) :: bf_sublayer_ptr
          integer(ikind)                             , intent(in)    :: gen_i
          integer(ikind)                             , intent(in)    :: gen_j
          integer(ikind)                             , intent(in)    :: gen_i_center
          integer(ikind)                             , intent(in)    :: gen_j_center
          integer(ikind), dimension(2)               , intent(in)    :: match_table


          logical :: compute_grdpt

          integer(ikind) :: loc_i
          integer(ikind) :: loc_j
          integer(ikind) :: loc_i_center
          integer(ikind) :: loc_j_center

          logical :: ierror
          

          loc_i        = gen_i        - match_table(1)
          loc_j        = gen_j        - match_table(2)
          loc_i_center = gen_i_center - match_table(1)
          loc_j_center = gen_j_center - match_table(2)

          
          compute_grdpt = bf_sublayer_ptr%update_grdpt_id_next_to_new_interior_pt(
     $         loc_i,
     $         loc_j,
     $         loc_i_center,
     $         loc_j_center)

          !in some tests, only the update of the geometry
          !of the computational is tested
          if(.not.debug_geometry_update) then

             if(compute_grdpt) then
                ierror = this%compute_newgrdpt(
     $               p_model,
     $               t,
     $               dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               [gen_i,gen_j])
             end if

          end if

        end subroutine update_grdpts_id_around_new_interior_pt_local        

      end module mainlayer_interface_grdpts_id_update_class
