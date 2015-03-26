      program test_bf_newgrdpt_procedure

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use bf_newgrdpt_procedure_module, only :
     $       get_config_id,
     $       get_newgrdpt_procedure

        use check_data_module, only :
     $       is_real_vector_validated

        use parameters_bf_layer, only :
     $       bc_pt,
     $       no_pt,
     $       BF_SUCCESS

        use parameters_constant, only :
     $       vector_x,
     $       vector_y,
     $       newgrdpt_test

        use parameters_input, only :
     $       ne,
     $       ic_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        test_loc = test_get_config_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_config_id: '',L1)', test_loc
        print '()'

        test_loc = test_get_newgrdpt_procedure(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_newgrdpt_procedure: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains

        function test_get_config_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                 :: config_id_test
          integer, dimension(3,3) :: grdpts_id
          integer                 :: config_id
          logical                 :: test_loc

          test_validated = .true.


          do config_id_test=0,255

             !input
             grdpts_id = get_grdpts_id_config(config_id_test)

             !output 
             config_id = get_config_id(2,2,grdpts_id)
             
             !validation
             test_loc = config_id.eq.config_id_test
             test_validated = test_validated.and.test_loc

             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test: '',I3, '' -> '',I3)', config_id, config_id_test
             end if

          end do          

        end function test_get_config_id

      
        function test_get_newgrdpt_procedure(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc
          integer :: k


          test_validated = .true.

          do k=0,255

             test_loc = test_get_newgrdpt_procedure_config(k,detailled)
             test_validated = test_validated.and.test_loc

          end do


        end function test_get_newgrdpt_procedure


        function test_get_newgrdpt_procedure_config(
     $     config_id,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: config_id
          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)    :: t
          real(rkind)    :: dt
          integer(ikind) :: i1
          integer(ikind) :: j1

          integer(ikind), dimension(2,2) :: bf_align0
          real(rkind), dimension(7)      :: bf_x_map0
          real(rkind), dimension(7)      :: bf_y_map0
          real(rkind), dimension(7,7,ne) :: bf_nodes0
          integer(ikind), dimension(2,2) :: bf_align1
          real(rkind), dimension(7)      :: bf_x_map1
          real(rkind), dimension(7)      :: bf_y_map1
          real(rkind), dimension(7,7,ne) :: bf_nodes1

          integer    , dimension(3,3)    :: grdpts_id

          type(pmodel_eq)                :: p_model
          
          integer                        :: ierror_test
          integer                        :: ierror
          real(rkind), dimension(ne)     :: newgrdpt_data
          real(rkind), dimension(ne)     :: newgrdpt


          test_validated = .true.

          if(ic_choice.ne.newgrdpt_test) then
             stop 'the test requires ic_choice=newgrdpt_test'
          end if

          
          !0.0) initialization
          t  = 0.0d0
          dt = 0.25d0
          i1 = 4
          j1 = 4

          bf_align0(1,1) = 0
          bf_align0(2,1) = 0
          bf_align1(1,1) = 0
          bf_align1(2,1) = 0

          call p_model%initial_conditions%ini_far_field()
          

          !0.1) initialize the nodes and maps
          call initialize_nodes_and_maps(
     $         bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_x_map1, bf_y_map1, bf_nodes1)

          !0.2) prepare the grdpts_id configuration
          grdpts_id = get_grdpts_id_config(config_id)

          !1) compute the new grid point
          call compute_newgrdpt(
     $         newgrdpt_data,
     $         ierror_test,
     $         grdpts_id,
     $         p_model,t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1)

          !2) determine the new gridpoint if the field is reflected
          !   in the x-direction
          call reflect_x(
     $         p_model,
     $         grdpts_id,
     $         bf_x_map0, bf_nodes0, 
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data)

          call compute_newgrdpt(
     $         newgrdpt,
     $         ierror,
     $         grdpts_id,
     $         p_model,t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1)

          if(ierror.eq.BF_SUCCESS) then
             test_loc = is_real_vector_validated(newgrdpt,newgrdpt_data,detailled)
          else
             test_loc = ierror.eq.ierror_test
          end if
          test_validated = test_validated.and.test_loc

          !validation
          if((.not.test_loc).and.detailled) then
             print '(''test_reflect_x('',I3,'',1): '',L1)', config_id, test_loc
          end if


          !5) determine the new gridpoint if the field is reflected
          !   in the y-direction
          call reflect_y(
     $         p_model,
     $         grdpts_id,
     $         bf_y_map0, bf_nodes0, 
     $         bf_y_map1, bf_nodes1,
     $         j1,
     $         newgrdpt_data)

          call compute_newgrdpt(
     $         newgrdpt,
     $         ierror,
     $         grdpts_id,
     $         p_model,t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1)

          if(ierror.eq.BF_SUCCESS) then
             test_loc = is_real_vector_validated(newgrdpt,newgrdpt_data,detailled)
          else
             test_loc = ierror.eq.ierror_test
          end if
          test_validated = test_validated.and.test_loc

          !validation
          if((.not.test_loc).and.detailled) then
             print '(''test_reflect_y('',I3,'',2): '',L1)', config_id, test_loc
          end if


          !6) determine the new gridpoint if the field is reflected
          !   in the x-direction
          call reflect_x(
     $         p_model,
     $         grdpts_id,
     $         bf_x_map0, bf_nodes0, 
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data)

          call compute_newgrdpt(
     $         newgrdpt,
     $         ierror,
     $         grdpts_id,
     $         p_model,t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1)

          if(ierror.eq.BF_SUCCESS) then
             test_loc = is_real_vector_validated(newgrdpt,newgrdpt_data,detailled)
          else
             test_loc = ierror.eq.ierror_test
          end if
          test_validated = test_validated.and.test_loc
          
          !validation
          if((.not.test_loc).and.detailled) then
             print '(''test_reflect_x('',I3,'',2): '',L1)', config_id, test_loc
          end if

        end function test_get_newgrdpt_procedure_config


        subroutine compute_newgrdpt(
     $     newgrdpt,
     $     ierror,
     $     grdpts_id,
     $     p_model,t,dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1)

          implicit none

          real(rkind)    , dimension(ne) , intent(out) :: newgrdpt
          integer                        , intent(out) :: ierror
          integer        , dimension(3,3), intent(in)  :: grdpts_id
          type(pmodel_eq)                , intent(in)  :: p_model
          real(rkind)                    , intent(in)  :: t
          real(rkind)                    , intent(in)  :: dt
          integer(ikind), dimension(2,2) , intent(in)  :: bf_align0
          real(rkind), dimension(7)      , intent(in)  :: bf_x_map0
          real(rkind), dimension(7)      , intent(in)  :: bf_y_map0
          real(rkind), dimension(7,7,ne) , intent(in)  :: bf_nodes0
          integer(ikind), dimension(2,2) , intent(in)  :: bf_align1
          real(rkind), dimension(7)      , intent(in)  :: bf_x_map1
          real(rkind), dimension(7)      , intent(in)  :: bf_y_map1
          real(rkind), dimension(7,7,ne) , intent(in)  :: bf_nodes1
          integer(ikind)                 , intent(in)  :: i1
          integer(ikind)                 , intent(in)  :: j1

          
          integer                 :: nb_procedures
          integer, dimension(4)   :: procedure_type
          integer, dimension(4)   :: gradient_type

          type(bf_newgrdpt)            :: bf_newgrdpt_used

          !1) ask for the procedures needed to compute
          !   the new grid point
          ierror = get_newgrdpt_procedure(
     $         2,2,grdpts_id,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type)

          !2) compute the new grid point
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model,t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         nb_procedures,
     $         procedure_type, gradient_type)

        end subroutine compute_newgrdpt


        subroutine initialize_nodes_and_maps(
     $     bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_x_map1, bf_y_map1, bf_nodes1)

          implicit none

          real(rkind), dimension(7)     , intent(out) :: bf_x_map0
          real(rkind), dimension(7)     , intent(out) :: bf_y_map0
          real(rkind), dimension(7,7,ne), intent(out) :: bf_nodes0
          real(rkind), dimension(7)     , intent(out) :: bf_x_map1
          real(rkind), dimension(7)     , intent(out) :: bf_y_map1
          real(rkind), dimension(7,7,ne), intent(out) :: bf_nodes1
          

          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0, 5.5d0, 6.5d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.5d0]
          bf_nodes0 = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.29d0, 1.42d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.26d0, 1.48d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.29d0, 1.24d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0, 1.37d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.45d0, 1.38d0,
     $         1.41d0, 1.26d0, 1.25d0, 1.32d0, 1.38d0, 1.34d0, 1.42d0,
     $         1.34d0, 1.22d0, 1.34d0, 1.22d0, 1.28d0, 1.44d0, 1.22d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.121d0,0.118d0,
     $         0.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.142d0,0.135d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.153d0,0.162d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.119d0,0.112d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.158d0,0.153d0,
     $         0.158d0, 0.137d0, 0.112d0, 0.169d0, 0.146d0, 0.127d0,0.139d0,
     $         0.138d0, 0.148d0, 0.186d0, 0.123d0, 0.176d0, 0.128d0,0.112d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0, 0.061d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.085d0, 0.082d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.015d0, 0.017d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.042d0, 0.045d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.023d0, 0.021d0,
     $          0.015d0, 0.04d0 , 0.012d0,  0.06d0, 0.075d0, 0.034d0, 0.073d0,
     $          0.084d0, 0.011d0,  0.09d0, 0.015d0, 0.052d0, 0.032d0, 0.043d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0,4.834d0,4.824d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0,4.875d0,4.865d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0,4.826d0,4.626d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0,4.952d0,4.572d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0, 4.95d0,
     $         4.71d0, 4.653d0, 4.962d0,  4.92d0, 4.712d0, 4.87d0, 4.92d0,
     $         4.86d0, 4.865d0, 4.840d0, 4.816d0, 4.813d0,4.876d0,4.826d0
     $         /),
     $         (/7,7,ne/))

          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0,4.5d0,5.5d0,6.5d0]
          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0, 1.0d0, 1.25d0, 1.50d0]
          bf_nodes1 = reshape((/
     $         1.50d0, 1.455d0, 1.48d0, 1.29d0, 1.45d0 ,1.35d0, 1.545d0,
     $         1.20d0, 1.350d0, 1.25d0, 1.25d0, 1.36d0 ,1.22d0, 1.035d0,
     $         1.49d0, 1.250d0, 1.40d0, 1.54d0, 1.256d0,1.41d0, 1.205d0,
     $         1.52d0, 1.28d0,  1.45d0, 1.52d0, 1.29d0 ,1.56d0, 1.82d0, 
     $         1.59d0, 1.26d0,  1.26d0, 1.37d0, 1.48d0 ,1.95d0, 1.22d0, 
     $         1.412d0, 1.261d0, 1.253d0, 1.322d0, 1.384d0, 1.345d0, 1.421d0,
     $         1.343d0, 1.223d0, 1.343d0, 1.221d0, 1.285d0, 1.442d0, 1.227d0,
     $         
     $         0.128d0, 0.150d0, 0.135d0,  0.14d0, 0.18d0  , 0.14d0 , 0.182d0,
     $         0.148d0, 0.150d0, 0.122d0, 0.136d0, 0.1742d0, 0.163d0, 0.184d0,
     $         0.142d0, 1.152d0, 0.136d0,  0.12d0, 0.1592d0, 0.21d0 , 0.124d0,
     $         0.136d0, 0.185d0, 0.196d0, 0.192d0, 0.1263d0, 0.129d0, 0.316d0,
     $         0.158d0,0.1265d0,  0.13d0, 0.182d0, 0.1426d0, 0.128d0, 0.581d0,
     $         0.148d0, 0.150d0, 0.122d0, 0.136d0, 0.1742d0, 0.163d0, 0.418d0,
     $         0.136d0, 0.185d0, 0.196d0, 0.192d0, 0.1263d0, 0.129d0, 0.163d0,
     $         
     $         0.006d0,	0.0600d0, 0.020d0, 0.052d0, 0.562d0, 0.002d0,   0.0060d0,
     $         0.000d0,	0.0028d0, 0.035d0, 0.4825d0,0.365d0, 0.000d0,	0.0082d0,
     $         0.020d0,	0.0030d0, 0.040d0, 0.085d0, 0.0269d0,0.014d0,	0.0300d0,
     $         0.0151d0, 0.042d0 ,0.0123d0, 0.063d0, 0.0715d0, 0.0434d0, 0.0732d0,
     $         0.0842d0, 0.0113d0,0.029d0, 0.0152d0, 0.0532d0, 0.0232d0, 0.0431d0,
     $         0.0115d0, 0.024d0 ,0.0132d0, 0.036d0, 0.0571d0, 0.0434d0, 0.0723d0,
     $         0.0428d0, 0.0131d0,0.029d0, 0.0125d0, 0.0253d0, 0.0232d0, 0.0341d0,
     $         
     $         4.876d0, 4.825d0, 4.862d0, 4.821d0, 4.932d0, 4.826d0, 4.812d0,
     $         4.890d0, 4.871d0, 4.892d0, 4.875d0, 4.926d0, 4.829d0, 4.857d0,
     $         4.865d0, 4.757d0, 4.895d0, 4.856d0, 4.818d0, 4.859d0, 4.865d0,
     $         4.825d0, 4.785d0, 4.825d0, 4.826d0, 4.865d0, 4.852d0, 4.862d0,
     $         4.852d0, 4.826d0, 4.815d0, 4.825d0, 4.291d0, 4.851d0, 4.952d0,
     $         4.715d0, 4.635d0, 4.926d0, 4.912d0, 4.721d0, 4.87d0,  4.982d0,
     $         4.846d0, 4.856d0, 4.804d0, 4.861d0, 4.831d0, 4.867d0, 4.862d0
     $         /),
     $         (/7,7,ne/))

        end subroutine initialize_nodes_and_maps


        subroutine reflect_x(
     $     p_model,
     $     grdpts_id,
     $     bf_x_map0, bf_nodes0, 
     $     bf_x_map1, bf_nodes1,
     $     i1,
     $     newgrdpt_data)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          integer    , dimension(3,3)   , intent(inout) :: grdpts_id
          real(rkind), dimension(7)     , intent(inout) :: bf_x_map0
          real(rkind), dimension(7,7,ne), intent(inout) :: bf_nodes0
          real(rkind), dimension(7)     , intent(inout) :: bf_x_map1
          real(rkind), dimension(7,7,ne), intent(inout) :: bf_nodes1
          integer(ikind)                , intent(inout) :: i1
          real(rkind), dimension(:)     , intent(inout) :: newgrdpt_data

          integer :: i
          integer :: i_r
          integer :: j
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          integer    , dimension(3,3)    :: grdpts_id_r
          real(rkind), dimension(7)      :: bf_x_map0_r
          real(rkind), dimension(7,7,ne) :: bf_nodes0_r
          real(rkind), dimension(7)      :: bf_x_map1_r
          real(rkind), dimension(7,7,ne) :: bf_nodes1_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_x) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the x_maps
          do i=1, size(bf_x_map0,1)

             i_r = size(bf_x_map0)-(i-1)

             bf_x_map0_r(i) = -bf_x_map0(i_r)
             bf_x_map1_r(i) = -bf_x_map1(i_r)

          end do
          bf_x_map0 = bf_x_map0_r
          bf_x_map1 = bf_x_map1_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_x) then

                do j=1, size(bf_nodes0,2)
                   do i=1, size(bf_nodes0,1)
                      i_r = size(bf_nodes0,1)-(i-1)
                      bf_nodes0_r(i,j,k) = - bf_nodes0(i_r,j,k)
                      bf_nodes1_r(i,j,k) = - bf_nodes1(i_r,j,k)
                   end do
                end do

             else

                do j=1, size(bf_nodes0,2)
                   do i=1, size(bf_nodes0,1)
                      i_r = size(bf_nodes0,1)-(i-1)
                      bf_nodes0_r(i,j,k) = bf_nodes0(i_r,j,k)
                      bf_nodes1_r(i,j,k) = bf_nodes1(i_r,j,k)
                   end do
                end do

             end if             

          end do
          bf_nodes0 = bf_nodes0_r
          bf_nodes1 = bf_nodes1_r


          !reflect the i1,j1
          i1 = size(bf_x_map0,1) - (i1-1)


          !reflect the newgrdpt_data_test
          do k=1, ne
             if(var_type(k).eq.vector_x) then
                newgrdpt_data(k) = - newgrdpt_data(k)
             end if
          end do


          !reflect the grdpts_id
          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)
                i_r = size(grdpts_id,1)-i+1
                grdpts_id_r(i,j) = grdpts_id(i_r,j)
             end do
          end do
          grdpts_id = grdpts_id_r

        end subroutine reflect_x


        subroutine reflect_y(
     $     p_model,
     $     grdpts_id,
     $     bf_y_map0, bf_nodes0, 
     $     bf_y_map1, bf_nodes1,
     $     j1,
     $     newgrdpt_data)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          integer    , dimension(3,3)   , intent(inout) :: grdpts_id
          real(rkind), dimension(7)     , intent(inout) :: bf_y_map0
          real(rkind), dimension(7,7,ne), intent(inout) :: bf_nodes0
          real(rkind), dimension(7)     , intent(inout) :: bf_y_map1
          real(rkind), dimension(7,7,ne), intent(inout) :: bf_nodes1
          integer(ikind)                , intent(inout) :: j1
          real(rkind), dimension(:)     , intent(inout) :: newgrdpt_data

          integer :: i
          integer :: j
          integer :: j_r
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          integer    , dimension(3,3)    :: grdpts_id_r
          real(rkind), dimension(7)      :: bf_y_map0_r
          real(rkind), dimension(7,7,ne) :: bf_nodes0_r
          real(rkind), dimension(7)      :: bf_y_map1_r
          real(rkind), dimension(7,7,ne) :: bf_nodes1_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_y) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the y_maps
          do j=1, size(bf_y_map0,1)

             j_r = size(bf_y_map0,1)-(j-1)

             bf_y_map0_r(j) = -bf_y_map0(j_r)
             bf_y_map1_r(j) = -bf_y_map1(j_r)

          end do
          bf_y_map0 = bf_y_map0_r
          bf_y_map1 = bf_y_map1_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_y) then

                do j=1, size(bf_nodes0,2)

                   j_r = size(bf_nodes0,2)-(j-1)

                   do i=1, size(bf_nodes0,1)
                      bf_nodes0_r(i,j,k) = - bf_nodes0(i,j_r,k)
                      bf_nodes1_r(i,j,k) = - bf_nodes1(i,j_r,k)
                   end do

                end do

             else

                do j=1, size(bf_nodes0,2)

                   j_r = size(bf_nodes0,2)-(j-1)

                   do i=1, size(bf_nodes0,1)
                      bf_nodes0_r(i,j,k) = bf_nodes0(i,j_r,k)
                      bf_nodes1_r(i,j,k) = bf_nodes1(i,j_r,k)
                   end do

                end do

             end if             

          end do
          bf_nodes0 = bf_nodes0_r
          bf_nodes1 = bf_nodes1_r


          !reflect the j1
          j1 = size(bf_y_map0,1) - (j1-1)


          !reflect the newgrdpt_data_test
          do k=1, ne
             if(var_type(k).eq.vector_y) then
                newgrdpt_data(k) = - newgrdpt_data(k)
             end if
          end do


          !reflect the grdpts_id
          do j=1, size(grdpts_id,2)
             j_r = size(grdpts_id,2)-j+1
             do i=1, size(grdpts_id,1)
                grdpts_id_r(i,j) = grdpts_id(i,j_r)
             end do
          end do
          grdpts_id = grdpts_id_r

        end subroutine reflect_y


        function get_grdpts_id_config(config_id)
     $     result(grdpts_id)

          implicit none

          integer                , intent(in) :: config_id
          integer, dimension(3,3)             :: grdpts_id


          integer, dimension(8) :: bin
          integer, dimension(8) :: bin_r
          integer               :: i
          integer               :: j
          

          bin = get_bin_decomposition(config_id)

          do i=1, size(bin,1)
             bin_r(i) = bin(size(bin,1)-i+1)
          end do

          grdpts_id(1:3,1) = bin_r(1:3)
          grdpts_id(1:3,2) = [bin_r(4),0,bin_r(5)]
          grdpts_id(1:3,3) = bin_r(6:8)

          do j=1,3
             do i=1,3
                select case(grdpts_id(i,j))
                  case(0)
                     grdpts_id(i,j) = no_pt
                  case(1)
                     grdpts_id(i,j) = bc_pt
                  case default
                     print '(''test_bf_newgrdpt_procedure'')'
                     print '(''get_grdpts_id_grdpts_id'')'
                     print '(''grdpts_id('',2I2,''): '',I3)', i,j,grdpts_id(i,j)
                     print '(''not recognized'')'
                     stop ''
                end select
             end do
          end do
          
        end function get_grdpts_id_config


        function get_bin_decomposition(config_id)
     $     result(bin)

          implicit none

          integer              , intent(in) :: config_id
          integer, dimension(8)             :: bin

          integer :: n
          integer :: i_start
          integer :: t
          integer :: i
          integer :: r
          integer :: b


          bin = [0,0,0,0,0,0,0,0]

          if(config_id.ne.0) then
             n = floor(log(real(config_id))/log(real(2)))+1
             i_start = 8-n+1

             t = config_id
             do i=1, n

                r = MODULO(t,2**(n-i))
                b = nint(real(t-r)/real(2**(n-i)))
                t = r

                bin(i_start+i-1) = b

             end do             

          end if

        end function get_bin_decomposition

      end program test_bf_newgrdpt_procedure
