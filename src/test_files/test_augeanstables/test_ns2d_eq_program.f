      program test_ns2d_eq_program

        use ns2d_parameters, only :
     $     viscous_r,
     $     Re, Pr, gamma,
     $     mach_infty

        use ns2d_prim_module, only :
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       qx_inviscid_x_flux,
     $       qy_inviscid_y_flux,
     $       qxy_transport,
     $       energy_inviscid_x_flux,
     $       energy_inviscid_y_flux,
     $       speed_of_sound,
     $       compute_jacobian_prim_to_cons,
     $       compute_jacobian_cons_to_prim,
     $       compute_cons_lodi_matrix_x,
     $       compute_cons_lodi_matrix_y

        use ns2d_fluxes_module, only :
     $       flux_x_mass_density, flux_y_mass_density,
     $       flux_x_momentum_x  , flux_y_momentum_x,
     $       flux_x_momentum_y  , flux_y_momentum_y,
     $       flux_x_total_energy, flux_y_total_energy

        use parameters_input, only : ne

        use parameters_kind, only : ikind, rkind

        use pmodel_eq_class, only : pmodel_eq

        use sd_operators_class, only : sd_operators


        implicit none

        real(rkind), dimension(5,5,4) :: nodes
        real(rkind)                   :: dx
        real(rkind)                   :: dy
        type(sd_operators)            :: s
        type(pmodel_eq)               :: p_model
        logical                       :: detailled
        logical                       :: test_validated

        detailled = .false.
        if(
     $       (.not.is_test_validated(viscous_r,-2.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(Re,10.0d0,detailled)).or.
     $       (.not.is_test_validated(Pr,1.0d0,detailled)).or.
     $       (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(mach_infty,1.0d0,detailled))) then

           print '(''the test requires: '')'
           print '(''viscous_r=-2/3'')'
           print '(''Re=10.0'')'
           print '(''Pr=1.0'')'
           print '(''gamma=5/3'')'
           print '(''mach_infty=1.0'')'
           stop ''

        end if


        !initialize the nodes
        call initialize_nodes(nodes,dx,dy)
        

        detailled = .false.


        !compute the initial conditions related quantities
        !for the NS equations and
        print '(''test ns2_initial_conditions'')'
        print '(''---------------------------'')'
        call test_ns2d_ic(p_model)
        print '(''---------------------------'')'
        print '()'


        !compute the quantities for the NS equations
        !and compare with the test data
        print '(''test ns2_prim'')'
        print '(''---------------------------'')'
        test_validated = test_ns2d_prim(nodes,detailled)
        if(detailled) print '(''---------------------------'')'
        print '(''test_validated: '',L3)', test_validated
        print '(''---------------------------'')'
        print '()'


        !compute the fluxes for the NS equations
        !and compare with the test data
        detailled = .false.

        print '(''test ns2_fluxes'')'
        print '(''---------------------------'')'
        test_validated = test_ns2d_fluxes(nodes,dx,dy,s,detailled)
        if(detailled) print '(''---------------------------'')'
        print '(''test_validated: '',L3)', test_validated
        print '(''---------------------------'')'
        print '()'


        !compute the eigenvalues and eigenmatrices
        !and compare with the test data
        detailled = .false.

        print '(''test ns2_eq'')'
        print '(''---------------------------'')'
        test_validated = test_ns2d_eq(nodes,p_model,detailled)
        if(detailled) print '(''---------------------------'')'
        print '(''test_validated: '',L3)', test_validated
        print '(''---------------------------'')'
        print '()'


        contains

        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated


        function is_matrix_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(ne,ne), intent(in) :: var
          real(rkind), dimension(ne,ne), intent(in) :: cst
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          integer :: i,j
          logical :: loc

          test_validated = .true.

          do j=1,ne
             do i=1,ne
                loc = is_test_validated(var(i,j),cst(i,j),detailled)
                test_validated = test_validated.and.loc
                if(detailled) print '(''var('',2I2,''): '',L1)', i,j,loc
             end do
          end do

        end function is_matrix_validated


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind)                  , intent(out) :: dx
          real(rkind)                  , intent(out) :: dy

          integer(ikind) :: j


          !space steps
          dx = 0.5
          dy = 0.6


          !mass
          nodes(1,1,1) = 0.5
          nodes(2,1,1) = 0.2
          nodes(3,1,1) = 1.2
          nodes(4,1,1) = 5.0
          nodes(5,1,1) = 0.6

          nodes(1,2,1) = 3.0
          nodes(2,2,1) = 4.2
          nodes(3,2,1) = 11.0
          nodes(4,2,1) = 10.6
          nodes(5,2,1) = 5.2

          nodes(1,3,1) = -14.2
          nodes(2,3,1) = 23
          nodes(3,3,1) = 9.8
          nodes(4,3,1) = 3.4
          nodes(5,3,1) = 9.1

          nodes(1,4,1) = 2.45
          nodes(2,4,1) = 0.2
          nodes(3,4,1) = 9.0
          nodes(4,4,1) = 5.4
          nodes(5,4,1) =-2.3

          nodes(1,5,1) = 3.6
          nodes(2,5,1) = 0.1
          nodes(3,5,1) = 6.3
          nodes(4,5,1) = 8.9
          nodes(5,5,1) = -4.23


          !momentum-x
          nodes(1,1,2) = 7.012
          nodes(2,1,2) =-6.323
          nodes(3,1,2) = 3.012
          nodes(4,1,2) = 4.5
          nodes(5,1,2) = 9.6
                    
          nodes(1,2,2) = 4.26
          nodes(2,2,2) = 4.23
          nodes(3,2,2) = 4.5
          nodes(4,2,2) = 7.56
          nodes(5,2,2) = 7.21
                    
          nodes(1,3,2) = 0.23
          nodes(2,3,2) = 7.23
          nodes(3,3,2) = 3.1
          nodes(4,3,2) = 8.9
          nodes(5,3,2) = 9.3
                    
          nodes(1,4,2) = 8.23
          nodes(2,4,2) = -3.1
          nodes(3,4,2) = 6.03
          nodes(4,4,2) = 6.25
          nodes(5,4,2) = 5.12
                    
          nodes(1,5,2) = 3.2
          nodes(2,5,2) = 8.12
          nodes(3,5,2) = 8.9
          nodes(4,5,2) = 4.2
          nodes(5,5,2) = 7.8


          !momentum-y
          nodes(1,1,3) = 7.1
          nodes(2,1,3) = 1.052
          nodes(3,1,3) = 1.23
          nodes(4,1,3) = 7.89
          nodes(5,1,3) = 8.0
                    
          nodes(1,2,3) = 8.362
          nodes(2,2,3) = 4.56
          nodes(3,2,3) = 9.6
          nodes(4,2,3) = 8.96
          nodes(5,2,3) = -3.23
                    
          nodes(1,3,3) = 2.53
          nodes(2,3,3) = -3.23
          nodes(3,3,3) = 7.25
          nodes(4,3,3) = 1.02
          nodes(5,3,3) = 9.26
                    
          nodes(1,4,3) = 8.965
          nodes(2,4,3) = 4.789
          nodes(3,4,3) = 4.56
          nodes(4,4,3) = 3.012
          nodes(5,4,3) = -1.45
                    
          nodes(1,5,3) = 6.26
          nodes(2,5,3) = 5.201
          nodes(3,5,3) = 2.03
          nodes(4,5,3) = 7.89
          nodes(5,5,3) = 9.889


          !total energy
          nodes(1,1,4) = 6.23
          nodes(2,1,4) = 4.12
          nodes(3,1,4) = -3.6
          nodes(4,1,4) = -6.52
          nodes(5,1,4) = 9.57
                    
          nodes(1,2,4) = -0.12
          nodes(2,2,4) = 8.2
          nodes(3,2,4) = 1.2
          nodes(4,2,4) = 7.89
          nodes(5,2,4) = 5.62
                    
          nodes(1,3,4) = -6.23
          nodes(2,3,4) = 6.201
          nodes(3,3,4) = 6.7
          nodes(4,3,4) = 4.12
          nodes(5,3,4) = 1.29
                    
          nodes(1,4,4) = 1.2
          nodes(2,4,4) = 7.958
          nodes(3,4,4) = 1
          nodes(4,4,4) = -5.62
          nodes(5,4,4) = 0.36
                    
          nodes(1,5,4) = 9.6
          nodes(2,5,4) = 6.12
          nodes(3,5,4) = 8.9
          nodes(4,5,4) = 8.95
          nodes(5,5,4) = 6.3

          print '()'
          print '(''mass_density'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,1)
          end do
          print '()'

          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,2)
          end do
          print '()'

          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,3)
          end do
          print '()'

          print '()'
          print '(''total energy'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,4)
          end do
          print '()'

        end subroutine initialize_nodes


        subroutine test_ns2d_ic(p_model)

          implicit none

          type(pmodel_eq), intent(in) :: p_model
          
          real(rkind) :: t,x,y

          print '(''get_mach_ux_infty: '', F8.3)', p_model%get_mach_ux_infty()
          print '(''get_mach_uy_infty: '', F8.3)', p_model%get_mach_uy_infty()
          print '(''get_u_in:  '', F8.3)', p_model%get_u_in(t,x,y)
          print '(''get_v_in:  '', F8.3)', p_model%get_v_in(t,x,y)
          print '(''get_T_in:  '', F8.3)', p_model%get_T_in(t,x,y)
          print '(''get_P_out: '', F8.3)', p_model%get_P_out(t,x,y)

        end subroutine test_ns2d_ic


        function test_ns2d_prim(nodes,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          real(rkind), dimension(14)    :: test_data
          real(rkind), dimension(ne,ne) :: test_data_jac_pv
          real(rkind), dimension(ne,ne) :: test_data_jac_vp
          real(rkind), dimension(ne,ne) :: test_data_lodi_matrix_x
          real(rkind), dimension(ne,ne) :: test_data_lodi_matrix_y
          logical                       :: loc
          integer(ikind)                :: i,j


          test_data(1)  =  9.8        !mass
          test_data(2)  =  3.1        !momentum-x
          test_data(3)  =  7.25       !momentum-y
          test_data(4)  =  6.7        !total energy
                        
          test_data(5)  =  0.31632653 !velocity-x
          test_data(6)  =  0.73979592 !velocity-y
          test_data(7)  =  2.35195578 !pressure
          test_data(8)  =  0.39999248 !temperature
          
          test_data(9)  =  3.33256802 !qx_inviscid_x_flux
          test_data(10) =  7.71547619 !qy_inviscid_y_flux
          test_data(11) =  2.29336734 !qxy_transport
          test_data(12) =  2.86337376 !energy_inviscid_x_flux
          test_data(13) =  6.69659994 !energy_inviscid_y_flux

          test_data(14) =  0.632449587!speed_of_sound
          
          !jac_prim_to_cons
          test_data_jac_pv(1,1) = 1.0
          test_data_jac_pv(2,1) = 0.0
          test_data_jac_pv(3,1) = 0.0
          test_data_jac_pv(4,1) = 0.0
                         
          test_data_jac_pv(1,2) =-0.032278217
          test_data_jac_pv(2,2) = 0.102040816
          test_data_jac_pv(3,2) = 0.0
          test_data_jac_pv(4,2) = 0.0
                         
          test_data_jac_pv(1,3) =-0.075489379
          test_data_jac_pv(2,3) = 0.0
          test_data_jac_pv(3,3) = 0.102040816
          test_data_jac_pv(4,3) = 0.0
                         
          test_data_jac_pv(1,4) = 0.215786825
          test_data_jac_pv(2,4) =-0.210884354
          test_data_jac_pv(3,4) =-0.493197279
          test_data_jac_pv(4,4) = 0.666666667

          !jac_cons_to_prim
          test_data_jac_vp(1,1) = 1.0
          test_data_jac_vp(2,1) = 0.0
          test_data_jac_vp(3,1) = 0.0
          test_data_jac_vp(4,1) = 0.0
                         
          test_data_jac_vp(1,2) = 0.316326531
          test_data_jac_vp(2,2) = 9.8
          test_data_jac_vp(3,2) = 0.0
          test_data_jac_vp(4,2) = 0.0
                         
          test_data_jac_vp(1,3) = 0.739795918
          test_data_jac_vp(2,3) = 0.0
          test_data_jac_vp(3,3) = 9.8
          test_data_jac_vp(4,3) = 0.0
                         
          test_data_jac_vp(1,4) = 0.323680237
          test_data_jac_vp(2,4) = 3.1
          test_data_jac_vp(3,4) = 7.25
          test_data_jac_vp(4,4) = 1.5

          !lodi_matrix_x
          test_data_lodi_matrix_x(1,1) =-0.075489379
          test_data_lodi_matrix_x(2,1) = 0.0
          test_data_lodi_matrix_x(3,1) = 0.102040816
          test_data_lodi_matrix_x(4,1) = 0.0
                    
          test_data_lodi_matrix_x(1,2) = 0.184205655
          test_data_lodi_matrix_x(2,2) = 0.210884354
          test_data_lodi_matrix_x(3,2) = 0.493197279
          test_data_lodi_matrix_x(4,2) =-0.666666667
                    
          test_data_lodi_matrix_x(1,3) = 0.415847409
          test_data_lodi_matrix_x(2,3) =-0.843333941
          test_data_lodi_matrix_x(3,3) =-0.493197279
          test_data_lodi_matrix_x(4,3) = 0.666666667
                    
          test_data_lodi_matrix_x(1,4) = 0.015726241
          test_data_lodi_matrix_x(2,4) = 0.421565233
          test_data_lodi_matrix_x(3,4) =-0.493197279
          test_data_lodi_matrix_x(4,4) = 0.666666667

          !lodi_matrix_y
          test_data_lodi_matrix_y(1,1) =-0.032278217
          test_data_lodi_matrix_y(2,1) = 0.102040816
          test_data_lodi_matrix_y(3,1) = 0.0
          test_data_lodi_matrix_y(4,1) = 0.0

          test_data_lodi_matrix_y(1,2) = 0.184205655
          test_data_lodi_matrix_y(2,2) = 0.210884354
          test_data_lodi_matrix_y(3,2) = 0.493197279
          test_data_lodi_matrix_y(4,2) =-0.666666667
                                
          test_data_lodi_matrix_y(1,3) = 0.683670448
          test_data_lodi_matrix_y(2,3) =-0.210884354
          test_data_lodi_matrix_y(3,3) =-1.125646866
          test_data_lodi_matrix_y(4,3) = 0.666666667
                                
          test_data_lodi_matrix_y(1,4) =-0.252096798
          test_data_lodi_matrix_y(2,4) =-0.210884354
          test_data_lodi_matrix_y(3,4) = 0.139252308
          test_data_lodi_matrix_y(4,4) = 0.666666667


          i = 3
          j = 3


          test_validated = .true.
          loc = is_test_validated(
     $         mass_density(nodes,i,j),
     $         test_data(1),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test mass: '', L3)', loc


          loc = is_test_validated(
     $         momentum_x(nodes,i,j),
     $         test_data(2),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test momentum-x: '', L3)', loc


          loc = is_test_validated(
     $         momentum_y(nodes,i,j),
     $         test_data(3),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test momentum-y: '', L3)', loc


          loc = is_test_validated(
     $         total_energy(nodes,i,j),
     $         test_data(4),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test total_energy: '', L3)', loc


          loc = is_test_validated(
     $         velocity_x(nodes,i,j),
     $         test_data(5),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test velocity_x: '', L3)', loc


          loc = is_test_validated(
     $         velocity_y(nodes,i,j),
     $         test_data(6),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test velocity_y: '', L3)', loc


          loc = is_test_validated(
     $         pressure(nodes,i,j),
     $         test_data(7),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test pressure: '', L3)', loc


          loc = is_test_validated(
     $         temperature(nodes,i,j),
     $         test_data(8),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test temperature: '', L3)', loc


          loc = is_test_validated(
     $         qx_inviscid_x_flux(nodes,i,j),
     $         test_data(9),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qx_inviscid_x_flux: '', L3)', loc


          loc = is_test_validated(
     $         qy_inviscid_y_flux(nodes,i,j),
     $         test_data(10),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qy_inviscid_y_flux: '', L3)', loc


          loc = is_test_validated(
     $         qxy_transport(nodes,i,j),
     $         test_data(11),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qxy_transport: '', L3)', loc


          loc = is_test_validated(
     $         energy_inviscid_x_flux(nodes,i,j),
     $         test_data(12),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test energy_inviscid_x_flux: '', L3)', loc


          loc = is_test_validated(
     $         energy_inviscid_y_flux(nodes,i,j),
     $         test_data(13),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test energy_inviscid_y_flux: '', L3)', loc

          loc = is_test_validated(
     $         speed_of_sound(nodes(i,j,:)),
     $         test_data(14),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test speed_of_sound: '', L3)', loc

          loc = is_matrix_validated(
     $         compute_jacobian_prim_to_cons(nodes(i,j,:)),
     $         test_data_jac_pv,
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test jac_prim_to_cons: '', L3)', loc

          loc = is_matrix_validated(
     $         compute_jacobian_cons_to_prim(nodes(i,j,:)),
     $         test_data_jac_vp,
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test jac_cons_to_prim: '', L3)', loc

          loc = is_matrix_validated(
     $         compute_cons_lodi_matrix_x(nodes(i,j,:)),
     $         test_data_lodi_matrix_x,
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test cons_lodi_matrix_x: '', L3)', loc

          loc = is_matrix_validated(
     $         compute_cons_lodi_matrix_y(nodes(i,j,:)),
     $         test_data_lodi_matrix_y,
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test cons_lodi_matrix_y: '', L3)', loc

        end function test_ns2d_prim


        function test_ns2d_fluxes(nodes,dx,dy,s,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          type(sd_operators)           , intent(in) :: s
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          real(rkind), dimension(8) :: test_data
          logical                   :: loc
          integer(ikind)            :: i,j


          test_data(1)  =  5.165      !flux_x_mass
          test_data(2)  =  5.039545704!flux_x_momentum_x
          test_data(3)  =  1.139892201!flux_x_momentum_y
          test_data(4)  =  3.210322639!flux_x_total_energy

          test_data(5)  =  8.425      !flux_y_mass
          test_data(6)  =  3.115780591!flux_y_momentum_x
          test_data(7)  =  6.840166285!flux_y_momentum_y
          test_data(8)  =  2.615690238!flux_y_total_energy

          i = 3
          j = 3


          !fluxes-x
          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_mass_density(nodes,s,i,j),
     $         test_data(1),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_mass: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_momentum_x(nodes,s,i,j,dx,dy),
     $         test_data(2),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_momentum_x: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_momentum_y(nodes,s,i,j,dx,dy),
     $         test_data(3),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_momentum_y: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_total_energy(nodes,s,i,j,dx,dy),
     $         test_data(4),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_total_energy: '', L3)', loc


          !fluxes-y
          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_mass_density(nodes,s,i,j),
     $         test_data(5),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_mass: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_momentum_x(nodes,s,i,j,dx,dy),
     $         test_data(6),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_momentum_x: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_momentum_y(nodes,s,i,j,dx,dy),
     $         test_data(7),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_momentum_y: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_total_energy(nodes,s,i,j,dx,dy),
     $         test_data(8),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_total_energy: '', L3)', loc

        end function test_ns2d_fluxes


        function test_ns2d_eq(nodes,p_model,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          type(pmodel_eq)              , intent(in) :: p_model
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          real(rkind), dimension(ne,ne) :: test_data
          real(rkind), dimension(ne)    :: eigenvalues
          real(rkind), dimension(ne,ne) :: eigenmatrix
          logical                       :: loc
          logical                       :: detailled_loc
          integer(ikind)                :: i,j
          
          i = 3
          j = 3

          test_validated = .true.


          !test eigenvalues x
          print '(''x_eigenvalues'')'

          test_data(1,1) =  0.3163265
          test_data(2,1) =  0.3163265
          test_data(3,1) = -0.3161230
          test_data(4,1) =  0.9487761

          detailled_loc = detailled
          eigenvalues = p_model%compute_x_eigenvalues(nodes(i,j,:))
          loc = test_eigenvalues(eigenvalues, test_data(:,1),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'


          !test eigenvalues y
          print '(''y_eigenvalues'')'

          test_data(1,1) =  0.7397959
          test_data(2,1) =  0.7397959
          test_data(3,1) =  0.1073463
          test_data(4,1) =  1.3722455

          detailled_loc = detailled
          eigenvalues = p_model%compute_y_eigenvalues(nodes(i,j,:))
          loc = test_eigenvalues(eigenvalues, test_data(:,1),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'


          !test x_left_eigenmatrix
          print '(''x_left_eigenmatrix'')'

          test_data(1,1) = -0.0754837
          test_data(2,1) =  0
          test_data(3,1) =  0.1020408
          test_data(4,1) =  0

          test_data(1,2) =  0.4605227
          test_data(2,2) =  0.5272207
          test_data(3,2) =  1.2330163
          test_data(4,2) = -1.6666980

          test_data(1,3) =  0.2079237
          test_data(2,3) = -0.4216669
          test_data(3,3) = -0.2465986
          test_data(4,3) =  0.3333333
          
          test_data(1,4) =  0.0078631
          test_data(2,4) =  0.2107826
          test_data(3,4) = -0.2465986
          test_data(4,4) =  0.3333333

          detailled_loc = detailled
          eigenmatrix = p_model%compute_x_lefteigenvector(nodes(i,j,:))
          loc = test_eigenmatrix(eigenmatrix, test_data, detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'


          !test x_right_eigenmatrix
          print '(''x_right_eigenmatrix'')'

          test_data(1,1) =  0
          test_data(2,1) =  1
          test_data(3,1) =  2.5000470
          test_data(4,1) =  2.5000470

          test_data(1,2) =  0
          test_data(2,2) =  0.3163265
          test_data(3,2) = -0.7903224
          test_data(4,2) =  2.3719848

          test_data(1,3) =  9.8
          test_data(2,3) =  0.7397959
          test_data(3,3) =  1.8495245
          test_data(4,3) =  1.8495245
         
          test_data(1,4) =  7.25
          test_data(2,4) =  0.3236802
          test_data(3,4) =  1.8090549
          test_data(4,4) =  2.8093766

          detailled_loc = detailled
          eigenmatrix = p_model%compute_x_righteigenvector(nodes(i,j,:))
          loc = test_eigenmatrix(eigenmatrix, test_data, detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'


          !test y_left_eigenmatrix
          print '(''y_left_eigenmatrix'')'

          test_data(1,1) = -0.0322782
          test_data(2,1) =  0.1020408
          test_data(3,1) =  0
          test_data(4,1) =  0

          test_data(1,2) =  0.4605227
          test_data(2,2) =  0.5272207
          test_data(3,2) =  1.2330163
          test_data(4,2) = -1.6666980

          test_data(1,3) =  0.3418352
          test_data(2,3) = -0.1054421
          test_data(3,3) = -0.5628234
          test_data(4,3) =  0.3333333
         
          test_data(1,4) = -0.1260483
          test_data(2,4) = -0.1054421
          test_data(3,4) =  0.0696261
          test_data(4,4) =  0.3333333

          detailled_loc = detailled
          eigenmatrix = p_model%compute_y_lefteigenvector(nodes(i,j,:))
          loc = test_eigenmatrix(eigenmatrix, test_data, detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'


          !test y_right_eigenmatrix
          print '(''y_right_eigenmatrix'')'

          test_data(1,1) =  0
          test_data(2,1) =  1
          test_data(3,1) =  2.50004700
          test_data(4,1) =  2.50004700

          test_data(1,2) =  9.8
          test_data(2,2) =  0.3163265
          test_data(3,2) =  0.7908311
          test_data(4,2) =  0.7908311 

          test_data(1,3) =  0
          test_data(2,3) =  0.7397959
          test_data(3,3) =  0.2683708
          test_data(4,3) =  3.4306782
         
          test_data(1,4) =  3.1
          test_data(2,4) =  0.3236802
          test_data(3,4) =  1.1394847
          test_data(4,4) =  3.4789468

          detailled_loc = detailled
          eigenmatrix = p_model%compute_y_righteigenvector(nodes(i,j,:))
          loc = test_eigenmatrix(eigenmatrix, test_data, detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '', L3)', loc
          print '()'

        end function test_ns2d_eq


        function test_eigenvalues(eigenvalues,test_data,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(ne), intent(in) :: eigenvalues
          real(rkind), dimension(ne), intent(in) :: test_data
          logical                   , intent(in) :: detailled
          logical                                :: test_validated

          integer :: i
          logical :: loc

          test_validated = .true.

          do i=1,ne

             loc = is_test_validated(eigenvalues(i), test_data(i), detailled)
             test_validated = test_validated.and.loc
             if(detailled) print '(''eigenvalue('',I2,''): '',L3)', i, loc

          end do

        end function test_eigenvalues


        function test_eigenmatrix(eigenmatrix,test_data,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(ne,ne), intent(in) :: eigenmatrix
          real(rkind), dimension(ne,ne), intent(in) :: test_data
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          integer :: i,j
          logical :: loc

          test_validated = .true.

          do j=1,ne
             do i=1,ne
                
                loc = is_test_validated(eigenmatrix(i,j), test_data(i,j), detailled)
                test_validated = test_validated.and.loc
                if(detailled) then
                   print '(''eigenmatrix('',I2,'','',I2,''): '',L3)', i, j, loc
                end if
                
             end do
          end do

        end function test_eigenmatrix

      end program test_ns2d_eq_program
