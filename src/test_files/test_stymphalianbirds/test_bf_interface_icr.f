      program test_bf_interface_icr

        use bf_interface_icr_class, only : bf_interface_icr

        implicit none

        type(bf_interface_icr) :: bf_interface_icr_used
        logical                :: detailled
        logical                :: test_validated
        
        detailled = .true.

        test_validated = test_get_central_grdpt(bf_interface_icr_used, detailled)
        print '(''test_bf_interface_icr: '',L1)', test_validated
        print '()'


        contains


        function test_get_central_grdpt(
     $       bf_interface_icr_used,
     $       detailled)
     $       result(test_validated)

          implicit none

          class(bf_interface_icr), intent(in) :: bf_interface_icr_used
          logical                , intent(in) :: detailled
          logical                             :: test_validated


          integer(ikind), dimension(2)  :: d_icoord
          real(rkind)   , dimension(2)  :: d_rcoord
          real(rkind)   , dimension(2)  :: velocity
          real(rkind)   , dimension(nx) :: interior_x_map
          real(rkind)   , dimension(ny) :: interior_y_map
          integer(ikind), dimension(2)  :: d_icoord_n
          real(rkind)   , dimension(2)  :: d_rcoord_n
          integer(ikind), dimension(2)  :: cpt_coord
          integer(ikind), dimension(2)  :: d_icoord_n_test
          real(rkind)   , dimension(2)  :: d_rcoord_n_test
          integer(ikind), dimension(2)  :: cpt_coord_test


          if(.not.is_test_validated(dt,2.0d0,.false.)) then
             print '(''test_bf_interface_icr'')'
             print '(''test_get_central_grdpt'')'
             print '(''the test requires dt=2.0'')'
             stop ''
          end if
          
          d_icoord            = [2,1]
          d_rcoord            = [0.5d0,0.0d0]
          velocity            = [3.0d0,1.0d0]
          interior_x_map(1:2) = [0.0d0,6.0d0]
          interior_y_map(1:2) = [0.0d0,3.0d0]
          d_icoord_n_test     = [3,1]
          d_rcoord_n_test     = [6.5d0,2.0d0]
          cpt_coord_test      = [6,2]

          cpt_coord = bf_interface_icr_used%get_central_grdpt(
     $         d_icoord,
     $         d_rcoord,
     $         velocity,
     $         interior_x_map,
     $         interior_y_map,
     $         d_icoord_n,
     $         d_rcoord_n)

          test_validated = .true.

          test_loc = d_icoord_n(1).eq.d_i_coord_n_test(1)
          if(detailled.and.(.not.test_loc)) then
             print '(''d_icoord_n(1): '',F6.2,'' -> '',F6.2)',
     $            d_icoord_n(1), d_icoord_n_test(1)
          end if
          test_validated = test_validated.and.test_loc

          test_loc = d_icoord_n(2).eq.d_icoord_n_test(2)
          if(detailled.and.(.not.test_loc)) then
             print '(''d_icoord_n(2): '',F6.2,'' -> '',F6.2)',
     $            d_icoord_n(2), d_icoord_n_test(2)
          end if
          test_validated = test_validated.and.test_loc

          test_loc = is_test_validated(d_rcoord_n(1),d_rcoord_n(1),detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_test_validated(d_rcoord_n(2),d_rcoord_n(2),detailled)
          test_validated = test_validated.and.test_loc          

        end function test_get_central_grdpt

      end program test_bf_interface_icr
