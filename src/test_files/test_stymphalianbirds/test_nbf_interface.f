      program test_nbf_interface

        use nbf_interface_class, only :
     $     nbf_interface

        implicit none




        

        contains


        function test_compute_newgrdpt(detailled)
     $       result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface) :: nbf_interface_used


          !initialize the interior nodes
          
          
          !initialize two buffer layers

          !initialize the links in the nbf_interface

          !test the computation of the new grdpt
          ! - using only a buffer layer
          ! - using the temporary tables
          

        end function test_compute_newgrdpt

      end program test_nbf_interface
