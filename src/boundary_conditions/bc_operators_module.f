      !> @file
      !> determine how the boudnary conditions
      !> are applied
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> determine whether the boundary conditions
      !> are applied on the fluxes, nodes or timedev
      !
      !> @date
      !> 03_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       module bc_operators_module

        use parameters_constant, only :
     $       bc_nodes_choice,
     $       bc_fluxes_choice,
     $       bc_timedev_choice,
     $       bc_flux_and_node_choice

        use parameters_input, only :
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice
        
        implicit none

        private
        public :: 
     $       shall_bc_on_nodes_be_applied,
     $       shall_bc_on_fluxes_be_applied,
     $       shall_bc_on_timedev_be_applied

        contains


        function shall_bc_on_nodes_be_applied()

          implicit none

          logical :: shall_bc_on_nodes_be_applied

          shall_bc_on_nodes_be_applied =
     $         (bc_N_type_choice.eq.bc_nodes_choice).or.
     $         (bc_S_type_choice.eq.bc_nodes_choice).or.
     $         (bc_E_type_choice.eq.bc_nodes_choice).or.
     $         (bc_W_type_choice.eq.bc_nodes_choice).or.
     $         (bc_N_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_S_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_E_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_W_type_choice.eq.bc_flux_and_node_choice)

        end function shall_bc_on_nodes_be_applied


        function shall_bc_on_fluxes_be_applied()
     
          implicit none

          logical :: shall_bc_on_fluxes_be_applied

          shall_bc_on_fluxes_be_applied =
     $         (bc_N_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_S_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_E_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_W_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_N_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_S_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_E_type_choice.eq.bc_flux_and_node_choice).or.
     $         (bc_W_type_choice.eq.bc_flux_and_node_choice)

        end function shall_bc_on_fluxes_be_applied


        function shall_bc_on_timedev_be_applied()

          implicit none

          logical :: shall_bc_on_timedev_be_applied

          shall_bc_on_timedev_be_applied =
     $       (bc_N_type_choice.eq.bc_timedev_choice).or.
     $       (bc_S_type_choice.eq.bc_timedev_choice).or.
     $       (bc_E_type_choice.eq.bc_timedev_choice).or.
     $       (bc_W_type_choice.eq.bc_timedev_choice)

        end function shall_bc_on_timedev_be_applied


      end module bc_operators_module
