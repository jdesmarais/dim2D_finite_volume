bf_layer_computations
---------------------

 - bf_compute_class: main object encapsulating the time integration procedures for the buffer layer

 - bf_layer_bc_procedure_module: for the determination of which procedure should be used to apply the boundary conditions
 - bf_layer_bc_sections_class: for the localization of the boundary sections in the buffer layer
 - bf_interior_bc_sections_module: for the localization of the boundary sections in the interior domain

 - bf_layer_newgrdpt_procedure_module: for the determination of which procedure should be used to compute the new grid point
 - bf_layer_newgrdpt_class: for the computation of the new grid point

 - bf_bc_crenel_module: for the removal of bc crenels in the grdpts_id
 - bf_suspicious_bc_interior_pt_module: for the transformation of bc_interior_pt into interior_pt if all grdpts are available


bf_layer_sync
-------------

 - bf_layer_extract_module: for the extraction of grid points from the interior and the buffer layer
 - bf_layer_exchange_module: for the exchange of the grid points between the interior and the buffer layers


