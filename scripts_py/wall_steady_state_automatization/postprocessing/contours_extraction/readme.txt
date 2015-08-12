============================================================
Author: J.L. Desmarais (desmaraisjulien@gmail.com)
Date  : 16/06/2015
============================================================
python files written to create and extract the bubble
contours from netcdf files

procedure sumary:
-----------------

library_nc_to_vtklines.py :
	The data*.nc files are plotted using visit (mass contours).
	The contours in visit are exported as vtk files which are
	then reconstructed as continuous lines by the library
	(library_vtklines_to_curve_positions.py). These data are
	analyzed to find the contact length betwen the wall and the
	bubble as well as the bubble volume.	

library_vtklines_to_curve_positions.py:
	The vtk files are analyzed to extract the polylines.

library_contours_graph.py:
	The contours extracted, the contact length and the bubble
	volume are plotted as function of time. They are compared
	to the spherical cap approximation.

extract_contours.py:
	This is the script which is used to run the postprocessing
	of the bubble steady state contours.

