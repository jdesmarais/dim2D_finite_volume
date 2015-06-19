#!/usr/bin/python -i

import sys
import inspect
import os
import getopt
sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit
import time

from library_vtklines_to_curve_positions import get_curve_coords_from_vtkfile
from math import sqrt

import numpy as np

# add the python files from sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from library_sm_lg_inputs import (get_mass_density_vapor,
                                  get_mass_density_liquid,
                                  get_interface_length)

from automatization_contours_csts import cv_r,we


def create_mg_progress(mg_progress):
    '''
    @description:
    print a message which is overwritten
    '''
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


def create_mg_final(mg_progress):
    '''
    @description:
    print a message which is not overwritten
    '''
    sys.stdout.write('%s' % mg_progress)
    sys.stdout.flush()
    print '\n'


def compute_volume(x_data,y_data):
    '''
    @description: compute the volume of
    the contours using trapezoid rule
    '''

    y_min = min(y_data)

    x_prev = x_data[0]
    y_prev = y_data[0]

    volume = 0.0
        
    for i in range(1,len(x_data)):
            
        x = x_data[i]
        y = y_data[i]
            
        dx = x - x_prev
        dy = y - y_prev
            
        volume += dx*(y + 0.5*dy)
            
        x_prev = x
        y_prev = y

    volume-=y_min*(x_data[-1]-x_data[0])

    return volume


def get_time(ncPath):
    '''
    @description:
    extract the time from the database
    '''

    md = visit.GetMetaData(ncPath)
    return float(md.times[0])


def get_plot_range():
    '''
    @description:
    extract plot range of the current plot
    '''
    
    minVar = visit.Query('Min')
    maxVar = visit.Query('Max')

    minVar = float(minVar.split('=')[1].split('(')[0])
    maxVar = float(maxVar.split('=')[1].split('(')[0])

    return [minVar,maxVar]


def get_nb_nodes():
    '''
    @description:
    extract the number of nodes
    '''

    nbNodes = visit.Query('NumNodes')
    nbNodes = int(nbNodes.split('is')[1].split('.')[0])

    return nbNodes


def get_domain_borders():
    '''
    @description:
    extract the [x_min,x_max,y_min,y_max] of the
    domain plotted    
    '''
    
    visit.AddPlot("Curve", "x", 1, 1)
    visit.DrawPlots()
    [x_min,x_max] = get_plot_range()
    nx = get_nb_nodes()
    visit.DeleteActivePlots()

    visit.AddPlot("Curve", "y", 1, 1)
    visit.DrawPlots()
    [y_min,y_max] = get_plot_range()
    ny = get_nb_nodes()
    visit.DeleteActivePlots()

    dx = (x_max-x_min)/(nx-1)
    dy = (y_max-y_min)/(ny-1)

    domain_borders = {}
    domain_borders['nx']=nx
    domain_borders['ny']=ny
    domain_borders['dx']=dx
    domain_borders['dy']=dy
    domain_borders['x_min']=x_min
    domain_borders['x_max']=x_max
    domain_borders['y_min']=y_min
    domain_borders['y_max']=y_max    

    return domain_borders
    

def get_var_range(var):
    '''
    @description: extract the min,max,mid values
    of a variable by first plotting the pseudocolor
    plot, extracting the min and max values and 
    computing the middle value
    '''
    
    visit.AddPlot("Pseudocolor",var, 1, 1)
    visit.DrawPlots()

    [minVar,maxVar] = get_plot_range()

    visit.DeleteActivePlots()

    midVar = 0.5*(minVar+maxVar)

    return [minVar,maxVar,midVar]


def get_mid_by_max_grad(domain_borders,phase_check=False):
    '''
    @description: extract the mass density at the
    interface between the liquid and vapor phases
    by extracting the mass at the location where
    the maximum mass denisty gradient alogn the wall
    is reached    
    '''

    visit.DefineScalarExpression("mass_grad_x", "gradient(mass)[0]")
    visit.AddPlot("Pseudocolor","mass_grad_x", 1, 1)

    visit.AddOperator("Box",1)
    BoxAtts = visit.BoxAttributes()
    BoxAtts.amount = BoxAtts.All
    BoxAtts.minx = domain_borders['x_min']
    BoxAtts.maxx = domain_borders['x_max']
    BoxAtts.miny = domain_borders['y_min']+1.5*domain_borders['dy']
    BoxAtts.maxy = domain_borders['y_min']+3.5*domain_borders['dy']
    BoxAtts.inverse = 0
    visit.SetOperatorOptions(BoxAtts,1)

    visit.DrawPlots()

    gradxMass = visit.Query('Max', use_actual_data=1)
    coord1 = float(gradxMass.split(',')[0].split('<')[1])
    coord2 = float(gradxMass.split(',')[1].split('>')[0])

    grad_x_Mass = float(gradxMass.split('=')[1].split('(')[0])

    visit.DeleteActivePlots()

    linearInterpolation=(grad_x_Mass!=0.0)

    if(grad_x_Mass!=0.0):       

    	if(linearInterpolation):
    	    # two points are needed for the linear interpolation
    	    # these are the points next to the node where the maximum
    	    # mass gradient along x is reached
    	    coord1_int0 = coord1-domain_borders['dx']
    	    coord1_int1 = coord1+domain_borders['dx']
    	    
    	    # we need to find the coordinates where the maximum
    	    # gradient is reached by d2mass/dx2 = 0 using linear
    	    # interpolation between the two previous points
    	    visit.DefineScalarExpression("mass_grad_xx", "gradient(mass_grad_x)[0]")
    	    visit.AddPlot("Pseudocolor","mass_grad_xx", 1, 1)
    	    visit.DrawPlots()
    	
    	    massGradXX_0 = visit.NodePick(
    	        coord=(coord1_int0,coord2,0),
    	        vars=("default","mass_grad_xx"))
    	    massGradXX_0 = massGradXX_0['mass_grad_xx']
    	
    	    massGradXX_1 = visit.NodePick(
    	        coord=(coord1_int1,coord2,0),
    	        vars=("default","mass_grad_xx"))
    	    massGradXX_1 = massGradXX_1['mass_grad_xx']
    	
    	    coord1 = coord1_int0 - massGradXX_0*(coord1_int1-coord1_int0)/(massGradXX_1-massGradXX_0)
    	    visit.DeleteActivePlots()
    	
    	    # we need to determine the mass at the coord1
    	    # using linear interpolation
    	    visit.AddPlot("Pseudocolor","mass", 1, 1)
    	    visit.DrawPlots()
    	
    	    mass_0 = visit.NodePick(
    	        coord=(coord1_int0,coord2,0),
    	        vars=("default","mass"))
            mass_0 = mass_0['mass']
    	
    	    mass_1 = visit.NodePick(
    	        coord=(coord1_int1,coord2,0),
    	        vars=("default","mass"))
            mass_1 = mass_1['mass']
    	
    	    midMass = mass_0 + (mass_1-mass_0)/(coord1_int1-coord1_int0)*(coord1-coord1_int0)
    	
    	    visit.DeleteActivePlots()
    	
    	else:
    	
    	    midMass = visit.NodePick(coord=(coord1,coord2,0), vars=("default","mass"))        
    	    midMass = midMass['mass']
    	
    	if(phase_check):
    	
    	    # definition of the kortweg energy
    	    visit.DefineScalarExpression("mass_grad_y", "gradient(mass)[1]")
    	    visit.DefineScalarExpression("mass_grad_squared", "mass_grad_x^2+mass_grad_y^2")
    	    visit.DefineScalarExpression("we", str(we))
    	    visit.DefineScalarExpression("korteweg_energy", "0.5/we*mass_grad_squared")
    	
    	    # definition of the kinetic energy
    	    visit.DefineScalarExpression("velocity_x", "momentum_x/mass")
    	    visit.DefineScalarExpression("velocity_y", "momentum_y/mass")
    	    visit.DefineScalarExpression("kinetic_energy", "0.5*mass*(velocity_x^2+velocity_y^2)")
    	
    	    # definition of the temperature
    	    visit.DefineScalarExpression("cv_r", str(cv_r))
    	    visit.DefineScalarExpression("temperature", "3/(8*cv_r)*(1/mass*(energy-kinetic_energy-korteweg_energy)+3*mass)")
    	
    	    # draw the temperature and extract the temperature at the maximum gradient point
    	    visit.AddPlot("Pseudocolor","temperature", 1, 1)
    	    visit.DrawPlots()
    	    midTemperature = visit.NodePick(coord=(coord1,coord2,0), vars=("default","temperature"))
    	    midTemperature = midTemperature['temperature']
    	    visit.DeleteActivePlots()
    	
    	    # deduce the mass densities of the liquid and vapor phases at this temperature
    	    mass_vap = get_mass_density_vapor(midTemperature)
    	    mass_liq = get_mass_density_liquid(midTemperature)
    	
    	    # check whether the mass density extracted is close enough from the
    	    # mid mass density
    	    mid_mass_c = 0.5*(mass_vap+mass_liq)
    	    check = abs((midMass-mid_mass_c)/(mass_liq-mass_vap)) < 0.2
    	
    	    # check with the norm of the mass density graident
    	    visit.AddPlot("Pseudocolor","mass_grad_y", 1, 1)
    	    visit.DrawPlots()
    	    gradyMass = visit.NodePick(coord=(coord1,coord2,0), vars=("default","mass_grad_y"))
    	    grad_y_Mass = gradyMass['mass_grad_y']
    	    visit.DeleteActivePlots()
    	
    	    gradMassNorm = sqrt(grad_x_Mass**2+grad_y_Mass**2)
    	
    	    # check by comparing the maximum gradient of the mass density
    	    interface_lgh = get_interface_length(we,midTemperature)
    	    mid_gradMass = (mass_liq-mass_vap)/interface_lgh
    	    check = abs((gradMassNorm - mid_gradMass)/mid_gradMass) < 0.5
    	
    	else:
    	
    	    check = True
    	
    	visit.RemoveOperator(0, 1)
    	visit.DeleteActivePlots()
    	visit.ClearPickPoints()

    else:
        midMass = 0.0
        check = False

    return [midMass,check]


def add_contours(var,contourBorders):
    '''
    @description:
    add one single contour
    '''

    visit.AddPlot("Contour", var, 1, 1)

    ContourAtts = visit.ContourAttributes()
    ContourAtts.changedColors = ()
    ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ContourAtts.colorTableName = "Default"
    ContourAtts.invertColorTable = 0
    ContourAtts.legendFlag = 1
    ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    ContourAtts.lineWidth = 1
    ContourAtts.singleColor = (255, 0, 0, 255)
    ContourAtts.contourNLevels = 1
    ContourAtts.contourValue = ()
    ContourAtts.contourPercent = ()
    ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
    ContourAtts.minFlag = 1
    ContourAtts.maxFlag = 1
    ContourAtts.min = contourBorders[0]
    ContourAtts.max = contourBorders[1]
    ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
    ContourAtts.wireframe = 0
    visit.SetPlotOptions(ContourAtts)


def remove_boundary_pts(domain_borders,bc_size=2):
    '''
    @description:
    remove the boundary points from the
    computational domain
    '''

    x_min_R = domain_borders['x_min']+bc_size*domain_borders['dx']
    x_max_R = domain_borders['x_max']-bc_size*domain_borders['dx']
    y_min_R = domain_borders['y_min'] #+bc_size*domain_borders['dy']
    y_max_R = domain_borders['y_max']-bc_size*domain_borders['dy']

    visit.AddOperator("Box",1)
    BoxAtts = visit.BoxAttributes()
    BoxAtts.amount = BoxAtts.All #Some
    BoxAtts.minx = x_min_R
    BoxAtts.maxx = x_max_R
    BoxAtts.miny = y_min_R
    BoxAtts.maxy = y_max_R
    BoxAtts.inverse = 0
    visit.SetOperatorOptions(BoxAtts,1)


def add_reflection_x(x_reflection):

    visit.AddOperator("Reflect", 1)
    ReflectAtts = visit.ReflectAttributes()
    ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
    ReflectAtts.useXBoundary = 0
    ReflectAtts.specifiedX = x_reflection
    ReflectAtts.useYBoundary = 1
    ReflectAtts.specifiedY = 0
    ReflectAtts.useZBoundary = 1
    ReflectAtts.specifiedZ = 0
    ReflectAtts.reflections = (1, 1, 0, 0, 0, 0, 0, 0)
    visit.SetOperatorOptions(ReflectAtts, 1)


def generate_contours(var,domain_borders,contourBorders,reflection=False):
    '''
    @description: generate the contours
    corresponding to the isosurfaces
    of 'var' defined by contourBorders
    '''

    bc_size=2

    # plot the contour
    add_contours(var,contourBorders)

    # remove the boundary
    remove_boundary_pts(domain_borders,bc_size=bc_size)

    # if reflection added, add a reflection operator
    if(reflection):
        x_reflection = domain_borders['x_min']+bc_size*domain_borders['dx']
        add_reflection_x(x_reflection)

    visit.DrawPlots()


def export_contours_to_vtk(vtkPath):
    '''
    @description: export the data of the plot
    to vtk format
    '''

    outFile = os.path.basename(vtkPath)
    outDir  = os.path.dirname(vtkPath)

    ExportDBAtts = visit.ExportDBAttributes()
    ExportDBAtts.allTimes = 0
    ExportDBAtts.db_type = "VTK"
    ExportDBAtts.db_type_fullname = "VTK_1.0"
    ExportDBAtts.filename = outFile
    ExportDBAtts.dirname = outDir
    ExportDBAtts.variables = ()
    ExportDBAtts.opts.types = (0)
    visit.ExportDatabase(ExportDBAtts)


def generate_vtklines(ncPath,
                      contourRootPath,
                      var='mass',
                      contourPer=0.15,
                      contourType='mass',
                      phase_check=False,
                      time=0,
                      reflection=True):
    '''
    @description: generate the vtkfiles from a netcdf file
    '''

    visit.OpenDatabase("localhost:"+ncPath, 0)

    vtkPath = 'contours'

    domain_borders = get_domain_borders()

    if(contourType=='mass'):

        #extract [min,max,mid] for the mass
        [minVar,maxVar,midVar] = get_var_range(var)
        contourMin = midVar
        contourMax = midVar

    elif(contourType=='gradient'):

        #extract [min,max,mid] for the Korteweg energy (0.5/We*(dm)^2)
        visit.DefineScalarExpression("mass_grad_x", "gradient(mass)[0]")
        visit.DefineScalarExpression("mass_grad_y", "gradient(mass)[1]")
        visit.DefineScalarExpression("mass_grad_squared", "mass_grad_x^2+mass_grad_y^2")
        [minVar,maxVar,midVar] = get_var_range("mass_grad_squared")
        contourMin = maxVar*(1.0-contourPer)
        contourMax = contourMin
        var = 'mass_grad_squared'

    elif(contourType=='wall_max_gradient'):
                
        [contourMin,generateContours] = get_mid_by_max_grad(domain_borders,
                                                            phase_check=phase_check)
        contourMax = contourMin

    else:

        print 'library_nc_to_vtklines'
        print 'generate_vtklines'
        print 'contourType not recognized'
        sys.exit(2)

    # choose to plot the contours corresponding to 
    # a fixed mass density or corresponding to the 
    # maximum gradient of the mass density
    if(generateContours):

        generate_contours(var,domain_borders,[contourMin,contourMax],reflection=reflection)

        # export the contours to a vtk file
        export_contours_to_vtk(vtkPath)

        # get the coordinates of the contour
        if(reflection):
            
            graph_data0 = get_curve_coords_from_vtkfile(
                os.path.basename(vtkPath)+'.0.vtk')
            
            graph_data1 = get_curve_coords_from_vtkfile(
                os.path.basename(vtkPath)+'.1.vtk')
            
            graph_data = []
            graph_data.append([])
            graph_data.append([])
            graph_data[0] = graph_data1[0]+graph_data0[0]
            graph_data[0].append(-graph_data1[0][0])
            graph_data[1] = graph_data1[1]+graph_data0[1]
            graph_data[1].append(graph_data1[1][0])
            
        else:
            graph_data = get_curve_coords_from_vtkfile(
                os.path.basename(vtkPath)+'.vtk')
            
        outfile=contourRootPath+str(time)+'.curve'
        out = open(outfile, 'w')
        for (x,y) in zip(graph_data[0],graph_data[1]):
            out.write("%f %f\n" % (x,y))
        out.close()

        if(reflection):
            os.remove(os.path.basename(vtkPath)+'.0.vtk',)
            os.remove(os.path.basename(vtkPath)+'.1.vtk')
        else:
            sys.rm(os.path.basename(vtkPath)+'.vtk')

        #get the volume
        volume = compute_volume(graph_data[0][:],graph_data[1][:])

    else:
        #no contours
        graph_data='None'

        #no volume
        volume = 0.0

    visit.DeleteActivePlots()

    return [graph_data, volume, contourMin, domain_borders]


def generate_time_contour_data(ncRootPath,
                               contourRootPath,
                               timeRange=[0,1,1],
                               var='mass',
                               contourPer=0.15,
                               contourType='mass',
                               phase_check=False,
                               reflection=True):

    '''
    @description: extract the surface contact between the bubble
    and the wall as well as the volume of the bubble
    '''
    
    
    # create matrices of size timeRange[1]-timeRange[0]
    # that will contain the time, contact length and volume
    # of the bubble
    nt = int(float(timeRange[1]-timeRange[0])/float(timeRange[2]))

    file_id     = np.empty([nt])
    time        = np.empty([nt])
    contact_lgh = np.empty([nt])
    volume      = np.empty([nt])

    mg_progress = 'generating contour files: ...'
    create_mg_progress(mg_progress)
    vtkRootPath = contourRootPath


    # extract tha data as functions of time
    for i in range(0,nt):

        t = i*timeRange[2]

        ncPath  = ncRootPath+str(t)+'.nc'
        vtkPath = vtkRootPath+str(t)

        if(not os.path.isfile(ncPath)):
            break

        # write the file id
        file_id[i] = t

        # extract the time
        time[i] = get_time(ncPath)

        # extract the graph data
        [graph_data_t,
         volume_t,
         contour_t,
         domain_borders] = generate_vtklines(
            ncPath,
            contourRootPath,
            var=var,
            contourPer=contourPer,
            contourType=contourType,
            phase_check=phase_check,
            time=t,
            reflection=reflection)
        
        # save the volume at t
        volume[i]= volume_t
        
        # determine the contact length
        if(graph_data_t!='None'):

            # interpolate x=f(y) to find x such that y=0 for
            # the contact length
            y = np.array(graph_data_t[1][:])
            i_min, = np.unravel_index(y.argmin(),y.shape)
            y_min = y[i_min]
            if(y_min>(domain_borders['y_min'])): #+2*domain_borders['dy']
                contact_lgh[i]=0.0
            else:
                
                contact_lgh[i] = abs(graph_data_t[0][i_min])

        else:

            contact_lgh[i]=0.0
            

        mg_progress = 'generating contour files: '+str(i+1)+' / '+str(nt)
        create_mg_progress(mg_progress)

    mg_progress = 'generating contour files: done              '
    create_mg_final(mg_progress)


    # write the volume(t) on an output file
    out = open(os.path.dirname(contourRootPath)+'/volume.txt', 'w')
    for (i,t,v) in zip(file_id,time,volume):
        out.write("%f %f %f\n" % (i,t,v))
    out.close()


    # write the length(t) on an output file
    out = open(os.path.dirname(contourRootPath)+'/contact_lgh.txt', 'w')
    for (i,t,l) in zip(file_id,time,contact_lgh):
        out.write("%f %f %f\n" % (i,t,l))
    out.close()
