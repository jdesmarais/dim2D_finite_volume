import sys
import os
import getopt
sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit
import shutil

from library_messages import\
    print_mg_progress,\
    print_mg_final

import time
import math

import numpy as np


def define_dim_variables(cv_r=3.05,we=83634):
    """
    Define the extra variables for displaying relevant quantities
    when analysing the results of Diffuse Interface Model simulations
    """

    visit.DefineScalarExpression("mass_grad_x", "gradient(mass)[0]")
    visit.DefineScalarExpression("mass_grad_y", "gradient(mass)[1]")
    visit.DefineScalarExpression("mass_grad_squared", "mass_grad_x^2+mass_grad_y^2")
    visit.DefineScalarExpression("cv_r", str(cv_r))
    visit.DefineScalarExpression("nordstrom", "pressure+mass^2*(-3+2*mass)")
    visit.DefineScalarExpression("pressure", "8*mass*temperature/(3-mass)-3*mass^2")
    visit.DefineScalarExpression("kinetic_energy", "0.5*mass*(velocity_x^2+velocity_y^2)")
    visit.DefineScalarExpression("velocity_x", "momentum_x/mass")
    visit.DefineScalarExpression("velocity_y", "momentum_y/mass")
    visit.DefineVectorExpression("velocity", "{velocity_x,velocity_y}")
    visit.DefineScalarExpression("we", str(we))
    visit.DefineScalarExpression("korteweg_energy", "0.5/we*mass_grad_squared")
    visit.DefineScalarExpression("temperature", "3/(8*cv_r)*(1/mass*(energy-kinetic_energy-korteweg_energy)+3*mass)")
    visit.DefineScalarExpression("norm_temp_grad", "sqrt(gradient(temperature)[0]^2+gradient(temperature)[1]^2)")
    visit.DefineScalarExpression("norm_velocity", "sqrt(velocity_x^2+velocity_y^2)")
    visit.DefineScalarExpression("div_velocity", "gradient(velocity_x)[0]+gradient(velocity_y)[1]")

    return


def set_default_pseudocolor_param(PseudocolorAtts):
    """
    Set the default values for the parameters of
    the pseudocolor plot    
    """
    
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.lightingFlag = 1
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 1

    return


def set_default_contour_param(ContourAtts):
    """
    Set the default values for the parameters of
    the contour plot
    """
    
    ContourAtts.defaultPalette.GetControlPoints(0).colors = (255, 0, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(0).position = 0
    ContourAtts.defaultPalette.GetControlPoints(1).colors = (0, 255, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(1).position = 0.034
    ContourAtts.defaultPalette.GetControlPoints(2).colors = (0, 0, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(2).position = 0.069
    ContourAtts.defaultPalette.GetControlPoints(3).colors = (0, 255, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(3).position = 0.103
    ContourAtts.defaultPalette.GetControlPoints(4).colors = (255, 0, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(4).position = 0.138
    ContourAtts.defaultPalette.GetControlPoints(5).colors = (255, 255, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(5).position = 0.172
    ContourAtts.defaultPalette.GetControlPoints(6).colors = (255, 135, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(6).position = 0.207
    ContourAtts.defaultPalette.GetControlPoints(7).colors = (255, 0, 135, 255)
    ContourAtts.defaultPalette.GetControlPoints(7).position = 0.241
    ContourAtts.defaultPalette.GetControlPoints(8).colors = (168, 168, 168, 255)
    ContourAtts.defaultPalette.GetControlPoints(8).position = 0.276
    ContourAtts.defaultPalette.GetControlPoints(9).colors = (255, 68, 68, 255)
    ContourAtts.defaultPalette.GetControlPoints(9).position = 0.31
    ContourAtts.defaultPalette.GetControlPoints(10).colors = (99, 255, 99, 255)
    ContourAtts.defaultPalette.GetControlPoints(10).position = 0.345
    ContourAtts.defaultPalette.GetControlPoints(11).colors = (99, 99, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(11).position = 0.379
    ContourAtts.defaultPalette.GetControlPoints(12).colors = (40, 165, 165, 255)
    ContourAtts.defaultPalette.GetControlPoints(12).position = 0.414
    ContourAtts.defaultPalette.GetControlPoints(13).colors = (255, 99, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(13).position = 0.448
    ContourAtts.defaultPalette.GetControlPoints(14).colors = (255, 255, 99, 255)
    ContourAtts.defaultPalette.GetControlPoints(14).position = 0.483
    ContourAtts.defaultPalette.GetControlPoints(15).colors = (255, 170, 99, 255)
    ContourAtts.defaultPalette.GetControlPoints(15).position = 0.517
    ContourAtts.defaultPalette.GetControlPoints(16).colors = (170, 79, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(16).position = 0.552
    ContourAtts.defaultPalette.GetControlPoints(17).colors = (150, 0, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(17).position = 0.586
    ContourAtts.defaultPalette.GetControlPoints(18).colors = (0, 150, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(18).position = 0.621
    ContourAtts.defaultPalette.GetControlPoints(19).colors = (0, 0, 150, 255)
    ContourAtts.defaultPalette.GetControlPoints(19).position = 0.655
    ContourAtts.defaultPalette.GetControlPoints(20).colors = (0, 109, 109, 255)
    ContourAtts.defaultPalette.GetControlPoints(20).position = 0.69
    ContourAtts.defaultPalette.GetControlPoints(21).colors = (150, 0, 150, 255)
    ContourAtts.defaultPalette.GetControlPoints(21).position = 0.724
    ContourAtts.defaultPalette.GetControlPoints(22).colors = (150, 150, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(22).position = 0.759
    ContourAtts.defaultPalette.GetControlPoints(23).colors = (150, 84, 0, 255)
    ContourAtts.defaultPalette.GetControlPoints(23).position = 0.793
    ContourAtts.defaultPalette.GetControlPoints(24).colors = (160, 0, 79, 255)
    ContourAtts.defaultPalette.GetControlPoints(24).position = 0.828
    ContourAtts.defaultPalette.GetControlPoints(25).colors = (255, 104, 28, 255)
    ContourAtts.defaultPalette.GetControlPoints(25).position = 0.862
    ContourAtts.defaultPalette.GetControlPoints(26).colors = (0, 170, 81, 255)
    ContourAtts.defaultPalette.GetControlPoints(26).position = 0.897
    ContourAtts.defaultPalette.GetControlPoints(27).colors = (68, 255, 124, 255)
    ContourAtts.defaultPalette.GetControlPoints(27).position = 0.931
    ContourAtts.defaultPalette.GetControlPoints(28).colors = (0, 130, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(28).position = 0.966
    ContourAtts.defaultPalette.GetControlPoints(29).colors = (130, 0, 255, 255)
    ContourAtts.defaultPalette.GetControlPoints(29).position = 1
    ContourAtts.defaultPalette.smoothing = ContourAtts.defaultPalette.None  # None, Linear, CubicSpline
    ContourAtts.defaultPalette.equalSpacingFlag = 1
    ContourAtts.defaultPalette.discreteFlag = 1
    ContourAtts.defaultPalette.categoryName = "Standard"
    ContourAtts.changedColors = ()
    ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ContourAtts.colorTableName = "Default"
    ContourAtts.invertColorTable = 0
    ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    ContourAtts.SetMultiColor(0, (255, 0, 0, 255))
    ContourAtts.SetMultiColor(1, (0, 255, 0, 255))
    ContourAtts.SetMultiColor(2, (0, 0, 255, 255))
    ContourAtts.SetMultiColor(3, (0, 255, 255, 255))
    ContourAtts.SetMultiColor(4, (255, 0, 255, 255))
    ContourAtts.SetMultiColor(5, (255, 255, 0, 255))
    ContourAtts.SetMultiColor(6, (255, 135, 0, 255))
    ContourAtts.SetMultiColor(7, (255, 0, 135, 255))
    ContourAtts.SetMultiColor(8, (168, 168, 168, 255))
    ContourAtts.SetMultiColor(9, (255, 68, 68, 255))
    ContourAtts.contourValue = ()
    ContourAtts.contourPercent = ()
    ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
    ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
    ContourAtts.wireframe = 0

    ContourAtts.contourNLevels = 1

    return


def plot_pseudocolor(filePath,
                     field,
                     colorTableName='Hotandcold',
                     fieldMin='None',
                     fieldMax='None',
                     legend=True,
                     scaleType='Linear'):
    """
    Plot a pseudocolor field

    filePath       (String)  : path for the file
    field          (String)  : name of the field plotted
    colorTableName (String)  : color table
    fieldMin       (Float)   : min for the colorTable
    fieldMax       (Float)   : max for the colorTable
    legend         (Boolean) : display the legend
    scaleType      (String)  : either Linear, Log or Skew
    """

    # open the file
    visit.OpenDatabase('localhost:'+filePath, 0)

    # add the field
    visit.AddPlot("Pseudocolor", field, 1, 1)

    # load the existing parameters for the field
    PseudocolorAtts = visit.PseudocolorAttributes()

    # set the default parameters for plotting the field
    set_default_pseudocolor_param(PseudocolorAtts)

    # set the user-defined parameters for plotting the field
    # Linear, Log, Skew
    dictScaleType = {'Linear':PseudocolorAtts.Linear,
                     'Log':PseudocolorAtts.Log,
                     'Skew':PseudocolorAtts.Skew}

    PseudocolorAtts.scaling = dictScaleType[scaleType]

    PseudocolorAtts.colorTableName = colorTableName

    if(fieldMin!='None'):
        PseudocolorAtts.minFlag = 1
        PseudocolorAtts.min = fieldMin
    else:
        PseudocolorAtts.minFlag = 0
        
    if(fieldMax!='None'):
        PseudocolorAtts.maxFlag = 1
        PseudocolorAtts.max = fieldMax
    else:
        PseudocolorAtts.maxFlag = 0
    
    if(legend):
        PseudocolorAtts.legendFlag = 1
    else:
        PseudocolorAtts.legendFlag = 0

    # apply the plot parameters
    visit.SetPlotOptions(PseudocolorAtts)

    return


def plot_contour(filePath,
                 field,
                 contourColor=(255, 0, 0, 255),
                 contourWidth=1,
                 contourValue='None',
                 legend=True):
    """
    Plot a contour field

    filePath       (String)  : path for the file
    field          (String)  : name of the field plotted
    contourColor   (Tuple)   : (R,G,B,max)
    contourWidth   (Integer) : line thickness
    contourValue   (Float)   : value corresponding to the contour
    legend         (Boolean) : display the legend
    """

    # open the file
    visit.OpenDatabase('localhost:'+filePath, 0)

    # add the field
    visit.AddPlot("Contour", field, 1, 1)    
    
    # load the existing parameters for the field
    ContourAtts = visit.ContourAttributes()

    # set the default parameters for plotting the field
    set_default_contour_param(ContourAtts)

    # set the user-defined parameters for plotting the field
    ContourAtts.singleColor = contourColor
    ContourAtts.lineWidth   = contourWidth

    if(contourValue!='None'):
        ContourAtts.minFlag = 1
        ContourAtts.maxFlag = 1
        ContourAtts.min = contourValue
        ContourAtts.max = contourValue

    if(legend):
        ContourAtts.legendFlag = 1
    else:
        ContourAtts.legendFlag = 0

    # set the contour options
    visit.SetPlotOptions(ContourAtts)

    return


def set_2D_view(x_min,x_max,y_min,y_max):
    """
    Set the 2D view by constraining
    [x_min,x_max] x [y_min,y_max]
    """

    View2DAtts = visit.View2DAttributes()
    View2DAtts.windowCoords = (x_min, x_max, y_min, y_max)
    View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 1
    visit.SetView2D(View2DAtts)

    return


def launch_visit(window=False):
    """
    Launch visit w/ or w/o window

    window (Boolean) : display window
    """

    if(window):
        visit.Launch()
    else:
        visit.LaunchNowin()

    #visit.SuppressMessages(1)

    return


def set_visit_save_options(outputDir='None',
                           multipleWindows=True,
                           windowArray=[1,1]):
    """
    Set the save options for visit

    outputDir       (String)            : directory path where pictures are saved
    multipleWindows (Boolean)           : whether multiple windows are saved
    windowArray     ([Integer,Integer]) : how windows are arranged
    """

    # import the save options
    SaveWindowAtts = visit.SaveWindowAttributes()

    # set the default save options
    SaveWindowAtts.fileName = "visit"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0

    # set the user-defined save options
    if(outputDir=='None' or not os.path.isdir(outputDir)):
        SaveWindowAtts.outputToCurrentDirectory = 1
        SaveWindowAtts.outputDirectory = "."
    else:
        SaveWindowAtts.outputToCurrentDirectory = 0
        SaveWindowAtts.outputDirectory = outputDir

    if(multipleWindows):
        SaveWindowAtts.saveTiled = 1
    else:
        SaveWindowAtts.saveTiled = 0

    width  = windowArray[0]*1024
    height = windowArray[1]*1024

    SaveWindowAtts.width = width
    SaveWindowAtts.height = height


    # set the save options
    visit.SetSaveWindowAttributes(SaveWindowAtts)

    return


def save_visit_windows():
    """
    Save the graphs generated by visit
    """
    visit.SaveWindow()


def display_DIM_soundwave_with_error(dirNc,
                                     timestep,
                                     soundMin,
                                     soundMax,
                                     errorMin,
                                     errorMax,
                                     bubbleValue,
                                     maxBfLayers=10,
                                     borders='None',
                                     legend=True,
                                     windows=[1,2]):
    """
    Display two windows:
     - the first window displays the sound waves
       as the velocity divergence field with color
       table Grays and the bubble contour as a thick red
       line
     - the second window displays the local error compared
       to the large domain simulation

    dirNc       ([String]) : paths to the folders with .nc files
    timestep    (Integer)  : time step plotted
    soundMin    (Float)    : minimum for the sound wave field
    soundMax    (Float)    : maximum for the sound wave field
    bubbleValue (Float)    : mass value identifying the bubble contour
    maxBfLayers (Integer)  : maximum number of buffer layers per cardinal pt
    borders     (Integerx4): [x_min,x_max,y_min,y_max] for the view
    legend      (Boolean)  : display graph legend
    """
    

    # first window
    #--------------------
    visit.SetActiveWindow(windows[0])

    nbActivePlots=[0,0]
    
    fieldPseudoColor='div_velocity'
    colorTableName='Greys'
    fieldMin=soundMin
    fieldMax=soundMax

    fieldContour='mass'
    contourWidth=3
    contourValue=bubbleValue
    

    # plot the velocity divergence field + mass contours
    # for the interior field
    filePath = os.path.join(dirNc,'sm_domain','data'+str(timestep)+'.nc')

    plot_pseudocolor(filePath,
                     fieldPseudoColor,
                     colorTableName=colorTableName,
                     fieldMin=fieldMin,
                     fieldMax=fieldMax,
                     legend=legend)
    
    plot_contour(filePath,
                 fieldContour,
                 contourColor=(255, 0, 0, 255),
                 contourWidth=contourWidth,
                 contourValue=contourValue,
                 legend=legend)

    nbActivePlots[0]+=2


    # plot the velocity divergence field + mass contours
    # for the buffer layers
    for cardinalPt in ['N','S','E','W']:

        for bufferNb in range(1,maxBfLayers):

            filePath = os.path.join(dirNc,'sm_domain',
                                    cardinalPt+'_'+str(bufferNb)+'_'+
                                    str(timestep)+'.nc')

            if(os.path.isfile(filePath)):
                
                plot_pseudocolor(filePath,
                                 fieldPseudoColor,
                                 colorTableName=colorTableName,
                                 fieldMin=fieldMin,
                                 fieldMax=fieldMax,
                                 legend=legend)
                
                plot_contour(filePath,
                             fieldContour,
                             contourColor=(255, 0, 0, 255),
                             contourWidth=contourWidth,
                             contourValue=contourValue,
                             legend=legend)

                nbActivePlots[0]+=2

    if(borders!='None'):
        set_2D_view(borders[0],borders[1],borders[2],borders[3])

    visit.DrawPlots()


    # Second window
    #--------------------
    visit.SetActiveWindow(windows[1])

    # plot the error field for the interior
    filePath = os.path.join(dirNc,'error','error'+str(timestep)+'.nc')
    fieldPseudoColor = 'error_mass'
    colorTableName = 'hot'
    fieldMin = errorMin
    fieldMax = errorMax

    plot_pseudocolor(filePath,
                     fieldPseudoColor,
                     colorTableName=colorTableName,
                     fieldMin=fieldMin,
                     fieldMax=fieldMax,
                     legend=True,
                     scaleType='Log')

    nbActivePlots[1]+=1

    if(borders!='None'):
        set_2D_view(borders[0],borders[1],borders[2],borders[3])

    visit.DrawPlots()

    return nbActivePlots


def clean_DIM_soundwave_with_error(dirNc,
                                   timestep,
                                   nbActivePlots,
                                   maxBfLayers=10,
                                   windows=[1,2]):
    """
    Clean the plots + close the database
    
    dirNc         (String)            : path to the directory where the data files are saved
    timestep      (Integer)           : timestep integer
    nbActivePlots ([Integer,Integer]) : number of active plots per window
    maxBfLayers   (Integer)           : total number of buffer layer per cardinal point
    windows       ([Integer,Integer]) : windows where the plots should be cleaned
    """

    # clean the first window
    #------------------------------
    visit.SetActiveWindow(windows[0])
    # clean the fields
    for i in range(0,nbActivePlots[0]+1):
        visit.DeleteActivePlots()

    # close the database
    filePath = os.path.join(dirNc,'sm_domain','data'+str(timestep)+'.nc')
    visit.OpenDatabase(filePath)

    for cardinalPt in ['N','S','E','W']:

        for bufferNb in range(1,maxBfLayers):

            filePath = os.path.join(dirNc,'sm_domain',
                                    cardinalPt+'_'+str(bufferNb)+'_'+
                                    str(timestep)+'.nc')

            if(os.path.isfile(filePath)):
                visit.OpenDatabase(filePath)


    # clean the second window
    #------------------------------
    visit.SetActiveWindow(windows[1])

    #clean the field
    for i in range(0,nbActivePlots[1]+1):
        visit.DeleteActivePlots()

    #close database
    filePath = os.path.join(dirNc,'error','error'+str(timestep)+'.nc')
    visit.OpenDatabase(filePath)

    return


def get_soundwave_windows(i,nbDirs):
    """
    Return the window IDs for plotting the
    error and the waves

    i      : Integer           : directory ID
    nbDirs : Integer           : total number of directories

    windows: [Integer,Integer] : window ID for plotting the error and the waves
    """

    if(nbDirs==1):
        windows = [1,2]
    else:
        windows = [i+1,nbDirs+i+1]

    return windows


def postprocess_DIM_soundwave_with_error(
    dirsNc,
    timesteps,
    soundMin,
    soundMax,
    errorMin,
    errorMax,
    bubbleValue,
    borders,
    dirOutputPictures='None',
    check_dirOutputPictures_exist=True):
    """
    dirsNc                        ([String,...])              : directories where the Nc files are saved
    timesteps                     ([Integer,Integer,Integer]) : [t_min, t_max, t_step]
    soundMin                      (Float)                     : minimum to see the sound waves
    soundMax                      (Float)                     : maximum to see the sound waves
    errorMin                      (Float)                     : minimum when plotting the error
    errorMax                      (Float)                     : maximum when plotting the error
    bubbleValue                   (Float)                     : value to plot the mass contour of the bubble
    borders                       ([Float,Float,Float,Float]) : [x_min,x_max,y_min,y_max] to set the 2D-view
    dirOutputPictures             (String)                    : path for the output pictures
    check_dirOutputPictures_exist (Boolean)                   : check whether the output folder already exist
    """

    nbDirs=len(dirsNc)

    # create output dir to save the pictures
    if(dirOutputPictures=='None'):
        if(nbDirs==1):
            dirOutputPictures = os.path.join(dirsNc[0],'visit_wave')
        else:
            dirOutputPictures = 'visit_wave'

    if(check_dirOutputPictures_exist):
        if(os.path.isdir(dirOutputPictures)):
            print 'dir for visit picture already exists'
            print dirOutputPictures
            os._exit(1)
            #shutil.rmtree(dirOutputPictures)

        os.mkdir(dirOutputPictures)

    else:
        if(not os.path.isdir(dirOutputPictures)):
            print 'dir for visit picture does not exist'
            print dirOutputPictures
            os._exit(1)
    
    # determine how the windows should be organized
    # we need 2*nb of dirs since: 2 windows per dir:
    # one for error one for one for the density
    
    if(nbDirs>1):
        if(nbDirs<5):
            visit.SetWindowLayout(2*nbDirs)
            windowArray = [nbDirs,2]
            
        else:
            print 'the number of directories is too large'
            print '1<= nbDirs < 5'
            os._exit(1)
    else:
        windowArray = [2,1]
    

    # set the saving options
    set_visit_save_options(outputDir=dirOutputPictures,
                           multipleWindows=True,
                           windowArray=windowArray)


    # check that the folders and files exist
    for dirNc in dirsNc:

        if(not os.path.isdir(os.path.join(dirNc,'sm_domain'))):
            print dirNc+' : sm_domain does not exist'
            os._exit(1)
            
        if(not os.path.isdir(os.path.join(dirNc,'error'))):
            print dirNc+' : error does not exist'
            os._exit(1)
                
        if(not os.path.isfile(os.path.join(dirNc,'sm_domain','data'+str(timesteps[1])+'.nc'))):
            print dirNc+' : data files do not exist'
            os._exit(1)
                    
        if(not os.path.isfile(os.path.join(dirNc,'error','error'+str(timesteps[1])+'.nc'))):
            print dirNc+' : error files do not exist'
            os._exit(1)


    # load the extra variables needed when displaying
    # DIM simulations
    define_dim_variables()

    nbActivePlots = np.zeros([nbDirs,2],dtype=int)

    # create the pictures    
    for timestep in range(timesteps[0],timesteps[1]+1,timesteps[2]):

        # create the plots for the timestep
        for i in range(0,nbDirs):
            
            windows = get_soundwave_windows(i,nbDirs)

            nbActivePlots[i,:] = display_DIM_soundwave_with_error(dirsNc[i],
                                                                  timestep,
                                                                  soundMin,
                                                                  soundMax,
                                                                  errorMin,
                                                                  errorMax,
                                                                  bubbleValue,
                                                                  borders=borders,
                                                                  legend=False,
                                                                  windows=windows)


        # save all the plots as one bui picture
        save_visit_windows()

        # display progress
        print_mg_progress('files processed: '+str(timestep)+'/'+str(timesteps[1]))

        # clean the plots
        for i in range(0,nbDirs):

            windows = get_soundwave_windows(i,nbDirs)

            clean_DIM_soundwave_with_error(dirsNc[i],
                                           timestep,
                                           nbActivePlots[i,:],
                                           windows=windows)


if __name__=='__main__':

    # launch the graphical engine
    launch_visit(window=False)

    # set two windows
    visit.AddWindow()

    # set the options
    dirsNc = []

    for md in [0.05,0.1,0.2,0.3]:

        
        dirNc = os.path.join('/home','jdesmarais',
                             'projects','jcp2015_submission',
                             '20150509_dim2d_bb_trans_cv_r3.5_lin',
                             'dim2d_0.99_0.05_md'+str(md))
        dirsNc.append(dirNc)

        
    timesteps = [600,1017,1]
        
    soundMin =-0.003
    soundMax = 0.003
    errorMin = 1.0e-5
    errorMax = 5.0e-2
    bubbleValue = 1.012
    borders = [-0.5,0.5,-0.5,1.2]

    # postprocess
    postprocess_DIM_soundwave_with_error(
        dirsNc,
        timesteps,
        soundMin,
        soundMax,
        errorMin,
        errorMax,
        bubbleValue,
        borders,
        dirOutputPictures='test',
        check_dirOutputPictures_exist=False)

    print_mg_final('all files processed')
        
    input("PRESS ENTER TO CONTINUE.")
