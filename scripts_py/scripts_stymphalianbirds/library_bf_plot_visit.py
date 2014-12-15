import sys
import os
sys.path.insert(0,os.environ['VISIT_PYTHON_LIB'])
import visit
import fnmatch
import shutil
import subprocess


#list the files corresponding to the detectors of
#the simulation in the directory cdir at the time step
# given by the user
def get_dct_files(cdir,timestep):

    dct_file=[]

    for file in os.listdir(cdir):

        if file.endswith("detectors"+str(timestep)+".curve"):
            dct_file.append(file)

    print dct_file

    return dct_file


#list the files corresponding to the simulation
#in the directory cdir at the time step given by
#the user
def get_bf_files(cdir,timestep):

    bf_file=[]

    for file in os.listdir(cdir):

        if file.endswith("_"+str(timestep)+".nc"):
            bf_file.append(file)

    print bf_file

    return bf_file


#get the file corresponding to the simulation
#for the interior at the time step given by
#the user
def get_interior_file(cdir,timestep):

    data_file=[]

    for file in os.listdir(cdir):

        if file.endswith("data"+str(timestep)+".nc"):
            data_file.append(file)

    print data_file

    return data_file


#get the total number of timesteps written
def get_total_nb_timesteps(cdir):

    nbTimesteps = 0

    for file in os.listdir(cdir):
        if fnmatch.fnmatch(file, 'data[0-9]*.nc'):
            
            timestep = int(file.split('data')[1].split('.nc')[0])
            if(timestep>nbTimesteps):
                nbTimesteps = timestep

    return nbTimesteps
    

#open and plot the data "field" between min and max
#using the colorTable
def open_and_plot_pseudocolor(
    filePath,field,
    fieldMin,
    fieldMax,
    colorTableName,
    legendFlag):

    #open file
    visit.OpenDatabase(filePath, 0)

    #add the field to the plot
    visit.AddPlot("Pseudocolor", field, 1, 1)


    #configuration of the plot
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = fieldMin
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = fieldMax
    PseudocolorAtts.centering = PseudocolorAtts.Zonal  # Natural, Nodal, Zonal

    PseudocolorAtts.colorTableName = colorTableName
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
    PseudocolorAtts.legendFlag = legendFlag
    PseudocolorAtts.lightingFlag = 1

    visit.SetPlotOptions(PseudocolorAtts)


#set the 2D view: min and max
def set_2D_view(Xmin,Xmax,Ymin,Ymax):

    View2DAtts = visit.View2DAttributes()
    View2DAtts.windowCoords = (Xmin, Xmax, Ymin, Ymax)
    View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 1
    visit.SetView2D(View2DAtts)


#save view as png picture
def save_view_as_png(outputDirectory,fileName):

    SaveWindowAtts = visit.SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = outputDirectory
    SaveWindowAtts.fileName = fileName
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0

    visit.SetSaveWindowAttributes(SaveWindowAtts)

    visit.SaveWindow()


#open and plot the grdpts_id for a buffer layer file
def open_and_plot_pseudocolor_grdpts_id(filePath,legendFlag):
    open_and_plot_pseudocolor(filePath,'grdpts_id',0,3,"bf_layer",legendFlag)


#set the detector curve attributes
def set_dct_CurveAtts(detector_color=(255, 255, 0, 255)):

    CurveAtts = visit.CurveAttributes()
    CurveAtts.showLines = 1
    CurveAtts.lineStyle = CurveAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    CurveAtts.lineWidth = 0
    CurveAtts.showPoints = 1
    CurveAtts.symbol = CurveAtts.Point  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
    CurveAtts.pointSize = 5
    CurveAtts.pointFillMode = CurveAtts.Static  # Static, Dynamic
    CurveAtts.pointStride = 1
    CurveAtts.symbolDensity = 50
    CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
    CurveAtts.curveColor = detector_color
    CurveAtts.showLegend = 0
    CurveAtts.showLabels = 0
    CurveAtts.designator = ""
    CurveAtts.doBallTimeCue = 0
    CurveAtts.ballTimeCueColor = (0, 0, 0, 255)
    CurveAtts.timeCueBallSize = 0.01
    CurveAtts.doLineTimeCue = 0
    CurveAtts.lineTimeCueColor = (0, 0, 0, 255)
    CurveAtts.lineTimeCueWidth = 0
    CurveAtts.doCropTimeCue = 0
    CurveAtts.timeForTimeCue = 0
    CurveAtts.fillMode = CurveAtts.NoFill  # NoFill, Solid, HorizontalGradient, VerticalGradient
    CurveAtts.fillColor1 = (255, 0, 0, 255)
    CurveAtts.fillColor2 = (255, 100, 100, 255)
    CurveAtts.polarToCartesian = 0
    CurveAtts.polarCoordinateOrder = CurveAtts.R_Theta  # R_Theta, Theta_R
    CurveAtts.angleUnits = CurveAtts.Radians  # Radians, Degrees
    visit.SetPlotOptions(CurveAtts)


#open and plot a .curve file
def open_and_plot_curve(filePath,
                        N_detector_color=(255, 255, 0, 255),
                        S_detector_color=(255, 0, 255, 255),
                        W_detector_color=(0, 255, 0, 255),
                        E_detector_color=(0, 255, 255, 255)):

    visit.OpenDatabase(filePath, 0)

    visit.AddPlot("Curve", "N_detectors", 1, 1)
    set_dct_CurveAtts(detector_color=N_detector_colors)

    visit.AddPlot("Curve", "S_detectors", 1, 1)
    set_dct_CurveAtts(detector_color=S_detector_colors)

    visit.AddPlot("Curve", "W_detectors", 1, 1)
    set_dct_CurveAtts(detector_color=W_detector_colors)

    visit.AddPlot("Curve", "E_detectors", 1, 1)
    set_dct_CurveAtts(detector_color=E_detector_colors)


#plot all the buffer layer and the interior grdpts_id at a defined timestep
def plot_and_print_grdpts_id(
    dirInputNcFiles,
    dirOutputPictures,
    timestep):

    #get the buffer layer netcdf files
    bf_files = get_bf_files(dirInputNcFiles,timestep)


    #plot the interior grdpts_id
    open_and_plot_pseudocolor_grdpts_id(
        os.path.join(dirInputNcFiles,'interior_grdpts_id.nc'),
        1)

    #plot each buffer layer
    for bf_file in bf_files:
        open_and_plot_pseudocolor_grdpts_id(
            os.path.join(dirInputNcFiles,bf_file),
            0)

    #draw the plots
    visit.AddOperator("DualMesh", 1)
    visit.DrawPlots()

    #print the view in a png file
    save_view_as_png(".",os.path.join(dirOutputPictures,'grdpts_id'))

    #delete all the active plots
    for i in range(0,len(bf_files)+1):
        visit.DeleteActivePlots()


#plot all the buffer layer and the interior fields at a defined timestep
def plot_and_print_field(
    dirInputNcFiles,
    dirOutputPictures,
    timestep,
    field,
    fieldMin,
    fieldMax,
    plotDetectors=False):

    #get the buffer layer netcdf files
    bf_files = get_bf_files(dirInputNcFiles,timestep)

    #get the interior domain netcdf file
    interior_files = get_interior_file(dirInputNcFiles,timestep)

    #get the detectors file
    dct_files = get_dct_files(dirInputNcFiles,timestep)

    #plot the interior
    for interior_file in interior_files:
        open_and_plot_pseudocolor(
            os.path.join(dirInputNcFiles,interior_file),
            field,
            fieldMin,
            fieldMax,
            "hot",
            1)

    #plot the buffer layers
    for bf_file in bf_files:
        open_and_plot_pseudocolor(
            os.path.join(dirInputNcFiles,bf_file),
            field,
            fieldMin,
            fieldMax,
            "hot",
            0)

    #plot the detectors
    if(plotDetectors):
        for dct_file in dct_files:
            open_and_plot_curve(dct_file)

    #draw the plots
    #visit.AddOperator("DualMesh", 1)
    visit.DrawPlots()

    #print the view in a png file
    save_view_as_png(".",os.path.join(dirOutputPictures,field))

    #delete all the active plots
    for i in range(0,len(bf_files)+len(interior_files)+len(dct_files)):
        visit.DeleteActivePlots()

  
#print the buffer layer and the interior grdpts_id and
#create pictures for each time step
def create_grdpts_id_pictures(
    dirInputNcFiles,
    xmin,
    xmax,
    ymin,
    ymax,
    step=1):


    #determine the total number of timesteps
    nbTimesteps = get_total_nb_timesteps(dirInputNcFiles)

    #create a folder where the pictures will be saved
    dirOutputPictures=dirInputNcFiles+'/visit_grdpts_id'

    if(os.path.isdir(dirOutputPictures)):
        shutil.rmtree(dirOutputPictures)

    os.mkdir(dirOutputPictures)


    #set the 2D view: [xmin,xmax]x[ymin,ymax]
    set_2D_view(xmin,xmax,ymin,ymax)


    #loop over the time steps, extract the grdpts_id
    #and create a picture from the view
    for timestep in range(0,nbTimesteps+1,step):

        plot_and_print_grdpts_id(
            dirInputNcFiles,
            dirOutputPictures,
            timestep)


#print the interior and the buffer layer field and create
#pictures for each time step
def create_field_pictures(
    dirInputNcFiles,
    field,
    fieldMin,
    fieldMax,
    xmin,
    xmax,
    ymin,
    ymax,
    step=1):


    #determine the total number of timesteps
    nbTimesteps = get_total_nb_timesteps(dirInputNcFiles)

    #create a folder where the pictures will be saved
    dirOutputPictures=dirInputNcFiles+'/visit_'+field

    if(os.path.isdir(dirOutputPictures)):
        shutil.rmtree(dirOutputPictures)

    os.mkdir(dirOutputPictures)


    #set the 2D view: [xmin,xmax]x[ymin,ymax]
    set_2D_view(xmin,xmax,ymin,ymax)


    #loop over the time steps, extract the grdpts_id
    #and create a picture from the view
    for timestep in range(0,nbTimesteps+1,step):

        plot_and_print_field(
            dirInputNcFiles,
            dirOutputPictures,
            timestep,
            field,
            fieldMin,
            fieldMax)


#create a movie out of the pictures
def create_movie(
    dirInputFiles,
    filePattern,
    movieFileName,
    movieRate=1):

    if(os.path.isfile(os.path.join(dirInputFiles,movieFileName))):
        os.remove(os.path.join(dirInputFiles,movieFileName))

    cmd="ffmpeg"
    cmd+=" -r "+str(movieRate)
    cmd+=" -i "+os.path.join(dirInputFiles,filePattern)
    cmd+=" "+movieFileName
    cmd+=" && mv "+movieFileName+" "+dirInputFiles
    subprocess.call(cmd, shell=True)
    
