# File:       vtklines_to_curves.py
# Programmer: Jeremy Meredith
# Date:       August 26, 2011
#
# Description:
#   Input a file containing a VTK_POLYDATA with only 2-point line segments
#   and either a point or cell scalar.  Coalesce the short line segments into
#   polylines.  Output each polyline as a function of scalar-value vs distance.
#
# Usage:
#   python vtklines_to_curves.py <inputfile.vtk> <outputfile.curve>
#
import sys, os, math
import matplotlib.pyplot as plt 
 
# ----------------------------------------------------------------------------
#  Parse the points, line segments and cell or point scalars from a VTK file
# ----------------------------------------------------------------------------
def ParseFile(fn):
    f = open(fn)
    try:
        npts = 0
        pts = []
        nlines = 0
        lines = []
        vals  = []
        name  = "values"
        ispoint = -1
        # modes: 0 = searching for POINTS
        #        1 = reading points
        #        2 = looking for LINES
        #        3 = reading lines
        #        4 = looking for CELL_DATA, POINT_DATA, or SCALAR
        #        5 = found scalar, skipping LOOKUP_TABLE
        #        6 = reading scalar values
        mode = 0
        for l in f:
            words = l.split()
            if len(words) == 0:
                None
            elif mode == 0:
                if words[0] == "POINTS":
                    npts = int(words[1])
                    if npts == 0:
                        raise "Couldn't extract number of points"
                    mode = 1
            elif mode == 1:
                for w in words:
                    pts.append(float(w))
                if len(pts) >= npts*3:
                    mode = 2
            elif mode == 2:
                if words[0] == "LINES":
                    nlines = int(words[1])
                    if nlines == 0:
                        raise "Couldn't extract number of lines"
                    nints  =  int(words[2])
                    if nints != 3*nlines:
                        raise "Expected only 2-point line segments"
                    mode = 3
            elif mode == 3:
                if int(words[0]) != 2:
                    raise "Expected only 2-point line segments"
                seg = [int(words[1]),int(words[2])]
                lines.append(seg)
                if len(lines) >= nlines:
                    mode = 4
            elif mode == 4:
                if words[0] == "CELL_DATA":
                    ispoint = 0
                elif words[0] == "POINT_DATA":
                    ispoint = 1
                elif words[0] == "SCALARS":
                    if ispoint < 0:
                        raise "Found scalars before CELL_DATA or POINT_DATA"
                    name = words[1]
                    mode = 5
            elif mode == 5:
                if words[0] == "LOOKUP_TABLE":
                    mode = 6
            elif mode == 6:
                for w in words:
                    vals.append(float(w))
                if ((ispoint==1 and len(vals) >= npts) or
                    (ispoint==0 and len(vals) >= nlines)):
                    mode = 7
            else: # mode == 7
                break
    finally:
        f.close()
    return (name, ispoint, pts, lines, vals)
 
# ----------------------------------------------------------------------------
#  Operations on points
# ----------------------------------------------------------------------------
def calcdist(pt0, pt1):
    dx = pt0[0] - pt1[0]
    dy = pt0[1] - pt1[1]
    dz = pt0[2] - pt1[2]
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    return dist
 
def calcmid(pt0, pt1):
    x = (pt0[0] + pt1[0]) / 2.
    y = (pt0[1] + pt1[1]) / 2.
    z = (pt0[2] + pt1[2]) / 2.
    return (x,y,z)
 
def getpoint(pts, i):
    return (pts[3*i + 0], pts[3*i + 1], pts[3*i + 2])
 

#-----------------------------------------------------------------------------
# Operations generating the curve as a polyline
#-----------------------------------------------------------------------------
def get_polyline_from_vtk(infile):
    '''
    @description
    Extract all the segments constituing the curve from
    the vtk file. The segments are them assembled to create
    the polyline
    '''

    # parse the file
    (name, ispoint, pts, lines, vals) = ParseFile(infile)
	 
    # for cell vars, put the point vals in a dictionary
    linevals = {}
    if not ispoint:
        for i in range(len(lines)):
            seg = lines[i]
            linevals[seg[0],seg[1]] = vals[i]
            linevals[seg[1],seg[0]] = vals[i]
            
    # group the line segments into polylines, first pass
    polylines = []
    for seg in lines:
        found = False
        for i in range(len(polylines)):
            # if we can stick it on either end of a polyline, do it
            if seg[0] == polylines[i][0]:
                polylines[i].insert(0,seg[1])
                found = True
            elif seg[1] == polylines[i][0]:
                polylines[i].insert(0,seg[0])
                found = True
            elif seg[0] == polylines[i][-1]:
                polylines[i].append(seg[1])
                found = True
            elif seg[1] == polylines[i][-1]:
                polylines[i].append(seg[0])
                found = True
            if found:
                break
        if found == False:
            polylines.append(seg)
	 
    # now group up the polylines until we can't do it anymore
    done = False
    while not done:
        done = True
        for i in range(len(polylines)):
            for j in range(i):
                polyi = polylines[i]
                polyj = polylines[j]
                # if either end of the polylines match, merge them
                if polyi[0] == polyj[0]:
                    polylines.remove(polyj)
                    polyi.reverse()
                    polyi.extend(polyj[1:])
                    done = False
                elif polyi[-1] == polyj[0]:
                    polylines.remove(polyj)
                    polyi.extend(polyj[1:])
                    done = False
                elif polyi[0] == polyj[-1]:
                    polylines.remove(polyi)
                    polyj.extend(polyi[1:])
                    done = False
                elif polyi[-1] == polyj[-1]:
                    polylines.remove(polyj)
                    polyj.reverse()
                    polyi.extend(polyj[1:])
                    done = False
	 
                if done == False:
                    break
            if done == False:
                break

    return [ispoint,polylines]


def get_curve_coordinates(ispoint, polylines, outfile=''):
    '''
    @description
    Extract the position of the points in the polyline to
    create the coordinate data of the curve:
     - data[0]: x-coordinates
     - data[1]: y-coordinates
    These data can optionally be written on an output file
    '''

    writeOutput=(outfile!='')

    if(writeOutput):
        out = open(outfile, 'w')

    data=[]
    data.append([])
    data.append([])

    for p in range(len(polylines)):
        polyline = polylines[p]

        if(writeOutput):
            out.write("# %s%04d\n" % (name, p))

        # note it's different between cell- and point-centered values
        if ispoint:
            pos = 0.0
            for i in range(len(polyline)-1):
                pt = getpoint(pts, polyline[i])

                if(writeOutput):
                    out.write("%f %f\n" % (pt[0], pt[1]))

                data[0].append(pt[0])
                data[1].append(pt[1])

        else:
            pos = 0.0
            lastpt = (0,0,0)
            for i in range(len(polyline)-1):
                pt0 = getpoint(pts, polyline[i])
                pt1 = getpoint(pts, polyline[i+1])
                pt = calcmid(pt0,pt1)

                if(writeOutput):
                    out.write("%f %f\n" % (pt[0], pt[1]))

                data[0].append(pt[0])
                data[1].append(pt[1])

    if(writeOutput):
        out.close()

    return data

 
# ----------------------------------------------------------------------------
#  Main routine
# ----------------------------------------------------------------------------
if __name__=="__main__":

    # get command-line arguments
    if len(sys.argv) != 3:
        print "Usage:",sys.argv[0],"<inputfile.vtk> <outputfile.curve>"
        sys.exit(1)
        
    infile = sys.argv[1]
    if not os.path.exists(infile):
        print "Error:",infile,"didn't exist."
        sys.exit(1)


    # extract the polyline form the vtk file
    [ispoint,polylines] = get_polyline_from_vtk(infile)


    # extract the coordinate data and write output file
    outfile = sys.argv[2]
    graph_data = get_curve_coordinates(ispoint,
                                       polylines,
                                       outfile=outfile)    

    # plot the curve
    width = 3

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8,4))
    ax = plt.subplot(111)
    
    plt.plot(
        graph_data[0],
        graph_data[1],
        '-',
        linewidth=width,
        color='black')
    
    plt.show()
