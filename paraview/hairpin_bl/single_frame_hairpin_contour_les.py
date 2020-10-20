#########################################################################
#
# $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/paraview/single_frame.py $
# $LastChangedDate: 2010-10-08 19:48:34 +0200 (Fri, 08 Oct 2010) $
# $LastChangedBy: ilak@MECH.KTH.SE $
# $LastChangedRevision: 1566 $
#
#########################################################################
#
# Volume plotting, rendering and export for a single file
#
# Milos Ilak, September 16, 2010
# 
#########################################################################

from paraview.simple import *
from paraview import servermanager as sm
from paraviewscripts.read_colormap import colormap

import string
import sys

#########################################################################
#
# Main program
# Arguments are:
# 1. Image file name
# 2. Trajectory data file
# 3. Global parameter file
#
#########################################################################

ifile = sys.argv[1]
tfile = sys.argv[2]
gfile = sys.argv[3]

#ifile = "frame_hairpin_stationary"
#tfile = "tmp_traj2.dat"
#gfile = "glob_stationary.dat"

#
# Read trajectory file
#
g = open(tfile,"rt")
all_lines = g.readlines()
g.close()

#
# Extract cmap color. Not changed for the same dataset.
#
firstline = string.splitfields(all_lines[0])
cmap = firstline[11]

#
# Create reader
#
sline = string.splitfields(all_lines[0])
infile = sline[10]
if infile == "none":
    dnsData = False
else:
    dnsData = True
    reader = LegacyVTKReader(FileNames=infile)
    reader.UpdatePipeline()

#    infile2 = string.replace(infile,".lam.",".u.")
#    reader2 = LegacyVTKReader(FileNames="scaling.u.vtk")
#    reader2.UpdatePipeline()
#    AppendAttributes2 = AppendAttributes( Input=[ reader, reader2 ] )   

#
# Open parameterfile
#
f = open(gfile,"rt")
lines = f.readlines()
f.close()
parallelProjection = int(lines[0].strip("\n"))
sampling = string.splitfields(lines[1].strip("\n"))
sx=int(sampling[0])
sy=int(sampling[1])
sz=int(sampling[2])
cutOff = string.splitfields(lines[2].strip("\n"))
xmin=int(cutOff[0])
xmax=int(cutOff[1])
ymin=int(cutOff[2])
ymax=int(cutOff[3])
zmin=int(cutOff[4])
zmax=int(cutOff[5])
bColor = string.splitfields(lines[3].strip("\n"))
br=int(bColor[0])
bg=int(bColor[1])
bb=int(bColor[2])
vSize = string.splitfields(lines[4].strip("\n"))
vx=int(vSize[0])
vy=int(vSize[1])
bMeshes = string.splitfields(lines[5].strip("\n"))
mColor = string.splitfields(lines[6].strip("\n"))
mr=int(mColor[0])
mg=int(mColor[1])
mb=int(mColor[2])
nstartimage = int(lines[7].strip("\n"))

#
# Get render view and reset camera (crucial!)
#
view=GetRenderView()
ResetCamera()
camera = GetActiveCamera()
Render()

#
# Now extract and show
#
if dnsData:
    myVOI=[xmin,xmax,ymin,ymax,zmin,zmax]
    mysub1 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
    mysub1.VOI=myVOI
    Show(mysub1)

if dnsData:
    #calculator2 = Calculator()
    #calculator2.AttributeMode = 'point_data'
    #calculator2.Function = 'lambda2/scal'
    #calculator2.ResultArrayName = 'lambda2scaled'

    contourFilter = Contour()
    contourFilter.ContourBy = 'lambda2'
    contourFilter.Isosurfaces = [-0.005]

    [colorspace,mytf,lut] = colormap(cmap)
    calculator1 = Calculator(contourFilter)
    calculator1.AttributeMode = 'point_data'
    # calculator1.Function = 'coordsY'
    calculator1.Function = 'vel_u'
    calculator1.ResultArrayName = 'ydistance'
    # a1_ufluct_PVLookupTable = GetLookupTableForArray( "vel_u", 1 )
    a1_ufluct_PVLookupTable = GetLookupTableForArray( "vel_u", 1, Discretize=0, VectorMode='Magnitude', RGBPoints=lut, ColorSpace=colorspace, LockScalarRange=0)

# ****************************
#
# Create contour plot
#
# ****************************
if dnsData:
    dp=GetDisplayProperties(reader)
    dp.Representation = 'Outline'
    dp.Visibility = 0

    dp=GetDisplayProperties(mysub1)
    dp.Representation = 'Outline'
    dp.Visibility = 0

    dp=GetDisplayProperties(contourFilter)
    dp.Representation = 'Surface'
    dp.Visibility = 0

    dp=GetDisplayProperties(calculator1)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'ydistance'
    dp.LookupTable = a1_ufluct_PVLookupTable
    dp.ColorAttributeType = 'POINT_DATA'

#
# Load grid
#

#
# Loop over all camera positions and generate corresponding images
#
for line in all_lines:
     
    sline = string.splitfields(line)
    cpos1 = float(sline[0])
    cpos2 = float(sline[1])
    cpos3 = float(sline[2])
    fpos1 = float(sline[3])
    fpos2 = float(sline[4])
    fpos3 = float(sline[5])
    zoom  = float(sline[6])
    vpos1 = float(sline[7])
    vpos2 = float(sline[8])
    vpos3 = float(sline[9])
    infile = sline[10]
#    cmap = sline[11]
    iline = int(sline[12])
    outfile = ifile + "_%04d"%(iline+nstartimage) + ".png"

    print "Input file : ", infile
    print "Output file : ", outfile
    print "Parameter file : ", gfile

# make sure that you start from exactly 1902 for this to work :-)


    if int(iline+nstartimage) > 1901:
        if (iline+nstartimage) <= 1960:
            myop = None
            myop = float(1.0-float(iline-1901)/60.0)

    if int(iline+nstartimage) > 1960:
        myop = None
        myop = 0.0

    if int(iline+nstartimage) <= 1901:
        myop = None
        myop = float(1.0)

    print "Opacity ",myop

    if bMeshes[0] != "none":
        for bMesh in bMeshes:
            reader=LegacyVTKReader(FileNames=bMesh)
            reader.UpdatePipeline()

            dp = GetDisplayProperties(reader)
            dp.Representation = 'Wireframe'
            dp.AmbientColor = [mr/255.0, mg/255.0, mb/255.0]
            dp.Opacity = float(myop)

    myop = None
    #
    # Need to sync these with the camera position
    #
    view.Background = [br, bg, bb]
    view.CameraFocalPoint = [fpos1, fpos2, fpos3]
    view.CameraPosition = [cpos1, cpos2, cpos3]
    view.ViewSize = [vx, vy] 

    #
    # Light settings
    #
    view.UseLight = 1
    view.LightSwitch = 0
    view.KeyLightIntensity =2.0 # Prev 1.5 Def 0.75 2.0
    view.FillLightWarmth = 0.4   # Prev 0.4 Def 0.4  1.0

    #
    # Turn on zoom
    #
    view.CameraViewUp = [vpos1, vpos2, vpos3]
    if parallelProjection:
        view.CameraParallelProjection = 1
        view.CameraParallelScale = float(zoom)
        # camera.ComputeViewPlaneNormal() # Alternatively to the CameraViewUp
        # view.CameraParallelScale = 30.0
        # view.CameraViewAngle = 30.0
    else:
        view.CameraParallelProjection = 0
        camera.Dolly(zoom)

    #
    # Make a screen dump
    #
    # view.UseOffscreenRenderingForScreenshots = 0
    view.WriteImage(outfile, "vtkPNGWriter",1)
