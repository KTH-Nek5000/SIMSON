#!/software/apps/paraview/3.8.1/bin/pvpython
#########################################################################
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
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
from read_colormap import colormap

import string
import sys

#########################################################################
#
# Main program
#
#########################################################################

#
# Read in arguments
#
infile = sys.argv[1]
outfile = sys.argv[2]
cmap = sys.argv[3]
cpos1 = float(sys.argv[4])
cpos2 = float(sys.argv[5])
cpos3 = float(sys.argv[6])
fpos1 = float(sys.argv[7])
fpos2 = float(sys.argv[8])
fpos3 = float(sys.argv[9])
zoom = float(sys.argv[10])
vpos1 = float(sys.argv[11])
vpos2 = float(sys.argv[12])
vpos3 = float(sys.argv[13])
gfile = sys.argv[14]

print "Input file : ", infile
print "Output file : ", outfile
print "Parameter file : ", gfile

#
# Create reader
#
if infile == "none":
    dnsData = False
else:
    dnsData = True
    reader = LegacyVTKReader(FileNames=infile)
    reader.UpdatePipeline()

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

#
# Now extract and show
#
if dnsData:
    myVOI=[xmin,xmax,ymin,ymax,zmin,zmax]
    mysub1 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
    mysub1.VOI=myVOI
    Show(mysub1)


# ****************************
#
# Set up the transfer function (an ugly hack for now)
#
# ****************************

#
# Start by getting the range of the data
#
#if dnsData:
#    min=mysub1.PointData[0].GetRange()[0]
#    max=mysub1.PointData[0].GetRange()[1]

#
# Read in colormap and transfer function
#
if dnsData:
    [colorspace,mytf,lut] = colormap(cmap)
    a1_lambda2_PiecewiseFunction = CreatePiecewiseFunction(Points=mytf)
    a1_lambda2_PVLookupTable = GetLookupTableForArray( "lambda2", 1, Discretize=0, VectorMode='Magnitude', RGBPoints=lut, ColorSpace=colorspace, LockScalarRange=0)


# ****************************
#
# Now create volume plot
#
# ****************************
if dnsData:
    dp=GetDisplayProperties(mysub1)
    dp.Representation = 'Volume'
    dp.ScalarOpacityFunction = a1_lambda2_PiecewiseFunction
    dp.ColorArrayName = 'lambda2'
    dp.LookupTable = a1_lambda2_PVLookupTable
    dp.ColorAttributeType = 'POINT_DATA'


# ****************************
#
# Mess with the camera
#
# ****************************

#
# Get render view and reset camera (crucial!)
#
view=GetRenderView()
ResetCamera()

#
# Need to sync these with the camera position
#
view.Background = [br, bg, bb]
view.CameraFocalPoint = [fpos1, fpos2, fpos3]
view.CameraPosition = [cpos1, cpos2, cpos3]
view.ViewSize = [vx, vy] 

#
# Turn on zoom
#
camera = GetActiveCamera()
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
# Load grid
#
if bMeshes[0] != "none":
    for bMesh in bMeshes:
        reader=LegacyVTKReader(FileNames=bMesh)
        reader.UpdatePipeline()

        dp=GetDisplayProperties()
        dp.Representation = 'Wireframe'
        dp.AmbientColor = [mr/255.0, mg/255.0, mb/255.0]

#
# Draw image and render it
#
Show()
Render()
    
#
# Make a screen dump
#
WriteImage(outfile)
