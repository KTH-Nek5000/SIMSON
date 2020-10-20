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
# Contour plots of velocity, vorticity and wavemaker for the jet
#
# Milos Ilak, October 25, 2010
# 
#########################################################################

from paraview.simple import *
from paraview import servermanager as sm
from paraviewscripts import *

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

options = sm.vtkPVOptions()
options.SetUseOffscreenRendering(1)
#
# Create reader
#
if infile == "none":
    dnsData = False
else:
    dnsData = True
    reader=LegacyVTKReader(FileNames=infile)
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

    mysub2 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
    mysub2.VOI=[0, 816, 20, 121, 0, 287]

#
# Special case for lambda2 here
#
    mysub3 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
    mysub3.VOI=[300, 816, 0, 121, 0, 287]
    
#
# Get render view and reset camera (crucial!)
#
view=GetRenderView()
view.ViewSize = [vx, vy] 
ResetCamera()

# ****************************
#
# Mess with the camera
#
# ****************************

#
# Need to sync these with the camera position
#
view.Background = [br, bg, bb]
view.CameraFocalPoint = [fpos1, fpos2, fpos3]
view.CameraPosition = [cpos1, cpos2, cpos3]
view.UseLight = 1
view.CameraViewAngle = 3.0

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
# Spanwise vorticity volume plot
#

SetActiveSource(reader)
Contour1 = Contour( PointMergeMethod="Uniform Binning" )
Contour1.PointMergeMethod = "Uniform Binning"
Contour1.ContourBy = ['POINTS', 'vel_u']
Contour1.Isosurfaces = [-0.0001]
dpw = GetDisplayProperties(Contour1)
dpw.DiffuseColor = [0.0, 0.33, 1.0]
dpw.Opacity = float(1.0)

## #dp3=GetDisplayProperties(mysub3)
SetActiveSource(mysub1)
Contour2 = Contour( PointMergeMethod="Uniform Binning" )
Contour2.PointMergeMethod = "Uniform Binning"
Contour2.ContourBy = ['POINTS', 'vel_v']
Contour2.Isosurfaces = [0.477]


DataRepresentation13 = Show()
dp33 = GetDisplayProperties(Contour2)
DataRepresentation13.Texture = []
dp33.EdgeColor = [1.0, 1.0, 1.0] 
dp33.DiffuseColor = [0.6, 0.0 ,0.0]
dp33.Opacity = 1.0 

SetActiveSource(mysub1)
Contour3 = Contour( PointMergeMethod="Uniform Binning" )
Contour3.PointMergeMethod = "Uniform Binning"
Contour3.ContourBy = ['POINTS', 'vort_z']
Contour3.Isosurfaces = [-0.9]

DataRepresentation13 = Show()
dp33 = GetDisplayProperties(Contour3)
DataRepresentation13.Texture = []

dp33.EdgeColor = [1.0, 1.0, 1.0]
dp33.DiffuseColor = [0.0, 0.33, 0.0]
dp33.Opacity = 0.3 

#
# Add contour plot of wavemaker
#
wavemaker= LegacyVTKReader( FileNames=['/home/x_milil/nb2/low_r/modesr0.7/r0.7overlap.u.vtk'] )

Contour7 = Contour( PointMergeMethod="Uniform Binning" )
Contour7.PointMergeMethod = "Uniform Binning"
Contour7.ContourBy = ['POINTS', 'vel_u']
Contour7.Isosurfaces = [31215470430.0]

dpw = GetDisplayProperties(Contour7)
dpw.DiffuseColor = [1.0, 1.0, 1.0]
dpw.Opacity = 1.0 #0.59999999999999998
#
# Load grid
#
if bMeshes[0] != "none":
    for bMesh in bMeshes:
        reader=LegacyVTKReader(FileNames=bMesh)
        reader.UpdatePipeline()

        dp=GetDisplayProperties()
        dp.Representation = 'Wireframe'
#
# Add background plane
#
Plane1=Plane()
Plane1.Origin = [0.0, 0.0, -15.0]
Plane1.Point1 = [75.0, 0.0, -15.0]
Plane1.Point2 = [0.0, 0.0, 15.0]

dpplane = GetDisplayProperties(Plane1)
dpplane.DiffuseColor = [0.35, 0.35, 0.35]
dpplane.Opacity = 0.5
#
# Add jet orifice
#
Disk1=Disk()
Disk1.OuterRadius = 1.5
Disk1.InnerRadius = 0.0
Disk1.CircumferentialResolution = 100
Disk1.RadialResolution = 100

dpdisk = GetDisplayProperties(Disk1)
dpdisk.Position = [9.375, 0.0, 0.0]
dpdisk.Orientation = [90.0, 0.0, 0.0]
dpdisk.DiffuseColor = [0.0, 0.0, 0.0]
dpdisk.EdgeColor = [0.0, 0.0, 0.0]
dpdisk.Representation = 'Surface With Edges'
#
# Draw image and render it
#
Render()
    
#
# Make a screen dump
#
WriteImage(outfile)
