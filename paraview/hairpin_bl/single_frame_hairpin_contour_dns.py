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

#
# Read parts, if part is inluded in file name read two files 
#

infile = "none"

if infile == "none":
    dnsData = False
else:
    dnsData = True
    

    print "reading field 1"
    infile1 = "l2s_merged_07.vtk"
    # infile1 = "t150931_00001-01000.lam.vtk"
    reader1 = LegacyVTKReader(FileNames=infile1)
    reader1.UpdatePipeline()

    print "reading field 2"
    infile2 = "u_merged_07.vtk"
    # infile2 = "t150931_00001-01000.u.vtk"
    reader2 = LegacyVTKReader(FileNames=infile2)
    reader2.UpdatePipeline()
    AppendAttributes1 = AppendAttributes( Input=[ reader1, reader2 ] )

    print "reading field 3"
    infile3 = "l2s_merged_714.vtk"
    # infile3 = "t150931_07001-08000.lam.vtk"
    reader3 = LegacyVTKReader(FileNames=infile3)
    reader3.UpdatePipeline()

    print "reading field 4"
    infile4 = "u_merged_714.vtk"
    # infile4 = "t150931_07001-08000.u.vtk"
    reader4 = LegacyVTKReader(FileNames=infile4)
    reader4.UpdatePipeline()
    AppendAttributes2 = AppendAttributes( Input=[ reader3, reader4 ] )

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
#view.OrientationAxesVisibility = 0
#view.CenterAxesVisibility = 0
ResetCamera()
camera = GetActiveCamera()
Render()

#
# Now extract and show
#
#if dnsData:
#    myVOI=[xmin,xmax,ymin,ymax,zmin,zmax]
#    mysub1 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
#    mysub1.VOI=myVOI
#    Show(mysub1)

#    SetActiveSource(reader2)
#    mysub2 = ExtractSubset( IncludeBoundary=0, SampleRateK=sz, SampleRateJ=sy, SampleRateI=sx )
#    mysub2.VOI=myVOI
#    Show(mysub2)


if dnsData:

    [colorspace,mytf,lut] = colormap(cmap)

    #    contourFilter1 = Contour(mysub1)
    #    contourFilter1 = Contour(reader1)
    contourFilter1 = Contour(AppendAttributes1)
    contourFilter1.ContourBy = 'lambd_2'
    contourFilter1.Isosurfaces = [-0.005]

    calculator1 = Calculator(contourFilter1)
    calculator1.AttributeMode = 'point_data'
    calculator1.Function = 'vel_u'
    calculator1.ResultArrayName = 'ydistance'
    a1_ufluct_PVLookupTable = GetLookupTableForArray('vel_u', 1, Discretize=0, VectorMode='Magnitude', RGBPoints=lut, ColorSpace=colorspace, LockScalarRange=0)

    #    contourFilter2 = Contour(reader3)
    contourFilter2 = Contour(AppendAttributes2)
    contourFilter2.ContourBy = 'lambd_2'
    contourFilter2.Isosurfaces = [-0.005]

    calculator2 = Calculator(contourFilter2)
    calculator2.AttributeMode = 'point_data'
    calculator2.Function = 'vel_u'
    calculator2.ResultArrayName = 'ydistance'
    a1_ufluct_PVLookupTable = GetLookupTableForArray('vel_u', 1, Discretize=0, VectorMode='Magnitude', RGBPoints=lut, ColorSpace=colorspace, LockScalarRange=0)

# ****************************
#
# Create contour plot
#
# ****************************
if dnsData:
    dp=GetDisplayProperties(reader1)
    dp.Representation = 'Outline'
    dp.Visibility = 0

    dp=GetDisplayProperties(reader2)
    dp.Representation = 'Outline'
    dp.Visibility = 0

    dp=GetDisplayProperties(reader3)
    dp.Representation = 'Outline'
    dp.Visibility = 0

    dp=GetDisplayProperties(reader4)
    dp.Representation = 'Outline'
    dp.Visibility = 0


#    dp=GetDisplayProperties(mysub1)
#    dp.Representation = 'Outline'
#    dp.Visibility = 0

#    dp=GetDisplayProperties(mysub2)
#    dp.Representation = 'Outline'
#    dp.Visibility = 0

    dp=GetDisplayProperties(contourFilter1)
    dp.Representation = 'Surface'
    dp.Visibility = 0

    dp=GetDisplayProperties(contourFilter2)
    dp.Representation = 'Surface'
    dp.Visibility = 0
#    dp.LookupTable = a1_ufluct_PVLookupTable
#    dp.ColorAttributeType = 'POINT_DATA'

    dp=GetDisplayProperties(calculator1)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'ydistance'
    dp.LookupTable = a1_ufluct_PVLookupTable
    dp.ColorAttributeType = 'POINT_DATA'

    dp=GetDisplayProperties(calculator2)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'ydistance'
    dp.LookupTable = a1_ufluct_PVLookupTable
    dp.ColorAttributeType = 'POINT_DATA'

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

    #
    # Need to sync these with the camera position
    #
    view.Background = [br, bg, bb]
    view.CameraFocalPoint = [fpos1, fpos2, fpos3]
    view.CameraPosition = [cpos1, cpos2, cpos3]
    view.ViewSize = [vx, vy] 
    view.UseLight = 1
    view.LightSwitch =0
    view.KeyLightIntensity = 2.0 # Prev 1.5 Def 0.75
    view.FillLightWarmth = 0.4   # Prev 0.4 Def 0.4

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
    # Draw image and render it
    #
    # Show()
    # Render()

    #
    # Make a screen dump
    #
    # view.UseOffscreenRenderingForScreenshots = 0
    # WriteImage(outfile)
    view.WriteImage(outfile, "vtkPNGWriter",1)
