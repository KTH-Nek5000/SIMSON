#!/usr/bin/env python
############################################################################
#
# $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/test/fltype_-2/test.py $
# $LastChangedDate: 2010-09-03 15:20:32 +0200 (Fri, 03 Sep 2010) $
# $LastChangedBy: mattias@MECH.KTH.SE $
# $LastChangedRevision: 1507 $
#
############################################################################
#
# Merge multiple vtk datafiles
# See option --help for more info
#
# Written by Mattias Chevalier 2011
#
############################################################################

import sys
import os
import string
#import hashlib
#import shutil
import glob
import time
from array import array
import numpy
import base64

############################################################################
# Functions
############################################################################


#
# Add message about time date revision compiler and who compiled it to log file
#
# fname   Filename to be read
# tdata   Type of data (lambda2, vel_u ,vel_v, vel_w, omega_x, omega_y, omega_z, ufluct )
#def readBinaryVTK(fname,tdata,x,y,z,data):
#def readBinaryVTK(fname,tdata,xcoord,ycoord,zcoord):
def readBinaryVTK(fname,tdata,ifile):

    global xcoord,ycoord,zcoord,nx,ny,nz,fdatad

    # Use dir name 
    st  = "--- Read binary VTK file -------------------------------------------" + "\n"
    st += "Filename                      ::: " + fname + "\n"
    st += "Type of data to be retrieved  ::: " + tdata + "\n"
    st += "--------------------------------------------------------------------" + "\n"

    try:
        g = open(fname,"rb")

        # Read vtk header
        i = 0
        dbyte = 0
        vtkHeader = ""
        while dbyte != chr(10):
            i = i + 1
            dbyte = g.read(1)
#            print "Byte:",i," --- ",dbyte
#            if i > 60:
#                break
            vtkHeader += dbyte
        print vtkHeader
        
        # Read vtk title
        dbyte = 0
        vtkTitle = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkTitle += dbyte
        print vtkTitle

        # Read vtk title
        dbyte = 0
        vtkBinary = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkBinary += dbyte
        print vtkBinary

        # Read vtk grid type (RECTILINEAR_GRID, STRUCTURED_POINTS)
        dbyte = 0
        vtkGridType = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkGridType += dbyte
        print vtkGridType
        if vtkGridType == "DATASET STRUCTURED_POINTS" + "\n":
            uniform_y = True
            print "Structured_points - uniform wall normal grid resolution\n"
        else:
            uniform_y = False
            print "Rectilinear_grid - non-uniform wall normal grid resolution\n"

        # Read DIMENSIONS (nx,ny,nz)
        dbyte = 0
        vtkDimensions = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkDimensions += dbyte
        print vtkDimensions
        dimstr = string.splitfields(vtkDimensions)
        dim = [string.atoi(dimstr[1]),string.atoi(dimstr[2]),string.atoi(dimstr[3])]

        if uniform_y:
            # Read origin
            dbyte = 0
            vtkOrigin = ""
            while dbyte != chr(10):
                dbyte = g.read(1)
                vtkOrigin += dbyte
            print vtkOrigin
            origstr = string.splitfields(vtkOrigin)
            orig = [float(dimstr[1]),float(dimstr[2]),float(dimstr[3])]

            # Read spacing
            dbyte = 0
            vtkSpacing = ""
            while dbyte != chr(10):
                dbyte = g.read(1)
                vtkSpacing += dbyte
            print vtkSpacing
            spacingstr = string.splitfields(vtkSpacing)
            spacing = [float(spacingstr[1]),float(spacingstr[2]),float(spacingstr[3])]

        else:
            # Read X_COORDINATES
            dbyte = 0
            vtkXcoord = ""
            while dbyte != chr(10):
                dbyte = g.read(1)
                vtkXcoord += dbyte
            print vtkXcoord
            xstr = string.splitfields(vtkXcoord)
            print xstr
            nx = int(xstr[1])
            print "nx:" + str(nx)

            try:
                print "Reading x-coordinates..."
                dbyte2 = g.read(nx*4)
                fdatax = array('f')
                fdatax.fromstring(dbyte2)
#                fdatax.byteswap()
                xcoord = fdatax
#                fdata.fromfile(g,nx)
#                xcoord = struct.unpack('<f',g.read(nx*4))[nx]
                print "Done."
#                print fdatax
            except struct.error:
                print "xcoord error"

            # Read Y_COORDINATES
            dbyte = 0
            vtkYcoord = ""
            while dbyte != chr(10):
                dbyte = g.read(1)
                vtkYcoord += dbyte
            print vtkYcoord
            ystr = string.splitfields(vtkYcoord)
            print ystr
            ny = int(ystr[1])
            print "ny:" + str(ny)

            try:
                print "Reading y-coordinates..."
                dbytey = g.read(ny*4)
                fdatay = array('f')
                fdatay.fromstring(dbytey)
#                fdatay.byteswap()
                ycoord = fdatay
                print "Done."
#                print fdatay
            except struct.error:
                print "ycoord error"

            # Read Z_COORDINATES
            dbyte = 0
            vtkZcoord = ""
            while dbyte != chr(10):
                dbyte = g.read(1)
                vtkZcoord += dbyte
            print vtkZcoord
            zstr = string.splitfields(vtkZcoord)
            print zstr
            nz = int(zstr[1])
            print "nz:" + str(nz)

            try:
                print "Reading z-coordinates..."
                dbytez = g.read(nz*4)
                fdataz = array('f')
                fdataz.fromstring(dbytez)
#                fdataz.byteswap()
                zcoord = fdataz
                print "Done."
#                print fdataz
                
            except struct.error:
                print "zcoord error"

        # Read Metadata section header
        dbyte = 0
        vtkMetadataTitle = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetadataTitle += dbyte
#        print vtkMetadataTitle

        # Read Metadata TIME
        dbyte = 0
        vtkMetaTimeHeader = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetaTimeHeader += dbyte
#        print vtkMetaTimeHeader
        #dimstr = string.splitfields(vtkDimensions)
        #dim = [string.atoi(dimstr[1]),string.atoi(dimstr[2]),string.atoi(dimstr[3])]

        dbyte = 0
        vtkMetaTime = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetaTime += dbyte
#        print vtkMetaTime
        #dim = [string.atoi(dimstr[1]),string.atoi(dimstr[2]),string.atoi(dimstr[3])]

        # Read Metadata FILE
        dbyte = 0
        vtkMetaFileHeader = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetaFileHeader += dbyte
#        print vtkMetaFileHeader

#        vtkMetaFile = g.read(100)
#        print vtkMetaFile

        # Read Metadata PATH
        dbyte = 0
        vtkMetaPathHeader = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetaPathHeader += dbyte
#        print vtkMetaPathHeader

#        vtkMetaPath = g.read(100)
#        print vtkMetaPath

        # Read Metadata GENERATED
        dbyte = 0
        vtkMetaGeneratedHeader = ""
#        while dbyte != chr(10):
#            dbyte = g.read(1)
#            vtkMetaGeneratedHeader += dbyte
#        print vtkMetaGeneratedHeader

#        vtkMetaGenerated = g.read(100)
#        print vtkMetaGenerated

        # Read POINTDATA
        dbyte = 0
        vtkPointdata = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkPointdata += dbyte
        print vtkPointdata
        pointstr = string.splitfields(vtkPointdata)
        npoint = string.atoi(pointstr[1])

        # Read data header
        dbyte = 0
        vtkDataHeader = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkDataHeader += dbyte
        print vtkDataHeader

        # Read data header
        dbyte = 0
        vtkDataHeaderLookup = ""
        while dbyte != chr(10):
            dbyte = g.read(1)
            vtkDataHeaderLookup += dbyte
        print vtkDataHeaderLookup

        try:
            print "Reading lambda2..."
            dbytel = g.read(npoint*4)
#            fdatad = array('f')
            fdatad.fromstring(dbytel)
            fdatad.byteswap()

            print "Done."
#            print fdatad
        except struct.error:
            print "ycoord error"

        g.close()

    except:
        print "VTK file could not be read to the end ",

    print st

    file_content = []
    file_content.append(vtkHeader)
    file_content.append(vtkTitle)
    file_content.append(vtkBinary)
    file_content.append(vtkGridType)
    file_content.append(vtkDimensions)
    file_content.append(vtkXcoord)
    file_content.append(vtkYcoord)
    file_content.append(vtkZcoord)
    file_content.append(vtkPointdata)
    file_content.append(vtkDataHeader)
    file_content.append(vtkDataHeaderLookup)

    return file_content

#
# Write to file on legacy VTK format
#
# fname   Filename to be read
# tdata   Type of data (lambda2, vel_u ,vel_v, vel_w, omega_x, omega_y, omega_z, ufluct )
#def readBinaryVTK(fname,tdata,x,y,z,data):
def writeBinaryVTK(fname,data,xcoord,ycoord,zcoord,l2):

    # Use dir name 
    st  = "--- Write binary VTK file ------------------------------------------" + "\n"
    st += "Filename                      ::: " + fname + "\n"
    st += "--------------------------------------------------------------------" + "\n"

    g = open(fname,"wb")

    # vtkHeader
    g.write(data.pop(0))
    # vtkTitle
    g.write(data.pop(0))
    # vtkBinary
    g.write(data.pop(0))
    # vtkGridType
    g.write(data.pop(0))
    # vtkDimensions
    data.pop(0)
    g.write("DIMENSIONS  " + str(len(xcoord)) + " " + str(len(ycoord)) + "  " + str(len(zcoord)) + "\n")

    g.write("X_COORDINATES  " + str(len(xcoord)) + " float" + "\n")
    data.pop(0)
    g.write(xcoord)

    g.write("Y_COORDINATES  " + str(len(ycoord)) + " float" + "\n")
    data.pop(0)
    g.write(ycoord)

    g.write("Z_COORDINATES   " + str(len(zcoord)) + " float" + "\n")
    data.pop(0)
    g.write(zcoord)

    g.write("POINT_DATA  " + str(len(xcoord)*len(ycoord)*len(zcoord)) + "\n")
    data.pop(0)
    print "POINT_DATA  " + str(len(xcoord)*len(ycoord)*len(zcoord)) + "\n"

    # vtkDataHeader
    g.write(data.pop(0))
    # vtkDataHeaderLookup
    g.write(data.pop(0))

    g.write(l2)

    g.close()


#
# Write to file on XML VTK format
#
# fname   Filename to be read
# tdata   Type of data (lambda2, vel_u ,vel_v, vel_w, omega_x, omega_y, omega_z, ufluct )
#def readBinaryVTK(fname,tdata,x,y,z,data):
def writeBinaryXMLVTK(fname,data,xcoord,ycoord,zcoord,l2):

    # Use dir name 
    st  = "--- Read binary VTK file -------------------------------------------" + "\n"
    st += "Filename                      ::: " + fname + "\n"
    st += "--------------------------------------------------------------------" + "\n"

    g = open(fname,"wb")

    g.write('<?xml version="1.0"?>\n')
    g.write('<VTKFile type="RectilinearGrid" version="0.1" byte_order="BigEndian">\n')
    wholeExtentStr = '<RectilinearGrid WholeExtent="' + " 0 " + str(len(xcoord)-1) + " 0 " + str(len(ycoord)-1) + " 0 " + str(len(zcoord)-1) + '">\n'
    g.write(wholeExtentStr)
    pieceExtentStr = '<Piece Extent="' + " 0 " + str(len(xcoord)-1) + " 0 " + str(len(ycoord)-1) + " 0 " + str(len(zcoord)-1) + '">\n'
    g.write(pieceExtentStr)
    g.write('<PointData>\n')
    #g.write('<DataArray type="Float32" Name="lambda2" format="binary">\n')
    ##g.write(l2)
    #g.write(base64.standard_b64encode(l2))
    #g.write('</DataArray>\n')
    g.write('</PointData>\n')
    g.write('<CellData>\n')
    g.write('<DataArray type="Float32" Name="lambda2" format="binary">\n')
    #g.write(l2)
    g.write(base64.standard_b64encode(l2))
    g.write('</DataArray>\n')
    g.write('</CellData>\n')

    g.write('<Coordinates>\n')

    # x
    g.write('<DataArray type="Float32" Name="x" format="binary">\n')
    g.write(base64.standard_b64encode(xcoord))
    g.write('</DataArray>\n')

    # y
    g.write('<DataArray type="Float32" Name="y" format="binary">\n')
    #g.write(ycoord)
    g.write(base64.standard_b64encode(ycoord))
    g.write('</DataArray>\n')

    # z
    g.write('<DataArray type="Float32" Name="z" format="binary">\n')
    #g.write(zcoord)
    g.write(base64.standard_b64encode(zcoord))
    g.write('</DataArray>\n')

    g.write('</Coordinates>\n')

    g.write('</Piece>\n')
    g.write('</RectilinearGrid>\n')
    g.write('</VTKFile>\n')

    g.close()

#
# Read scaling data
#
# fname   Filename of scaling data
#
def readScalingData(fname):

    f = open(fname,"rt")
    lines = f.readlines()
    f.close()

    ind = 0
    sdata = array('f')
    for line in lines:

        sline = string.splitfields(line.strip("\n")," ")
        sdata.append(1.0/float(sline[4]))
        ind = ind + 1

    return sdata


############################################################################
# Main program
############################################################################

#
# If arguments supplied run only those test cases
#
print "Converting..."
args = sys.argv[1:]

#if sys.argv([]< 2):
#    print 'Too few arguments' 
#    print 'Give input data filename, start iteration, end iteration, output filename'    
#    break
#elif

#
# Read all datasets
#
#datafiles = ['field.0950.u.vtk','field.0950.u.vtk']
#datafiles = ['t150931_11001-12000.lam.vtk']
#datafiles = ['t150931_11001-12000.lam.vtk','t150931_12001-13000.lam.vtk','t150931_13001-14000.lam.vtk']
#datafiles = ['t150931_00001-01000.lam.vtk','t150931_01001-02000.lam.vtk','t150931_02001-03000.lam.vtk','t150931_03001-04000.lam.vtk','t150931_04001-05000.lam.vtk','t150931_05001-06000.lam.vtk','t150931_06001-07000.lam.vtk','t150931_07001-08000.lam.vtk','t150931_08001-09000.lam.vtk','t150931_09001-10000.lam.vtk','t150931_10001-11000.lam.vtk','t150931_11001-12000.lam.vtk','t150931_12001-13000.lam.vtk','t150931_13001-14000.lam.vtk']
#datafiles = ['t150931_00001-01000.u.vtk','t150931_01001-02000.u.vtk','t150931_02001-03000.u.vtk','t150931_03001-04000.u.vtk','t150931_04001-05000.u.vtk','t150931_05001-06000.u.vtk','t150931_06001-07000.u.vtk']
#datafiles = ['t150931_07001-08000.u.vtk','t150931_08001-09000.u.vtk','t150931_09001-10000.u.vtk','t150931_10001-11000.u.vtk','t150931_11001-12000.u.vtk','t150931_12001-13000.u.vtk','t150931_13001-14000.u.vtk']

datafiles = ['t150931_00001-01000.lam.vtk','t150931_01001-02000.lam.vtk','t150931_02001-03000.lam.vtk','t150931_03001-04000.lam.vtk','t150931_04001-05000.lam.vtk','t150931_05001-06000.lam.vtk','t150931_06001-07000.lam.vtk']
#datafiles = ['t150931_07001-08000.lam.vtk','t150931_08001-09000.lam.vtk','t150931_09001-10000.lam.vtk','t150931_10001-11000.lam.vtk','t150931_11001-12000.lam.vtk','t150931_12001-13000.lam.vtk','t150931_13001-14000.lam.vtk']

#output = 'us_merged_07.vtk'
#output = 'us_merged_714.vtk'
output = 'l2s_merged_07.vtk'
#output = 'l2s_merged_714.vtk'
#output = 'merged1.vtr'
nfiles = len(datafiles)

ifile = 0
xcoord = array('f')
ycoord = array('f')
zcoord = array('f')
xcoordglob = array('f')
nx = 0
ny = 0
nz = 0

# Read scaling data
sdata = readScalingData('l2scaling_DNS.data')

for datafile in datafiles:

    fdatad = array('f')
    ifile = ifile + 1

#    DNSdata = readBinaryVTK(datafile,"lambda2",xcoord,ycoord,zcoord)
    DNSdata = readBinaryVTK(datafile,"lambda2",ifile)
#    DNSdata = readBinaryVTK(datafile,"u_vel",ifile)

    if ifile == 1:
        nxg = nfiles * nx
        ny = ny
        nz = nz
#        l2 = numpy.arange(float(nz*ny*nxg*1.0)).reshape(nz,ny,nxg) #zeros(nz*ny*nxg)
        l2 = numpy.zeros((nz,ny,nxg), numpy.float32)
        print l2.dtype.name 
#        l2 = l2*0.0
#        l2 = numpy.reshape(nz,ny,nxg)
#        l2dataglob = [ [ [ 0.0 for z in range(0,nz)] for y in range(0,ny)] for x in range(0,nxg)]


    for iz in range(0, nz):
        for iy in range(0, ny):
            l2[iz][iy][(ifile-1)*nx:ifile*nx]=fdatad[iz*ny*nx+iy*nx:iz*ny*nx+iy*nx+nx]
            l2[iz][iy][(ifile-1)*nx:ifile*nx] = l2[iz][iy][(ifile-1)*nx:ifile*nx]*sdata[(ifile-1)*nx:ifile*nx]

        print "============== iz: " + str(iz) + " iy: " + str(iy) + " file: " + str(ifile)

    for ix in range(0, nx):
        xcoordglob.append(xcoord[ix])

#        l2dataglob[(ifile-1)*nx:ifile*nx-1][0:ny-1][0:nz-1]=l2data[0:nx-1][0:ny-1][0:nz-1]=l2data

#
# At the point where to create the actual data, loop over all loaded data and write to file
#
l2.byteswap(True)
DNSdata = writeBinaryVTK(output,DNSdata,xcoordglob,ycoord,zcoord,l2)
#DNSdata = writeBinaryXMLVTK(output,DNSdata,xcoordglob,ycoord,zcoord,l2)

print "Conversion finished."

#if __name__ == "__main__":
