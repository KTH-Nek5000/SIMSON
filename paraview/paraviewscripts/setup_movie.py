#!/software/apps/paraview/3.8.1/bin/pvpython
#########################################################################
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
#########################################################################
import vtk
import subprocess
import string
import array
import sys
import math
import circle

############################################################################
#
# Functions for setting up camera trajectories and other input parameter
# files.
#
############################################################################

#
# Extract movie globals from movie_??.dat parameter file
#
def compute_movie_globals(mfile):

    print " -- Reading movie parameter file : " + mfile

    st=""
    #
    # Open input movie file
    #
    f = open(mfile,"rt")
    lines = f.readlines()
    f.close()
    iline = 0

    #
    # Read header
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    mversion = float(lines[iline].strip("\n"))
    print "      %f"%(mversion)
    iline = iline + 1

    #
    # Title
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    title = lines[iline].strip("\n")
    print "      " + title
    iline = iline + 1

    #
    # Input directory
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    indir = lines[iline].strip("\n")
    print "      " + indir
    iline = iline + 1

    #
    # Scene files
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    sfiles = string.splitfields(lines[iline].strip("\n")," ")
    print "      " + str(sfiles)
    iline = iline + 1

    #
    # Single frame file
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    framefile = lines[iline].strip("\n")
    print "      " + framefile
    iline = iline + 1

    #
    # Trajectory and global movie parameter files
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    dfiles = string.splitfields(lines[iline].strip("\n")," ")
    print "      " + dfiles[0] + "  " + dfiles[1]
    tfile = dfiles[0]
    gfile = dfiles[1]
    iline = iline + 1

    #
    # Image file name
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    ifile = lines[iline].strip("\n")
    print "      " + ifile
    iline = iline + 1

    #
    # Parallel projection
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    parallelProjection = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + parallelProjection
    st += parallelProjection + "\n"

    #
    # Down sample
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    downSample = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + downSample
    st += downSample + "\n"

    #
    # Cut off
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    cutOff = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + cutOff
    st += cutOff + "\n"

    #
    # Background color
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    bColor = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + bColor
    st += bColor + "\n"

    #
    # View size
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    vSize = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + vSize
    st += vSize + "\n"

    #
    # Background meshes
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    mesh = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + mesh
    st += mesh + "\n"

    #
    # Background mesh color
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    mColor = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + mColor
    st += mColor + "\n"

    #
    # Initial image numbner
    #
    print "    " + lines[iline].strip("\n")
    iline = iline + 1
    nstartimage = lines[iline].strip("\n")
    iline = iline + 1
    print "       " + nstartimage
    st += nstartimage + "\n"

    #
    # Write output parameter file
    #
    g = open(gfile,'w')
    g.writelines(st)
    g.close()

    return indir,framefile,gfile,tfile,ifile,sfiles

#
# Function that reads scene file and generates camera position and focal
# point orbits data file which can be read from main program driver.py
#
def compute_orbit(fnames,gname,indir):

    st = ""
    for fname in fnames:

        print " -- Reading scene parameter file : " + fname

        #
        # Open input scene file
        #
        f = open(fname,"rt")
        lines = f.readlines()
        f.close()
        iline = 0

        #
        # Read header
        #
        print "    " + lines[iline].strip("\n")
        iline = iline + 1
        sversion = float(lines[iline].strip("\n"))
        print "      %f"%(sversion)
        iline = iline + 1

        #
        # Read number of frames
        #
        print "    " + lines[iline].strip("\n")
        iline = iline + 1
        nframe = int(lines[iline])
        print "      %d"%(nframe)
        if nframe == 1:
            print "Single plot!"
        #
        # Get colormap filename
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")
        iline = iline + 1
        colormap = lines[iline]
        colormap = colormap.strip("\n")
        print "      " + colormap

        #
        # Read camera positions, spline interpolate
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")

        iline = iline + 1
#        ncpos = int(lines[iline].strip("\n"))
        tline = string.splitfields(lines[iline].strip("\n"))
        ncpos = int(tline[0])
        cposClosed = int(tline[1])
        if ncpos == 1:
            circ = True
            print "    Camera moves on circle"
        else:
            circ = False
            print "    Camera follows interpolated based on point list"        

        cpos = vtk.vtkPoints()

        cposSplineX = vtk.vtkCardinalSpline()
        cposSplineY = vtk.vtkCardinalSpline()
        cposSplineZ = vtk.vtkCardinalSpline()

        if circ:
            cposSplineX.SetClosed(True)
            cposSplineY.SetClosed(True)
            cposSplineZ.SetClosed(True)
        else:
            cposSplineX.SetClosed(cposClosed)
            cposSplineY.SetClosed(cposClosed)
            cposSplineZ.SetClosed(cposClosed)

        if circ:
            iline = iline + 1

            xpos,ypos,zpos = circle.get_circ_points(lines[iline],nframe)

            for icpos in range(0, nframe):
                cpos.InsertPoint(icpos,xpos[icpos],ypos[icpos],zpos[icpos])
                cposSplineX.AddPoint(icpos, xpos[icpos])
                cposSplineY.AddPoint(icpos, ypos[icpos])
                cposSplineZ.AddPoint(icpos, zpos[icpos])  
        else:
                
            for icpos in range(0, ncpos):
                iline = iline + 1
                sline = string.splitfields(lines[iline].strip("\n"))
                cpos.InsertPoint(icpos,float(sline[0]),float(sline[1]),float(sline[2]))
                
                cposSplineX.AddPoint(icpos, float(sline[0]))
                cposSplineY.AddPoint(icpos, float(sline[1]))
                cposSplineZ.AddPoint(icpos, float(sline[2]))
                
                print "      %d" %(icpos) + " : " + sline[0] + " , " + sline[1] + " , "+ sline[2]
                    
        #
        # Read camera focus points, spline interpolate
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")

        iline = iline + 1
#        nfpos = int(lines[iline].strip("\n"))
        tline = string.splitfields(lines[iline].strip("\n"))
        nfpos = int(tline[0])
        fposClosed = int(tline[1])
        fpos = vtk.vtkPoints()

        fposSplineX = vtk.vtkCardinalSpline()
        fposSplineY = vtk.vtkCardinalSpline()
        fposSplineZ = vtk.vtkCardinalSpline()

        fposSplineX.SetClosed(fposClosed)
        fposSplineY.SetClosed(fposClosed)
        fposSplineZ.SetClosed(fposClosed)

        for ifpos in range(0, nfpos):
            iline = iline + 1
            sline = string.splitfields(lines[iline].strip("\n"))
            fpos.InsertPoint(ifpos,float(sline[0]),float(sline[1]),float(sline[2]))
            print "      %d" %(ifpos) + " : " + sline[0] + " , " + sline[1] + " , "+ sline[2]

            fposSplineX.AddPoint(ifpos, float(sline[0]))
            fposSplineY.AddPoint(ifpos, float(sline[1]))
            fposSplineZ.AddPoint(ifpos, float(sline[2]))

        #
        # Read zoom, spline interpolate
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")

        iline = iline + 1
#        nzoom = int(lines[iline].strip("\n"))
        tline = string.splitfields(lines[iline].strip("\n"))
        nzoom = int(tline[0])
        zoomClosed = int(tline[1])
        zoom = vtk.vtkPoints()
        zoomSplineX = vtk.vtkCardinalSpline()

        zoomSplineX.SetClosed(zoomClosed)

        for izoom in range(0, nzoom):
            iline = iline + 1
            zoom.InsertPoint(izoom,float(lines[iline]),0.0,0.0)
            print "      %d" %(izoom) + " : " + lines[iline].strip("\n")

            zoomSplineX.AddPoint(izoom, float(lines[iline]))

        #
        # Read camera view up vector
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")

        iline = iline + 1
#        nvpos = int(lines[iline].strip("\n"))
        tline = string.splitfields(lines[iline].strip("\n"))
        nvpos = int(tline[0])
        vposClosed = int(tline[1])

        vpos = vtk.vtkPoints()

        vposSplineX = vtk.vtkCardinalSpline()
        vposSplineY = vtk.vtkCardinalSpline()
        vposSplineZ = vtk.vtkCardinalSpline()

        vposSplineX.SetClosed(vposClosed)
        vposSplineY.SetClosed(vposClosed)
        vposSplineZ.SetClosed(vposClosed)

        for ivpos in range(0, nvpos):
            iline = iline + 1
            sline = string.splitfields(lines[iline].strip("\n"))
            vpos.InsertPoint(ivpos,float(sline[0]),float(sline[1]),float(sline[2]))

            vposSplineX.AddPoint(ivpos, float(sline[0]))
            vposSplineY.AddPoint(ivpos, float(sline[1]))
            vposSplineZ.AddPoint(ivpos, float(sline[2]))

            print "      %d" %(ivpos) + " : " + sline[0] + " , " + sline[1] + " , "+ sline[2]

        #
        # Read file name and data file animation parameters
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")

        iline = iline + 1
        nstart = int(lines[iline])
        print "      nstart : " + lines[iline].strip("\n")

        iline = iline + 1
        nstep = int(lines[iline])
        print "      nstep : " + lines[iline].strip("\n")

        iline = iline + 1
        nfield = int(lines[iline].strip("\n"))
        print "      nfield : " + lines[iline].strip("\n")

        iline = iline + 1
        nend = int(lines[iline])
        print "      nend : " + lines[iline].strip("\n")

        iline = iline + 1
        nrestart = int(lines[iline])
        print "      nrestart : " + lines[iline].strip("\n")

        iline = iline + 1
        sname = string.splitfields(lines[iline],"*")
#        print "Length of sname string" + str(len(sname))
#        if sname[0].strip("\n") == "none":
#            sname[0] = sname[0].strip("\n")
#            sname.extend(["\n"])
#        else:
#            sname[1] = sname[1].strip("\n")
#        print "      Name of field : " + str(sname)

        if len(sname) == 1:
            if sname[0].strip("\n") == "none":
                sname[0] = sname[0].strip("\n")
                sname.extend(["\n"])
            else:
               sname[0] = sname[0].strip("\n")
#               sname.extend(["\n"])
        else:
            sname[1] = sname[1].strip("\n")
        print "      Name of field : " + str(sname)

        #
        # Read spline parameter (not used at the moment)
        #
        iline = iline + 1
        print "    " + lines[iline].strip("\n")
        iline = iline + 1
        closedSpline = int(lines[iline].strip("\n"))

        #
        # Close spline curves
        #
#        if closedSpline:
#            fposSplineX.SetClosed(True)
#            fposSplineY.SetClosed(True)
#            fposSplineZ.SetClosed(True)
#            cposSplineX.SetClosed(True)
#            cposSplineY.SetClosed(True)
#            cposSplineZ.SetClosed(True)
#            zoomSplineX.SetClosed(True)
#            vposSplineX.SetClosed(True)
#            vposSplineY.SetClosed(True)
#            vposSplineZ.SetClosed(True)
#        else:
#            fposSplineX.SetClosed(False)
#            fposSplineY.SetClosed(False)
#            fposSplineZ.SetClosed(False)
#            cposSplineX.SetClosed(False)
#            cposSplineY.SetClosed(False)
#            cposSplineZ.SetClosed(False)
#            zoomSplineX.SetClosed(False)
#            vposSplineX.SetClosed(False)
#            vposSplineY.SetClosed(False)
#            vposSplineZ.SetClosed(False)
        #
        # Generate file name of data to be read
        #
        hname = []
        ifile = nstart - nstep
        for iframe in range(0, nframe):
            if ( iframe % nfield == 0 ):
                ifile = ifile + nstep
                if ifile > nend:
                    ifile = nrestart
                    print "Starting from first data field again..."
            if sname[0] == "none":
                nname = sname[0]
            elif len(sname) == 1:
                nname = indir + sname[0]
            else:
                nname = indir + sname[0] + "%(nnum)04d" % {'nnum':ifile} + sname[1]
            #    print "The name of the data file is " + nname
            hname.extend([nname])

        #
        # Generate a polyline for the spline
        #
        cposPoints = vtk.vtkPoints()
        fposPoints = vtk.vtkPoints()
        zoomPoints = vtk.vtkPoints()
        vposPoints = vtk.vtkPoints()

        #
        # Interpolate x, y and z by using the three spline filters and
        # create new points
        #
        if (nframe > 1):
            cpos=[0.0,0.0,0.0]
            fpos=[0.0,0.0,0.0]
            zoom=[0.0]
            vpos=[0.0,0.0,0.0]
            for iframe in range(0, nframe):
                if circ:
                    t = iframe
                else:
                    if cposClosed:
                        t = float(ncpos)/float(nframe)*iframe
                    else:
                        t = (ncpos-1.0)/(nframe-1.0)*iframe
                        
#            print ":::: %f" %(t) + " ncpos %d" %(ncpos)
            
                cpos[0] = cposSplineX.Evaluate(t)
                cpos[1] = cposSplineY.Evaluate(t)
                cpos[2] = cposSplineZ.Evaluate(t)
                cposPoints.InsertPoint(iframe, cpos[0], cpos[1], cpos[2])
                
                if fposClosed:
                    t = float(nfpos)/float(nframe)*iframe
                else:
                    t = (nfpos-1.0)/(nframe-1.0)*iframe
                fpos[0] = fposSplineX.Evaluate(t)
                fpos[1] = fposSplineY.Evaluate(t)
                fpos[2] = fposSplineZ.Evaluate(t)
                fposPoints.InsertPoint(iframe, fpos[0], fpos[1], fpos[2])
                    
                if zoomClosed:
                    t = float(nzoom)/float(nframe)*iframe
                else:
                    t = (nzoom-1.0)/(nframe-1.0)*iframe
                zoom[0] = zoomSplineX.Evaluate(t)
                zoomPoints.InsertPoint(iframe, zoom[0], 0.0, 0.0)

                if vposClosed:
                    t = float(nvpos)/float(nframe)*iframe
                else:
                    t = (nvpos-1.0)/(nframe-1.0)*iframe
                vpos[0] = vposSplineX.Evaluate(t)
                vpos[1] = vposSplineY.Evaluate(t)
                vpos[2] = vposSplineZ.Evaluate(t)
                vposPoints.InsertPoint(iframe, vpos[0], vpos[1], vpos[2])

                st += str(cpos[0]) + " " + str(cpos[1]) + " " + str(cpos[2]) + " "
                st += str(fpos[0]) + " " + str(fpos[1]) + " " + str(fpos[2]) + " "
                st += str(zoom[0]) + " "
                st += str(vpos[0]) + " " + str(vpos[1]) + " " + str(vpos[2]) + " "
                st += hname[iframe] + " " + colormap + "\n"
        else:
            ccpos = cpos.GetPoint(0)
            ffpos = fpos.GetPoint(0)
            zzoom = zoom.GetPoint(0)
            vvpos = vpos.GetPoint(0)
            st += str(ccpos[0]) + " " + str(ccpos[1]) + " " + str(ccpos[2]) + " "
            st += str(ffpos[0]) + " " + str(ffpos[1]) + " " + str(ffpos[2]) + " "
            st += str(zzoom[0]) + " "
            st += str(vvpos[0]) + " " + str(vvpos[1]) + " " + str(vvpos[2]) + " "
            st += hname[iframe] + " " + colormap + "\n"            

    #
    # Analyze trajectory file and count how many fields that are the same grouped together and add
    # that information in a new column. This is used to determine how many loops that should be
    # done before exit single_frame... and read in a new field
    #
    stlines = string.splitfields(st,"\n")
    nlines = len(stlines)
    print "Number of lines in trajectory file " + str(nlines)

    #    neqlines_list = []
    #    iline = 0
    #    icolumn = 10

    #    while iline < nlines - 1:

    #        split_string = string.splitfields(stlines[iline].strip("\n")," ")
    #        column_string = split_string[icolumn]

    #        neqlines = countLines(column_string,stlines[iline+1:],icolumn)
    #        iline = iline + neqlines
    #        neqlines_list.extend([neqlines])

    #print ":::" + neqlines_list

    #
    # Add index column in the rightmost column
    #
    for istr in range(0,nlines-1):
        stlines[istr] = stlines[istr] + " " + str(istr) + "\n"

    #
    # Write output file
    #
    g = open(gname,'w')
    g.writelines(stlines)
    g.close()

    #    return neqlines_list

#
# Returns the number of lines that have identical string in icolumn
#
def neqlines(gname):

    neqlines_list = []
    iline = 0
    icolumn = 10

    g = open(gname,'rt')
    stlines= g.readlines()
    g.close()
    nlines = len(stlines)

    while iline < nlines - 1:

        split_string = string.splitfields(stlines[iline].strip("\n")," ")
        column_string = split_string[icolumn]

        neqlines = countLines(column_string,stlines[iline+1:],icolumn)
        iline = iline + neqlines
        neqlines_list.extend([neqlines])

    return neqlines_list

#
# Returns the number of lines that have identical string in icolumn
#
def countLines(column_string,stlines,icolumn):

    neqlines = 1
    if len(stlines) >= 1:
        for istr in range(0,len(stlines)):
            split_string_next = string.splitfields(stlines[istr].strip("\n")," ")
            if column_string == split_string_next[icolumn]:
                neqlines = neqlines + 1
            else:
                break

    return neqlines


#
# Function that generates background vtk-format grid files with
# given box dimensions and grid spacings. Note, all input numbers real!
#
def gridgen(xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,bname):

    #
    # Create spacing and file header
    #
    xdim=abs(xmax-xmin)
    ydim=abs(ymax-ymin)
    zdim=abs(zmax-zmin)
    nx = int(xdim/dx) + 1
    ny = int(ydim/dy) + 1
    nz = int(zdim/dz) + 1

    xstr = ""
    for ix in range(0, nx):
        xstr = xstr + " " + str(xmin+float(ix)*dx)
    ystr = ""
    for iy in range(0, ny):
        ystr = ystr + " " + str(ymin+float(iy)*dy)
    zstr = ""
    for iz in range(0, nz):
        zstr = zstr + " " + str(zmin+float(iz)*dz)

    st  = "# vtk DataFile Version 2.0\n"
    st += "Grid: " + bname + "\n"
    st += "ASCII\n"
    st += "DATASET RECTILINEAR_GRID\n"

    #
    # Create xy file
    #
    stx = st
    stx+= "DIMENSIONS %d %d %d" %(nx,ny,2) + "\n"
    stx+= "X_COORDINATES %d" %(nx) + " double\n"
    stx+= xstr + "\n"
    stx+= "Y_COORDINATES %d" %(ny)  + " double\n"
    stx+= ystr + "\n"
    stx+= "Z_COORDINATES 2 double\n"
    stx+= str(zmin) + " " + str(zmax)+ "\n"

    bfile = "grid-xy-" + bname + ".vtk"
    g = open(bfile,'w')
    g.writelines(stx)
    g.close()

    #
    # Create xz file
    #
    sty = st
    sty+= "DIMENSIONS %d %d %d" %(nx,2,nz) + "\n"
    sty+= "X_COORDINATES %d" %(nx) + " double\n"
    sty+= xstr + "\n"
    sty+= "Y_COORDINATES 2 double\n"
    sty+= str(ymin) + " " + str(ymax)+ "\n"
    sty+= "Z_COORDINATES %d" %(nz) + " double\n"
    sty+= zstr + "\n"

    bfile = "grid-xz-" + bname + ".vtk"
    g = open(bfile,'w')
    g.writelines(sty)
    g.close()

    #
    # Create yz file
    #
    stz = st
    stz+= "DIMENSIONS %d %d %d" %(2,ny,nz) + "\n"
    stz+= "X_COORDINATES 2 double\n"
    stz+= str(xmin) + " " + str(xdim)+ "\n"
    stz+= "Y_COORDINATES %d" %(ny) + " double\n"
    stz+= ystr + "\n"
    stz+= "Z_COORDINATES %d" %(nz) + " double\n"
    stz+= zstr + "\n"

    bfile = "grid-yz-" + bname + ".vtk"
    g = open(bfile,'w')
    g.writelines(stz)
    g.close()
