#!/software/apps/paraview/3.8.1/bin/pvpython
#########################################################################
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
#########################################################################

import numpy
from numpy import sin as ns
from numpy import cos as nc
import string
import sys

############################################################################
#
# The transformation matrix Dinv can be found on p.80 of
# "Fundamentals of Astrodynamics" by Bate et al (found on Plan 6)
#
# There is probably a much simpler way to do this...
#
############################################################################

def get_circ_points(linein,nframe):

    params = string.splitfields(linein)
    
    rad = params[3]
    print "   Circle center: ", params[0],params[1],params[2]
    print "   Circle radius  ", params[3]
    if len(params) > 4:
        print "   Start angle (deg. clockwise, from pos. x-axis):", params[4]
        startang = numpy.pi/180.0*float(params[4])
        if len(params) > 5:
            print "   End angle (deg. clockwise, from pos. x-axis):", params[5]
            endang = numpy.pi/180.0*float(params[5])

            if len(params) > 6:
                if len(params) < 9:
                    sys.exit("Must give three components of the normal vector (or none, then a circle in (x,z) plane is assumed)!!!")

                else:

#
# Note the implicit transformation in defining the normal vector!!! 
# Going from SIMSON coordinates to proper Cartesian coordinates
#


                    normvec = numpy.array([float(params[6]),-float(params[8]),float(params[7])])
#                    print normvec
                    nn = normvec/numpy.linalg.norm(normvec)

#                    print nn
                    angl = numpy.arcsin(nn[2])
        
                    if abs(numpy.pi/2.0-abs(angl))< 1e-8:
                        print "Circle in (x,z) plane"
                        Dinv = numpy.eye(3)
                        angt = 0.0
                    else:
                        angt = numpy.arccos(nn[0]/nc(angl))
                        Dinv = numpy.matrix([[ns(angl)*nc(angt), -ns(angt), nc(angl)*nc(angt)],
                                [ns(angl)*ns(angt),  nc(angt), nc(angl)*ns(angt)],  
                                [-nc(angl)         ,     0    , ns(angl)        ]])
#                    print Dinv
            else:
                Dinv = numpy.eye(3)
                    
        else:
            endang = 2*numpy.pi*float(nframe-1)/float(nframe)
            Dinv = numpy.eye(3)
    else:
        startang = 0.0
        endang = 2*numpy.pi*float(nframe-1)/float(nframe)
        Dinv = numpy.eye(3)

    if nframe > 1:    
        angstep = (endang-startang)/float(nframe-1)
    else:
        angstep = 0.0
    xpos=numpy.empty(nframe)
    ypos=numpy.empty(nframe)
    zpos=numpy.empty(nframe)
        
    for icpos in range(0, nframe):

        camangle = startang + float(icpos)*angstep
        vec_circ = numpy.matrix([float(rad)*numpy.cos(camangle),float(rad)*numpy.sin(camangle),0.0])
        vec_xyz = Dinv*numpy.transpose(vec_circ)
#        print vec_xyz
#
# final transformation, back to SIMSON coordinate system
#
        vec_xyz = numpy.matrix([[1.0, 0.0 , 0.0],
                                [0.0, 0.0, -1.0],  
                                [0.0, 1.0, 0.0 ]])*vec_xyz
#        print vec_xyz
        xpos[icpos] = float(params[0]) + vec_xyz[0]
        ypos[icpos] = float(params[1]) + vec_xyz[1]
        zpos[icpos] = float(params[2]) + vec_xyz[2]

    return xpos,ypos,zpos
