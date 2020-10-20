#!/bin/tcsh
# ***********************************************************************
# 
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision: 510 $
#
# ***********************************************************************

foreach f (field.02??.u)
lambda2 $f $f:r
end
