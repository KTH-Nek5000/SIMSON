#!/bin/tcsh
# ***********************************************************************
# 
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision: 510 $
#
# ***********************************************************************

foreach f (*.tiff)
echo doing $f
convert -quality 100 $f $f:r.jpg
rm -f $f
end
