#!/bin/sh




i1=1
i2=2
while [ $i1 -lt $i2 ]
do
#    echo waiting...
    sleep 30
    if [ -e field.*.u.written ]
    then
        a=field.*.u.written
        b=`basename $a .written`
        echo copying $b...
#        mv -f $b /scratch/pschlatt/test_video/
	/pdc/vol/hsm/2.0/bin/hput $b
        rm -f $a
#        rm -f $b
    fi
done
