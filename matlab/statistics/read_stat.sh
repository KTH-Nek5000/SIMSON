#!/bin/sh


if [ 5 -ne "$#" ]; then
    echo "Usage: $0 xy.stat.file list coord scalar output"
else
    echo $1
#    echo $4 # scalar
#    echo "1" # averaging in channel flow
    echo "n" # filtering
#    echo "y"
#    echo 10
#    echo "0" # change T,1-T
    echo "n" # mean plots
    if [ -r $2 ]; then
	let "count=0"
	while read line; do
            let "count=count+1"
	    echo "$line"
	    echo "y" 
	    echo "8"
	    if [ $count -lt 10 ]; then
		echo "STAT0$count.data"
	    else
		echo "STAT$count.data"
	    fi
	    echo "0"
	    echo "3"
	    echo "1"
	    echo $3   # coordinate
	    echo "0"     # scaling
	    echo "1"
	done < $2
    fi
    echo "0"
#    echo "0"
fi | pxyst

paste -d " " STAT??.data | tr -s " " | cut -d " " -f 2,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61 > $5.10
rm -f STAT??.data

