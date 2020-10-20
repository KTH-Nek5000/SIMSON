for f in $1-?????
do
mv $f `echo $f | sed 's/'"$1"'/'"$2"'/'`
done

