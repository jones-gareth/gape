#!/bin/sh

rm */Fitted*

for d in *
do
    if test -d $d -a -f $d/GA_ranked_1.sdf
    then
	cd $d
	echo "processing $d"
	java com.cairn.molecule.IncrementalSuperpositonComparison \
	    GA_ranked_1.sdf xtal.sdf > result.txt
        cd ..
    fi
done

fgrep '**Best**' */result.txt > summary.txt


