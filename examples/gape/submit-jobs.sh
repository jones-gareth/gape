#!/bin/sh

for d in *
do
    if test -d $d -a -f $d/xtal.sdf
    then

	cd $d
	job="job_$d"
	job_name=`printf %15.15s omga_$d`
	test -f $job && qsub -q short -N $job_name -l nodes=1 ./$job
        cd ..

    fi
done
