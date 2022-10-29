#!/bin/sh

for d in *
do
    if test -d $d -a -f $d/xtal.sdf
    then
	cd $d

	job="job_$d"
        test -f $job && /bin/rm $job
        dir=`pwd`;
	echo "#!/bin/sh" >> $job
	echo "cd $dir || die" >> $job
        echo "java -server com.cairn.gape.Superposition ../superposition.conf omega.sdf" >> $job
	chmod +x $job

        cd ..
    fi
done

