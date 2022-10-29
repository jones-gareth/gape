#!/bin/sh

for d in *
do

   if test -d $d -a -f $d/xtal.sdf
   then
       cd $d
       rm omega*
       /home/packages/openeye/arch/redhat-RHEL5-x64/omega/2.3.2/omega2 -in xtal.sdf -maxconfs 1 -searchff mmff94s omega.sdf
       cd ..
   fi

done
