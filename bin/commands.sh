#!/bin/sh

if [ -L $0 ] ; then
    DIR=$(dirname $(readlink -f $0)) ;
else
    DIR=$(dirname $0) ;
fi ;
DIR=$DIR/..

uname=`uname -a`
cygwin=`echo $uname | grep -i -c cygwin`

CMD=`basename $0`
CMD=`echo $CMD | sed -e 's/\..*$//'`

echo "Install firectory is $DIR"

case $CMD in
conformermatch)
	CLASS="conformermatch.ConformerMatcher"
	;;
sdfparser)
	CLASS="conformermatch.RocsSdfParser"
	;;
conformeroverlay)
	CLASS="conformermatch.ConformerOverlay"
	;;
gape)
	CLASS="com.cairn.gape.Superposition"
	;;
pharmsearch)
	CLASS="com.cairn.gape.PharmSearch"
	;;
grips)
	CLASS="com.cairn.gape.RigidPharmSearch"
	;;
*)
	echo "Unknown command $CMD";
	exit
	;;
esac

lib="$DIR/target"
echo "running Java class $CLASS"
if [ "$cygwin" = '1' ]
then
    winlib=`cygpath -w $lib`
    java -server -Xmx1024M -ea -cp "$winlib\gape-1.0-SNAPSHOT-all.jar;" $CLASS "$@"
else 
    java -server -Xmx1024M -ea -cp "$lib/gape-1.0-SNAPSHOT-all.jar:" $CLASS "$@"
fi


