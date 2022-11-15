#!/bin/sh
NAME="CPTAC-CDAP-Reports"
VER=`cat ${NAME}.data/VERSION.txt`
DCCVER="$NAME-${VER}"
OS=`uname`
AR=`uname -m`
XX="linux-$AR"
rm -f dist/$NAME-${VER}.${XX}.tgz
./build.sh
if [ ! -f dist/$NAME-${VER}.${XX}.tgz ]; then
  echo "Problem!" 1>&2
  exit 1
fi
( cd dist; md5sum $DCCVER.*.tgz > $DCCVER.md5 ; touch $DCCVER.txt )
if [ "$1" ]; then 
  for comment in "$@"; do 
    echo "* $comment" >> dist/$DCCVER.txt
  done
fi
# gh release delete "$DCCVER" -y
# git push --delete origin "refs/tags/$DCCVER"
# git tag --delete "$DCCVER"
# gh release create -F dist/$DCCVER.txt "$DCCVER" dist/$DCCVER.*.tgz dist/$DCCVER.md5
for a in dist/$DCCVER.*.tgz; do
  a1=`basename $a`
  rclone copyto $a cptac-s3:cptac-cdap.georgetown.edu/$a1
  b1=`echo $a1 | sed "s/-$VER//"`
  aws --profile cptac s3api put-object --bucket cptac-cdap.georgetown.edu --key "$b1" --website-redirect-location "/$a1"
done
