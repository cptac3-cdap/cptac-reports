#!/bin/sh

VERSION=""
if [ "$1" != "" ]; then
  VERSION="-$1"
fi

WD=`readlink -f "$PWD"`

if [ ! -f ./update.sh -o `basename $WD` != "CPTAC-CDAP-Reports" ]; then
  echo "Please execute in CPTAC-CDAP-Reports directory" 1>&2
  exit 1
fi

wget -q -O - http://cptac-cdap.georgetown.edu.s3-website-us-east-1.amazonaws.com/CPTAC-CDAP-Reports$VERSION.linux-x86_64.tgz | tar xzf - -C ..
echo -n "Version: "
cat VERSION.txt
