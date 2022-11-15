#!/bin/sh
NAME=${1:-`basename $PWD`}
BASE=`dirname $0`
BASE=`readlink -f "$BASE"`
if [ ! -d mzIdentML ]; then
  ( cat | fmt ) 1>&2 <<EOF
  Directory mzIdentML is missing. 
  The mzIdentML directory should contain the mzid.gz files to be loaded.
EOF
  exit 1;
fi
if [ ! -f TARGET.txt ]; then
  ( cat | fmt ) 1>&2 <<EOF
  TARGET.txt is missing. 
  TARGET.txt should specify one of the following workflows:
`sed -n -e 's/^%.\([a-z0-9_]*\): *\\\/\1/p' $BASE/Makefile | awk '{print "    "$1}'`
EOF
  exit 1;
fi
if [ ! -f FDR.txt ]; then
  ( cat | fmt ) 1>&2 <<EOF
  FDR.txt is missing. FDR.txt should contain the max. spectral FDR
  (in percent) for a peptide to be considered for protein inference.
EOF
  exit 1;
fi 
if [ ! -f SPECCNT.txt ]; then
  ( cat | fmt ) 1>&2 <<EOF
  SPECCNT.txt is missing. SPECCNT.txt should contain the
  minimum spectral count for a peptide to be considered for protein inference.
EOF
  exit 1;
fi 
if [ ! -f "$NAME".sample.csv ]; then
  ( cat | fmt ) 1>&2 <<EOF 
  $NAME.sample.csv is missing. $NAME.sample.csv contains the experimental
  design that maps samples and labels to PSM files.
EOF
  exit 1;
fi 
# exec >>"$NAME.log" 2>&1 
make -k -f $BASE/Makefile BASE=$BASE "$NAME".`cat TARGET.txt` >>"$NAME.log" 2>&1 &
