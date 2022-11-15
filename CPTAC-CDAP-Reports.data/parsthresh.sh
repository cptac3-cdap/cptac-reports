#!/bin/sh
NAME=${2:-`basename $PWD`}
BASE=`dirname $0`
BASE=`readlink -f "$BASE"`
if [ ! -f TARGET.txt ]; then
  ( cat | fmt ) 1>&2 <<EOF
  TARGET.txt is missing. 
  TARGET.txt should specify one of the following workflows:
`sed -n -e 's/^%.\([a-z0-9_]*\): *\\\/\1/p' $BASE/Makefile | awk '{print "    "$1}'`
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
SPECCNT=`cat SPECCNT.txt`
UNIQ=2
exec >>"$NAME.log" 2>&1 
for t in ${1}; do
  # Ensure that temporary PSMDb is created in the local directory (and partition/disk)
  TMPDIR="$PWD"; export TMPDIR
  make -f $BASE/Makefile BASE=$BASE "$NAME".humgenegrp.pars-${t}-2-${SPECCNT}.txt
  fgrep ".pars-" "$NAME.log" | awk '$2 == "-1"' | tr '-' ' ' | awk 'NF == 11 {print $2,$7,$8,$10} NF == 8 {print "1.0",$4,$5,$7}' | sort -n
done
