#!/bin/sh
set -x
PROG0=cptac_peptideion_tmt10
DIR=`dirname $0`
if [ -f "$DIR/$PROG0.py" ]; then
    PROG="$DIR/$PROG0.py"
elif [ -f "$DIR/$PROG0.sh" ]; then
    PROG="$DIR/$PROG0.sh"
elif [ -f "$DIR/$PROG0" ]; then
    PROG="$DIR/$PROG0"
else
    echo "Can't find cptac_peptideion_tmt10." 1>&2
    exit 1
fi

FIRST=1
for AAGRP in A B C D E F G H I J K L M N O P Q R S T U V W X Y Z; do
    if [ $FIRST = 1 ]; then
        "$PROG" $@ --peptidentermaas "$AAGRP" --allratios
	FIRST=0
    else
        "$PROG" $@ --peptidentermaas "$AAGRP" --append --allratios
    fi
done
