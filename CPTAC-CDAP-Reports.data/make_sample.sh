#!/bin/sh
DIR=`dirname "$0"`
find -L . -name "*.mzid.gz" | python27 $DIR/make_sample.py $1
