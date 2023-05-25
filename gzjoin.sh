#!/bin/sh

for f in `find ProteinResources -name "*.gz.00"`; do
  BASE=`echo "$f"  | sed 's/\.gz\.00$//'`
  echo "$BASE"
  cat "$BASE".gz.[0-9][0-9] | gunzip -c > "$BASE"
  # rm -f "$BASE".gz.[0-9][0-9]
done
