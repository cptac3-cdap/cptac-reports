#!/bin/sh

for f in `find "ProteinResources" -size +50M`; do
  echo "$f"
  gzip -9 -c "$f" | split -d -b 40M - "$f".gz.
  # rm -f "$f"
done
