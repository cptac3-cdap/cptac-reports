#!/bin/sh
NAME="CPTAC-CDAP-Reports"
VER=`cat ${NAME}.data/VERSION.txt`
OS=`uname`
AR=`uname -m`
XX="linux-$AR"
mkdir -p build dist
PYFILES="
loader1.py sample1.py labeldecoy.py removedecoy.py genemap1.py
protmap1.py cptac_peptide_summary1.py cptac_assay_annotate1.py parsnip1.py
mayu.py cptac_protein_summary.py proremove1.py peptrypterm1.py
cptac_spectral_count.py cptac_precursor_area.py cptac_tmt10.py cptac_peptideion_tmt10.py
protortho1.py peptrypterm1.py psmdb-shell.py make_sample.py progrpfilt1.py
cptac_spectral_count.py cptac_precursor_area.py pepremap1.py summary1.py
"
./build.py $NAME ${NAME}-${VER}.${XX} $PYFILES
mv build/$NAME-${VER}.${XX}.tgz dist

