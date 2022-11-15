#!/bin/sh

help() {
  CODE=0
  if [ "$1" != "" ]; then
    CODE=1
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

CPTAC-CDAP-Reports/pdcdia_summary.sh [ options ] <base>.elib.peptides.txt

Options:

  -s    Species. One of Human, Mouse. Default: Human.
  -h    Help text.

EOF
  exit $CODE
}

SPECIES="Human"
while getopts ":s:h" o ; do
        case $o in
                s ) SPECIES="$OPTARG";;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac
done

shift $(($OPTIND - 1))

PEPTIDES="$1"
if [ "$PEPTIDES" = "" ]; then
  help "DIA matrix file required."
fi

case "$PEPTIDES" in
  *.elib.peptides.txt) ;;
  *) help "DIA matrix filename must be: <base>.elib.peptides.txt." ;;
esac

case "$SPECIES" in
  Human|Mouse) ;;
  *) help "Bad species string: $SPECIES. Should be one of Human, Mouse." ;;
esac

PROG=`readlink -f "$0"`
DIR=`dirname "$PROG"`
# DIR=`dirname $0`
EXT=""
if [ -f "$DIR/loader1.py" ]; then
  EXT=".py"
fi

BASE=`basename "$PEPTIDES" .elib.peptides.txt`
PEPTIDES=$BASE.elib.peptides.txt

if [ ! -f "$PEPTIDES" ]; then
  help "Bad input file: $PEPTIDES."
fi

if [ "$SPECIES" = "Human" ]; then
  UNIPROT="$DIR/ProteinDatabases/uniprot_human_reference_isoforms_20180717.fasta"
  REFSEQ="$DIR/ProteinDatabases/refseq_human_20180717.fasta"
  GENEMAP="$DIR/ProteinMapping/prhuman2gene-2019-01-22.csv"
elif [ "$SPECIES" = "Mouse" ]; then
  UNIPROT="$DIR/ProteinDatabases/uniprot_mouse_reference_isoforms_20180717.fasta"
  REFSEQ="$DIR/ProteinDatabases/refseq_mouse_20180717.fasta"
  GENEMAP="$DIR/ProteinMapping/prmouse2gene-2019-01-22.csv"
fi

rm -f ${BASE}*.psm ${BASE}.*.pars.*

# Load the PSMs to the PSMDb database
date
echo "Import peptide-sample DIA matrix ..."
set -x
$DIR/loader1${EXT} --alignments None -o $BASE.psm $PEPTIDES || exit 1;
set +x
echo "Done."

# Label each spectrum file with its own name as sample name
date
echo "Annotate spectra filenames with sample identifiers..."
set -x
$DIR/sample1${EXT} -d $BASE.psm --self || exit 1;
set +x
echo "Done."

# Remap to add current CPTAC UniProt and RefSeq protein sequences, keep old RefSeq accessions
date
echo "Remap peptide sequences to current CPTAC UniProt and RefSeq protein sequences..."
set -x
$DIR/pepremap1${EXT} -d $BASE.psm -o $BASE.remap.psm -s "$UNIPROT" --pracc UniProt --orphans -U  || exit 1;
$DIR/pepremap1${EXT} -d $BASE.remap.psm -s "$REFSEQ" --pracc RefSeqAcc -U  || exit 1;
set +x
echo "Done."

# Add genes from gene remap file. 
date
echo "Add CPTAC genes and assocaite proteins with genes..."
set -x
$DIR/genemap1${EXT} -d $BASE.remap.psm -o $BASE.gene.psm -g $GENEMAP || exit 1;
set +x
echo "Done."

# parsimony - "value" is set to 0.0 anyway, so use the standard filtering.
date
echo "Parsimony - no filtering of specific peptides (PSMs), two unshared"
echo "peptides per gene, leave peptide ions observsed in the least samples"
echo "uncovered..."
set -x
( $DIR/parsnip1${EXT} -d $BASE.gene.psm --mode Filter,Matrix,Stats --bygene --pracc Gene -U 2 --noalignments --pepweight Count 2>&1 | tee $BASE.gene.pars.log ) || exit 1
set +x
echo "Done."

# quant - just sum the abundance values from the psms
date
echo "Sum the DIA abundance across peptide ions from each gene..."
set -x
$DIR/cptac_precursor_area${EXT} -d $BASE.gene.pars.psm --quant "psm:abundance" --bygene -o $BASE.gene.pars.dia.tsv || exit 1;
set +x
echo "Done."
