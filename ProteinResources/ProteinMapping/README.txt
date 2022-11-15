
Proteins accessions read from the input file(s) can be mapped another
protein database's accessions, before parsimony analysis and spectrum
counting. UniProt or IPI accessions might be mapped to genes, IPI
accessions to UniProt, etc. Input protein accessions missing from the mapping
are dropped, which may result in peptide identifications being lost. 

Protein mappings are formatted as CSV or Excel (.xls) files. The first
column is the input protein accessions, the second is the new protein
accession, and the third is the new protein accession's description
(optional). The new protein accession's column header may be used in
the output report. Input protein accessions may map to more than one
new protein accession.

This is particularly useful when used with protein groupings based on a
particular databases accessions (e.g. map from IPI to SwissProt to use UniProt
based GO protein groupings), or to collapse counts for multiple protein
isoforms to a single gene.

Example (from ipihuman2sp.csv):

Accession,Protein
IPI00000001,O95793
IPI00000005,P01111
IPI00000006,P01112
IPI00000012,Q6XR72
IPI00000013,O60911
IPI00000015,Q08170
IPI00000023,P18507
IPI00000024,Q08174
IPI00000026,Q96NX5
...

Currently provided files include:

sp2id.csv       - Maps UniProt accessions to UniProt IDs. 

sp2orggene.csv  - Maps UniProt accessions to Organism/Gene pairs. 

spisoform2sp.csv - Maps SwissProt variant accessions to SwissProt accessions. 

ipihuman2sp.csv - Maps IPI Human accession to SwissProt accessions. 

ipihuman2gene.csv - Maps IPI Human accessions to NCBI gene ids.
