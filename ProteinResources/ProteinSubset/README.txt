
Proteins accessions read from the input file(s) can be retained or dropped
according to protein subsets. This might be used to focus attention on a
specific protein class (such as plasma membrane proteins). This will, by
design, lead to loss of peptide identifications from the results. Parsimony
analysis is still applied. 

Protein mappings are formatted as CSV or Excel (.xls) files. The first
column is the input protein accessions to retain. The column headings are
required but arbitrary.

Example (from ipipm.csv)

Accession
IPI00000006
IPI00000041
IPI00000070
IPI00000149
IPI00000190
IPI00000203
IPI00000357
IPI00000425
IPI00000513
...

Currently provided files include:

ipipm.csv       - Lists IPI Human accessions associated with gene onotology cellular
                  compoenent plasma membrane using a custom GO Slim. Only
                  experimental associations between UniProt and GO
                  are used.

ipipmall.csv    - Lists IPI Human accessions associated with gene onotology
                  cellular compoenent plasma membrane using a custom GO
                  Slim. All (inferred and experimental) associations
                  between UniProt and GO are used.

ipiuniqpm.csv   - Lists IPI Human accessions associated _only_ with gene
                  onotology cellular compoenent plasma membrane using a
                  custom GO Slim. Only experimental associations between
                  UniProt and GO are used.

ipiuniqpmall.csv - Lists IPI Human accessions associated _only_ with
                  gene onotology cellular compoenent plasma membrane
                  using a custom GO Slim. All (inferred and experimental)
                  associations between UniProt and GO are used.

sppm.csv        - Lists UniProt Human accessions associated with gene
                  onotology cellular compoenent plasma membrane using a
                  custom GO Slim. Only experimental associations between
                  UniProt and GO are used.

sppmall.csv     - Lists UniProt Human accessions associated with gene
                  onotology cellular compoenent plasma membrane using
                  a custom GO Slim. All (inferred and experimental)
                  associations between UniProt and GO are used.

spuniqpm.csv    - Lists UniProt Human accessions associated _only_
                  with gene onotology cellular compoenent plasma membrane
                  using a custom GO Slim. Only experimental associations
                  between UniProt and GO are used.

spuniqpmall.csv  - Lists UniProt Human accessions associated _only_ with
                  gene onotology cellular compoenent plasma membrane
                  using a custom GO Slim. All (inferred and experimental)
                  associations between UniProt and GO are used.
