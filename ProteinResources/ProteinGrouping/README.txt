
Proteins are grouped after parsimony analysis, but before final counts are tallied. 

Protein groups are formatted as CSV files or as Excel (.xls) files. The
first column are the protein accessions, the second is the group
accession or name, while the third is the group description (optional).
Column headers are required, but arbitrary. The group accession/name
column header is used in the output report. Protein accessions may be
associated with more than one group.

Example (from sp2gocc.csv):

Accession,CellularComponent,Description
A0AVT1,GO:0005737,cytoplasm
A0EJG6,GO:0005737,cytoplasm
A0JLT2,GO:0005634,nucleus
A1L3X0,GO:0005783,endoplasmic reticulum
A1L491,GO:0005634,nucleus
A1L491,GO:0005737,cytoplasm
A2A274,GO:0005634,nucleus
A2A274,GO:0005739,mitochondrion
A2A288,GO:0005634,nucleus
...

Currently provided files include:

sp2goallcc.csv  - UniProt Human accessions associated with gene ontology
                  cellular components: cytoplasm, nucleous, endoplasmic
                  reticulum, mitochondrion, and Golgi apparatus using a
                  custom GO Slim. All (inferred and experimental)
                  UniProt to GO associations are used.

sp2gocc.csv     - UniProt Human accessions associated with gene ontology
                  cellular components: cytoplasm, nucleous, endoplasmic
                  reticulum, mitochondrion, and Golgi apparatus using a
                  custom GO Slim. Only experimental associations between
                  UniProt and GO are used.

ipi2goallcc.csv - IPI Human accessions associated with gene ontology
                  cellular components: cytoplasm, nucleous, endoplasmic
                  reticulum, mitochondrion, and Golgi apparatus using a
                  custom GO Slim. All (inferred and experimental)
                  UniProt to GO associations are used.

ipi2gocc.csv    - IPI Human accessions associated with gene ontology
                  cellular components: cytoplasm, nucleous, endoplasmic
                  reticulum, mitochondrion, and Golgi apparatus using a
                  custom GO Slim. Only experimental associations between
                  UniProt and GO are used.

up2membrane_kw.csv - Human UniProt accessions associated with
                  SwissProt keyword (CellMembrane, Membrane,
                  Transmembrane, InnerMitochondrialMembrane,
                  OuterMitochondrialMembrane) and subcellular location
                  (CellMembrane, CellMembrane:Experimental, Membrane,
                  Membrane:Experimental) annotations.
