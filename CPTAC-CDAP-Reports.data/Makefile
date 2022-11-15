
DIR=$(BASE)
COMMON=
ACC=RefSeqUniProtIsoformBroadArtifact
DECOYPREFIX=XXX_
GENEMAP:=$(BASE)/ProteinMapping/prhumanmouse2gene.csv
GENEMAP1:=$(BASE)/ProteinMapping/prhuman2gene.csv
ORTHMAP:=$(BASE)/ProteinMapping/mousegene2humangene.csv
ORTHMAP1:=$(BASE)/ProteinMapping/rsmouse2rshuman.csv
ASSAYS:=$(BASE)/ProteinAnnotation/CPTAC_Assays.csv
SPECFDRTHRESH:=$(shell cat FDR.txt)
SPECCOUNTTHRESH:=$(shell cat SPECCNT.txt)
NGENES:=$(shell awk -F, '{print $$2}' $(GENEMAP) | sort -u | wc -l)
NGENES1:=$(shell awk -F, '{print $$2}' $(GENEMAP1) | sort -u | wc -l)

.SECONDARY:
.PRECIOUS: %.psm %.xlsx %.csv %.tsv
.IGNORE: %.psm

%.xenograft_label_free: \
     %.xenogene.orth.pars.gene.sum.tsv \
     %.xenogene.orth.pars.gene.cnt.tsv \
     %.xenogene.orth.pars.gene.area.tsv \
     %.xenogenegrp.pept.sum.tsv \

	@date
	@echo "Done"

%.human_label_free: \
     %.humgene.pars.gene.sum.tsv \
     %.humgene.pars.gene.cnt.tsv \
     %.humgene.pars.gene.area.tsv \
     %.humgenegrp.pept.sum.tsv \

	@date
	@echo "Done"


%.xenograft_itraq: \
     %.xenogenegrp.orth.pars.gene.sum.tsv \
     %.xenogenegrp.orth.pars.gene.itraq.tsv \
     %.xenogenegrp.pept.sum.tsv \

	@date
	@echo "Done"

%.phospho_xenograft_itraq: \
     %.xenogenegrp.pept.sum.tsv \
     %.xenogenegrp.orth.pars.gene.sum.tsv \
     %.xenogenegrp.orth.pars.undoorth.rs.porth.ppars.undoporth.tryp.phosphosite.itraq.tsv \
     %.xenogenegrp.orth.pars.undoorth.rs.porth.ppars.undoporth.tryp.phosphopep.itraq.tsv \

	@date
	@echo "Done"

%.glyco_xenograft_itraq: \
     %.xenogenegrp.pept.sum.tsv \
     %.xenogenegrp.orth.pars.gene.sum.tsv \
     %.xenogenegrp.orth.pars.undoorth.rs.porth.ppars.undoporth.tryp.glycosite.itraq.tsv \
     %.xenogenegrp.orth.pars.undoorth.rs.porth.ppars.undoporth.tryp.glycopep.itraq.tsv \

	@date
	@echo "Done"


%.human_itraq: \
     %.humgenegrp.pars.gene.sum.tsv \
     %.humgenegrp.pars.gene.itraq.tsv \
     %.humgenegrp.pept.sum.tsv \

	@date
	@echo "Done"

%.phospho_human_itraq: \
     %.humgenegrp.pept.sum.tsv \
     %.humgenegrp.pars.gene.sum.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphosite.itraq.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphopep.itraq.tsv \

	@date
	@echo "Done"

%.glyco_human_itraq: \
     %.humgenegrp.pept.sum.tsv \
     %.humgenegrp.pars.gene.sum.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.glycosite.itraq.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.glycopep.itraq.tsv \

	@date
	@echo "Done"

%.xenograft_tmt10: \
     %.xenogenegrp.orth.pars.nodecoy.gene.sum.tsv \
     %.xenogenegrp.orth.pars.nodecoy.gene.tmt10.tsv \
     %.xenogenegrp.orth.pars.nodecoy.gene.tmt10.nc.tsv \
     %.xenogenegrp.nodecoy.pept.sum.tsv \
     %.xenogenegrp.orth.pars.gene.sum.tsv \
     %.xenogenegrp.orth.pars.gene.tmt10.tsv \
     %.xenogenegrp.orth.pars.gene.tmt10.nc.tsv \
     %.xenogenegrp.pept.sum.tsv \

	@date
	@echo "Done"


%.human_tmt10: \
     %.humgenegrp.pars.gene.sum.tsv \
     %.humgenegrp.pars.gene.tmt10.tsv \
     %.humgenegrp.pars.gene.tmt10.nc.tsv \
     %.humgenegrp.pept.sum.tsv \
     %.humgenegrp.pars.nodecoy.gene.sum.tsv \
     %.humgenegrp.pars.nodecoy.gene.tmt10.tsv \
     %.humgenegrp.pars.nodecoy.gene.tmt10.nc.tsv \
     %.humgenegrp.nodecoy.pept.sum.tsv \

	@date
	@echo "Done"

%.phospho_human_tmt10: \
     %.humgenegrp.pars.gene.sum.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphosite.tmt10.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphosite.tmt10.nc.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphopep.tmt10.tsv \
     %.humgenegrp.pars.rs.ppars.tryp.phosphopep.tmt10.nc.tsv \
     %.humgenegrp.pept.sum.tsv \
     %.humgenegrp.pars.nodecoy.gene.sum.tsv \
     %.humgenegrp.pars.nodecoy.rs.ppars.tryp.phosphosite.tmt10.tsv \
     %.humgenegrp.pars.nodecoy.rs.ppars.tryp.phosphosite.tmt10.nc.tsv \
     %.humgenegrp.pars.nodecoy.rs.ppars.tryp.phosphopep.tmt10.tsv \
     %.humgenegrp.pars.nodecoy.rs.ppars.tryp.phosphopep.tmt10.nc.tsv \
     %.humgenegrp.nodecoy.pept.sum.tsv \

	@date
	@echo "Done"

%.psm: %.sample.csv
	@date
	find -L mzIdentML -name "*.mzid.gz" -exec $(DIR)/loader1 -t 1 -d $@ --pracc $(ACC) -l 7 -q $(COMMON) {} \;
	@date
	$(DIR)/sample1   -d $@ -s $*.sample.csv $(COMMON)
	@date
	$(DIR)/labeldecoy $@ $(DECOYPREFIX)

%.humgene.psm: %.psm
	@date
	$(DIR)/protmap1  -d $< -o $@ -p $(GENEMAP1) --decoyprefix $(DECOYPREFIX) $(COMMON)
	@date
	$(DIR)/cptac_assay_annotate1 $@ $(ASSAYS)

%.xenogene.psm: %.psm
	@date
	$(DIR)/protmap1  -d $< -o $@ -p $(GENEMAP) --decoyprefix $(DECOYPREFIX) $(COMMON)
	@date
	$(DIR)/cptac_assay_annotate1 $@ $(ASSAYS)

%.orth.psm: %.psm
	@date
	$(DIR)/protortho1  -d $< -o $@ --bygene --orthmap $(ORTHMAP) --decoyprefix $(DECOYPREFIX) $(COMMON)

%.porth.psm: %.psm
	@date
	$(DIR)/protortho1  -d $< -o $@ --orthmap $(ORTHMAP1) --decoyprefix $(DECOYPREFIX) $(COMMON)

%.undoorth.psm: %.psm
	@date
	$(DIR)/protortho1  -d $< -o $@ --bygene --undo --decoyprefix $(DECOYPREFIX) $(COMMON)

%.undoporth.psm: %.psm
	@date
	$(DIR)/protortho1  -d $< -o $@ --undo --decoyprefix $(DECOYPREFIX) $(COMMON)

%.humgenegrp.psm: %.psm
	@date
	$(DIR)/genemap1  -d $< -o $@ -g $(GENEMAP1) --decoyprefix $(DECOYPREFIX) $(COMMON)
	@date
	$(DIR)/cptac_assay_annotate1 $@ $(ASSAYS)

%.xenogenegrp.psm: %.psm
	@date
	$(DIR)/genemap1  -d $< -o $@ -g $(GENEMAP) --decoyprefix $(DECOYPREFIX) $(COMMON)
	@date
	$(DIR)/cptac_assay_annotate1 $@ $(ASSAYS)

%.pars.psm: %.psm
	@date
	$(DIR)/parsnip1  -d $< -U 2 -T $(SPECFDRTHRESH) -C $(SPECCOUNTTHRESH) --tieresolution Consistent --pepweights FDR --pracc OrthoGene --bygene --mode Stats,Matrix,Filter $(COMMON) --itersolve 50 --ttotal 15 --noalignments --extension pars > $*.pars.log 
	@date
	$(DIR)/mayu --bygene --summary $(NGENES) $@

define PARAM_PARS
%.pars-$1-$2-$3.txt: %.psm
	@date
	$(DIR)/parsnip1 -d $$< -T $1 -U $2 -C $3 --tieresolution Consistent --pepweights FDR --pracc OrthoGene --mode Stats,Dump $(COMMON) --bygene --itersolve 50 --ttotal 15 --noalignments --extension pars-$1-$2-$3 >& $$*.pars-$1-$2-$3.log
	@date
	$(DIR)/mayu --bygene --accessions $$@ --summary $(NGENES) $$<
endef
threshs := $(shell awk -v START=0.00 -v STOP=1.00 -v STEP=0.001 'BEGIN {for (i=START;i<=(STOP+STEP);i+=STEP) {printf("%.4f ",i)}}' | sed 's/\.*00* / /g')
$(foreach t,$(threshs),$(eval $(call PARAM_PARS,$(t),2,2)))
$(foreach t,$(threshs),$(eval $(call PARAM_PARS,$(t),2,1)))
# $(foreach t,$(threshs),$(eval $(call PARAM_PARS,$(t),1,2)))

%.rs.psm: %.psm
	@date
	$(DIR)/proremove1 -d $< -o $@ -f source -r "RefSeq" $(COMMON)

%.up.psm: %.psm
	@date
	$(DIR)/proremove1 -d $< -o $@ -f source -r "UniProt" $(COMMON)

%.ppars.psm: %.psm
	@date
	$(DIR)/parsnip1  -d $< -U 2 -T $(SPECFDRTHRESH) -C $(SPECCOUNTTHRESH) --tieresolution Consistent --pepweights FDR --pracc OrthoRefSeqUniProtIsoform --mode Stats,Matrix,Filter $(COMMON) --itersolve 50 --noalignments --extension ppars > $*.ppars.log 

%.nodecoy.psm: %.psm
	@date
	$(DIR)/removedecoy -d $< -o $@

%.gene.cnt.tsv: %.psm
	@date
	$(DIR)/cptac_spectral_count  -d $< -o $@ $(COMMON)

%.gene.area.tsv: %.psm
	@date
	$(DIR)/cptac_precursor_area  -d $< -o $@ $(COMMON)

%.gene.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_itraq  -d $< -c POOL --bygene -o $@ $(COMMON)

%.gene.tmt10.tsv: %.psm
	@date
	$(DIR)/cptac_tmt10  -d $< --bygene -o $@ $(COMMON)

%.gene.tmt10.nc.tsv: %.psm
	@date
	$(DIR)/cptac_tmt10  -d $< --bygene --nocorrection -o $@ $(COMMON)

%.tryp.psm: %.psm
	@date
	$(DIR)/peptrypterm1 -d $< --annotate -o $@ $(COMMON)

%.phosphosite.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< -c POOL --report PhosphositeCombination --localized --specific -o $@ $(COMMON)

%.phosphopep.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< -c POOL --report Phosphopeptide --localized --specific -o $@ $(COMMON)

%.phosphosite.tmt10.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_tmt10 -d $< --report PhosphositeCombination --localized --specific -o $@ $(COMMON)

%.phosphopep.tmt10.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_tmt10 -d $< --report Phosphopeptide --localized --specific -o $@ $(COMMON)

%.phosphosite.tmt10.nc.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_tmt10 -d $< --report PhosphositeCombination --localized --nocorrection --specific -o $@ $(COMMON)

%.phosphopep.tmt10.nc.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_tmt10 -d $< --report Phosphopeptide --localized --nocorrection --specific -o $@ $(COMMON)

%.glycosite.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< -c POOL --report GlycositeCombination --noties --specific -o $@ $(COMMON)

%.glycopep.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< -c POOL --report Glycopeptide --noties --specific -o $@ $(COMMON)

%.pepform.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< --report PeptideForm --localized -o $@ $(COMMON)

%.pept.itraq.tsv: %.psm
	@date
	$(DIR)/cptac_peptideion_itraq -d $< --report Peptide --localized -o $@ $(COMMON)

%.gene.sum.tsv: %.psm
	@date
	$(DIR)/cptac_protein_summary --bygene -d $< -o $@ $(COMMON)

%.pept.sum.tsv: %.psm
	@date
	$(DIR)/cptac_peptide_summary1  -d $< -o $@ $(COMMON)

%.export: 
	mkdir -p $*
	-cp *.gene.cnt.tsv $*/$(notdir $*).spectral_counts.tsv
	-cp *.gene.area.tsv $*/$(notdir $*).precursor_area.tsv
	-cp *.gene.sum.tsv $*/$(notdir $*).summary.tsv
	-cp *.gene.itraq.tsv $*/$(notdir $*).itraq.tsv
	-cp *.pept.sum.tsv $*/$(notdir $*).peptides.tsv
