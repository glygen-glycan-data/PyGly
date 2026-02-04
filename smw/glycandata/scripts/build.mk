ROOT=/data/projects/GlyGen
GLYRES=$(ROOT)/PyGly/scripts/glyres.py

DATESTAMP=$(shell date "+%Y%m%d")
TIMESTAMP_CMD = date +'%Y-%m-%d %H:%M:%S'
SHELL = /bin/bash -o pipefail
.DELETE_ON_ERROR:

GLYCANDATA=$(ROOT)/PyGly/smw/glycandata
SCRIPTS=$(GLYCANDATA)/scripts
DATA=../data

CACHE=cache

GLYGEN_ACCESSIONS=$(DATA)/glygen_accessions.txt
GLYGEN_NEW_ACCESSIONS=$(DATA)/glygen_new_accessions.txt
EXTRA_ACCESSIONS=$(DATA)/extra_accessions.txt

ACCESSIONS=$(GLYGEN_ACCESSIONS) $(GLYGEN_NEW_ACCESSIONS) $(EXTRA_ACCESSIONS)

GNOME_RAW=$(DATA)/gnome_subsumption_raw.txt

GLYTOUCAN_REPLACED=$(DATA)/glytoucan_replaced.txt
GLYTOUCAN_ARCHIVED=$(DATA)/glytoucan_archived.txt

UNICARBDATA=$(DATA)/uc2gtc.txt $(DATA)/uc2pubmed.txt $(DATA)/uckbcomp2glytoucan.txt
UNICARBDATA1=$(DATA)/uc2gtc.txt $(DATA)/uc2taxa.txt $(DATA)/uckbcomp2glytoucan.txt

.PHONY: all startall clean

.DEFAULT_GOAL: all

LOADGTC_DONE=.loadgtc.$(DATESTAMP).done
$(LOADGTC_DONE):
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc.py $(CACHE) $(ACCESSIONS)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2GLYCONNECTCOMP_DONE=.loadgtc2glyconnectcomp.$(DATESTAMP).done
$(LOADGTC2GLYCONNECTCOMP_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2glyconnectcomp.py $(CACHE)
	touch $@ 
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADUNICARB_DONE=.loadunicarb.$(DATESTAMP).done
$(LOADUNICARB_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadunicarb.py $(CACHE) $(UNICARBDATA)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

UNICARBTAXID_DONE=.unicarbkb_taxid.$(DATESTAMP).done
$(UNICARBTAXID_DONE): $(LOADUNICARB_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./unicarbkb_taxid.py $(UNICARBDATA1) > $(DATA)/unicarbkb_taxid.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

GLYGENDSTAXID_DONE=.glygends_taxid.$(DATESTAMP).done
$(GLYGENDSTAXID_DONE): 
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./glygends_taxid.py > $(DATA)/glygends_taxid.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADTAXID_DONE=.loadtaxid.$(DATESTAMP).done
$(LOADTAXID_DONE): $(UNICARBTAXID_DONE) $(GLYGENDSTAXID_DONE) $(LOADGTC_DONE) $(LOADGTC2GLYCONNECTCOMP_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadtaxid.py $(CACHE) $(DATA)/glygends_taxid.txt $(DATA)/unicarbkb_taxid.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2PUBCHEM_DONE=.loadgtc2pubchem.$(DATESTAMP).done
$(LOADGTC2PUBCHEM_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2pubchem.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2GLYCOSHAPE_DONE=.loadgtc2glycoshape.$(DATESTAMP).done
$(LOADGTC2GLYCOSHAPE_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2glycoshape.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2MATRIXDB_DONE=.loadgtc2matrixdb.$(DATESTAMP).done
$(LOADGTC2MATRIXDB_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2matrixdb.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2PSIMOD_DONE=.loadgtc2psimod.$(DATESTAMP).done
$(LOADGTC2PSIMOD_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2psimod.py $(CACHE) $(DATA)/psimod2glytoucan.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGLYGEN_DONE=.loadglygen.$(DATESTAMP).done
$(LOADGLYGEN_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadglygen.py $(CACHE) $(DATA)/glygen_accessions.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGTC2GLYCOMEDBEXP_DONE=.loadgtc2glycomedbexport.$(DATESTAMP).done
$(LOADGTC2GLYCOMEDBEXP_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgtc2glycomedbexport.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGDB2GOG_DONE=.loadgdb2gog.$(DATESTAMP).done
$(LOADGDB2GOG_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgdb2gog.py $(CACHE) $(DATA)/gdb2gog.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADEDLAB_DONE=.loadedlab.$(DATESTAMP).done
$(LOADEDLAB_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadedlab.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADGWB_DONE=.loadgwb.$(DATESTAMP).done
$(LOADGWB_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadgwb.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADMOTIF_DONE=.loadmotif.$(DATESTAMP).done
$(LOADMOTIF_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadmotif.py $(CACHE)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADSUBSUBMP_DONE=.loadsubsump.$(DATESTAMP).done
$(LOADSUBSUBMP_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadsubsump.py $(CACHE) $(GNOME_RAW) $(ACCESSIONS)
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADSPECIES_DONE=.loadspecies.$(DATESTAMP).done
$(LOADSPECIES_DONE): $(LOADGTC_DONE) $(LOADTAXID_DONE) $(LOADSUBSUBMP_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadspecies.py $(CACHE) $(GNOME_RAW) 
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADCLASS_DONE=.loadclassification.$(DATESTAMP).done
$(LOADCLASS_DONE): $(LOADGTC_DONE) $(LOADMOTIF_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadclassification.py $(CACHE) $(GNOME_RAW) 
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADTISS_DONE=.loadtissue.$(DATESTAMP).done
$(LOADTISS_DONE): $(LOADGTC_DONE) $(LOADSPECIES_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadtissue.py $(CACHE) $(GNOME_RAW) 
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

LOADNAMES_DONE=.loadnames.$(DATESTAMP).done
$(LOADNAMES_DONE): $(LOADGTC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $@"
	./loadnames.py $(CACHE) UniCarbKB EdwardsLab $(DATA)/uckbcomp2glytoucan.txt 
	./loadnames.py $(CACHE) ShortUniCarbKB EdwardsLab $(DATA)/shortuckbcomp2glytoucan.txt
	./loadnames.py $(CACHE) Byonic EdwardsLab $(DATA)/byonic2glytoucan.txt
	./loadnames.py $(CACHE) ShortComp EdwardsLab $(DATA)/shortcomp2glytoucan.txt
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

ALLDONE=$(LOADGTC_DONE) $(LOADTAXID_DONE) $(LOADGTC2GLYCONNECTCOMP_DONE) \
	$(LOADGTC2MATRIXDB_DONE) $(LOADGTC2PUBCHEM_DONE) \
	$(LOADGTC2GLYCOSHAPE_DONE) $(LOADUNICARB_DONE) $(LOADGTC2PSIMOD_DONE) \
	$(LOADGLYGEN_DONE) $(LOADGTC2GLYCOMEDBEXP_DONE) \
	$(LOADGDB2GOG_DONE) $(LOADEDLAB_DONE) $(LOADGWB_DONE) \
	$(LOADSPECIES_DONE) $(LOADSUBSUBMP_DONE) $(LOADTISS_DONE) \
	$(LOADCLASS_DONE) $(LOADNAMES_DONE) $(LOADMOTIF_DONE) \
	$(UNICARBTAXID_DONE) $(GLYGENDSTAXID_DONE)

all: startall $(ALLDONE)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

clean:
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	rm -f $(ALLDONE)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

startall:
	@echo "[`$(TIMESTAMP_CMD)`]: Start all"
