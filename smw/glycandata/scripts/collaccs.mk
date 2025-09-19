ROOT=/data/projects/GlyGen
GLYRES=$(ROOT)/PyGly/scripts/glyres.py

DATESTAMP=$(shell date "+%Y%m%d")
TIMESTAMP_CMD = date +'%Y-%m-%d %H:%M:%S'
SHELL = /bin/bash -o pipefail
.DELETE_ON_ERROR:

GLYCANDATA=$(ROOT)/PyGly/smw/glycandata
SCRIPTS=$(GLYCANDATA)/scripts
DATA=../data

TAXIDS=9606 10090 10116 10114 111108 11116 11103 3052230 \
	   63746 694009 2697049 7227 4932 9823 44689 9031 3702 \
	   9913 7955 10029

GLYGEN_ACCESSIONS=$(DATA)/glygen_accessions.txt
GLYGEN_FORCE_ACCESSIONS=$(DATA)/glygen_force_accessions.txt
GLYGEN_EXCLUDE_ACCESSIONS=$(DATA)/glygen_exclude_accessions.txt
GLYGEN_REQ_ACCESSIONS=$(DATA)/glygen_req_accessions.txt
GLYGEN_MANUAL_ACCESSIONS=$(DATA)/glygen_manual_accessions.txt
GLYGEN_NEW_ACCESSIONS=$(DATA)/glygen_new_accessions.txt
EXTRA_ACCESSIONS=$(DATA)/extra_accessions.txt

GNOME_RAW=$(DATA)/gnome_subsumption_raw.txt

BASECOMP_TXT=$(DATA)/basecomplist.txt
BASECOMP_LOG=$(DATA)/basecomplist.log
BASECOMP_DONE=$(DATA)/basecomplist.$(DATESTAMP).done

UCKBDATA_DONE=$(DATA)/uckbdata.$(DATESTAMP).done

GLYTOUCAN_REPLACED=$(DATA)/glytoucan_replaced.txt
GLYTOUCAN_ARCHIVED=$(DATA)/glytoucan_archived.txt
GLYTOUCAN_ALLACC=$(DATA)/glytoucan_allacc.txt
GLYCOSMOS_ALLACC=$(DATA)/glycosmos_allacc.txt
GLYTOUCANACC_DONE=$(DATA)/glytoucanacc.$(DATESTAMP).done

.PHONY: all startall clean

.DEFAULT_GOAL: all

all: startall $(BASECOMP_DONE) $(UCKBDATA_DONE) $(GLYTOUCANACC_DONE) \
     $(GLYGEN_ACCESSIONS) $(GLYGEN_NEW_ACCESSIONS) $(EXTRA_ACCESSIONS)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

startall:
	@echo "[`$(TIMESTAMP_CMD)`]: Start all"

clean: 
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	rm -f $(BASECOMP_DONE) $(UCKBDATA_DONE) $(GLYTOUCANACC_DONE) \
	      $(GLYGEN_FORCE_ACCESSIONS) $(GLYGEN_ACCESSIONS) \
		  $(GLYGEN_NEW_ACCESSIONS) $(EXTRA_ACCESSIONS)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

$(GLYGEN_EXCLUDE_ACCESSIONS):
	touch $@

$(GLYGEN_FORCE_ACCESSIONS):
	touch $@

$(GLYGEN_ACCESSIONS): $(GLYGEN_EXCLUDE_ACCESSIONS) $(GLYGEN_FORCE_ACCESSIONS) $(GLYTOUCANACC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	$(GLYRES) GlyGen allglycans | \
	fgrep -v -f $(GLYGEN_EXCLUDE_ACCESSIONS) | \
	fgrep -v -f $(GLYTOUCAN_ARCHIVED) | \
	sort -u > $@
	
$(GLYGEN_REQ_ACCESSIONS): $(GLYGEN_MANUAL_ACCESSIONS) $(GLYGEN_FORCE_ACCESSIONS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cp -f $< $@
	$(GLYRES) GlycoMotifDevNoCache allmotifs GGM | awk '$$1 ~ /^GGM.000/ {print $$2}' >> $@
	$(GLYRES) GlyTouCanNoCache bytaxa $(TAXIDS) | awk '{print $$1}' >> $@
	$(GLYRES) GlyCosmosNoCache bytaxa $(TAXIDS) | awk '{print $$1}' >> $@
	$(GLYRES) GlyGenSourceFile allsourcegtc | awk '{print $$2}' >> $@

$(GLYGEN_NEW_ACCESSIONS): $(GLYGEN_REQ_ACCESSIONS) $(GLYGEN_EXCLUDE_ACCESSIONS) $(GLYTOUCANACC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cat $(GLYGEN_REQ_ACCESSIONS) | \
	fgrep -v -f $(GLYGEN_ACCESSIONS) | \
	fgrep -v -f $(GLYGEN_EXCLUDE_ACCESSIONS) | \
	fgrep -v -f $(GLYTOUCAN_ARCHIVED) | \
	sort -u > $@

$(EXTRA_ACCESSIONS):  $(GNOME_RAW) $(GLYGEN_ACCESSIONS) $(GLYGEN_NEW_ACCESSIONS) $(GLYTOUCANACC_DONE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	./get_extra_accessions.py $(GNOME_RAW) $(GLYGEN_ACCESSIONS) $(GLYGEN_NEW_ACCESSIONS) | \
	fgrep -v -f $(GLYTOUCAN_ARCHIVED) | \
	sort -u > $@

$(BASECOMP_DONE): 
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	./getbasecomplist.py \* >$(BASECOMP_TXT) 2>$(BASECOMP_LOG)
	cd $(DATA); ./splitbasecomp.sh
	touch $@

$(UCKBDATA_DONE):
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	./uckbdata.sh
	touch $@

$(GLYTOUCANACC_DONE):
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	$(GLYRES) GlyCosmosNoCache replaced | bash -c '(read -r; echo "$$REPLY"; sort -u)' > $(GLYTOUCAN_REPLACED)
	$(GLYRES) GlyCosmosNoCache archived | bash -c '(read -r; echo "$$REPLY"; sort -u)' > $(GLYTOUCAN_ARCHIVED)
	$(GLYRES) GlyTouCanNoCache allaccessions | sort -u > $(GLYTOUCAN_ALLACC)
	$(GLYRES) GlyCosmosNoCache allaccessions | sort -u > $(GLYCOSMOS_ALLACC)
	touch $@
