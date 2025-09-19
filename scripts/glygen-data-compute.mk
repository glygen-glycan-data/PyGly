
ROOT=/data/projects/GlyGen
WORKERS=yion:8,bion:8,proton:8,maldi:8

DATESTAMP=$(shell date "+%Y%m%d")
TIMESTAMP_CMD = date +'%Y-%m-%d %H:%M:%S'
SHELL = /bin/bash -o pipefail
.DELETE_ON_ERROR:

.PHONY: all startall stage1 stage2 stage3 stage4 stage5 glycomotif_nocache \
        stage1start stage2start stage3start stage4start stage5start

GLYCOMOTIF=$(ROOT)/PyGly/smw/glycomotif
GLYCOMOTIFSCR=$(GLYCOMOTIF)/scripts
GLYCOMOTIFDATA=$(GLYCOMOTIF)/data

GLYCOMOTIFDEV_ALIGMENTS=$(GLYCOMOTIFDATA)/computealignments-dev.$(DATESTAMP).tsv
GLYCOMOTIFDEV_ALIGMENTS_LOG=$(GLYCOMOTIFDATA)/computealignments-dev.$(DATESTAMP).log
GLYCOMOTIFDEV_CLASSIFICATIONS=$(GLYCOMOTIFDATA)/computeclass-dev.$(DATESTAMP).tsv
GLYCOMOTIFDEV_RDF=$(GLYCOMOTIFDATA)/glycomotifdev_motifalign.$(DATESTAMP).rdf.gz

all: startall stage1 stage2 stage3 stage4 stage5
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete"

startall:
	@echo "[`$(TIMESTAMP_CMD)`]: Start all"

stage1start:
	@echo "[`$(TIMESTAMP_CMD)`]: Start stage1"

stage2start:
	@echo "[`$(TIMESTAMP_CMD)`]: Start stage2"

stage3start:
	@echo "[`$(TIMESTAMP_CMD)`]: Start stage3"

stage4start:
	@echo "[`$(TIMESTAMP_CMD)`]: Start stage4"

stage5start:
	@echo "[`$(TIMESTAMP_CMD)`]: Start stage5"

glycomotif_nocache:
	- cd $(GLYCOMOTIFSCR); \
	  rm -f .gtc.cache.* .gco.cache.*;

$(GLYCOMOTIFDEV_ALIGMENTS):
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR); \
        ./computealignments2.py -v --workers $(WORKERS) -o "$@" > $(GLYCOMOTIFDEV_ALIGMENTS_LOG) 2>&1
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOMOTIFDEV_CLASSIFICATIONS): $(GLYCOMOTIFDEV_ALIGMENTS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR); \
        ./computeclassification.py < $< > $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOMOTIFDEV_RDF): $(GLYCOMOTIFDEV_ALIGMENTS) $(GLYCOMOTIFDEV_CLASSIFICATIONS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR) && \
	./alignments_tsv2rdf.py $^ | gzip -9 -c > $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

GLYCOMOTIFPROD_ALIGMENTS=$(GLYCOMOTIFDATA)/computealignments-prod.$(DATESTAMP).tsv
GLYCOMOTIFPROD_ALIGMENTS_LOG=$(GLYCOMOTIFDATA)/computealignments-prod.$(DATESTAMP).log
GLYCOMOTIFPROD_CLASSIFICATIONS=$(GLYCOMOTIFDATA)/computeclass-prod.$(DATESTAMP).tsv
GLYCOMOTIFPROD_RDF=$(GLYCOMOTIFDATA)/glycomotif_motifalign.$(DATESTAMP).rdf.gz

$(GLYCOMOTIFPROD_ALIGMENTS): $(GLYCOMOTIFDEV_ALIGMENTS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR); \
        ./computealignments2.py --smwenv PROD -v --workers $(WORKERS) -o "$@" > $(GLYCOMOTIFPROD_ALIGMENTS_LOG) 2>&1
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

#FIXME WHEN CLASSIFICATIONS ARE POPULATED IN PROD
$(GLYCOMOTIFPROD_CLASSIFICATIONS): $(GLYCOMOTIFPROD_ALIGMENTS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR); \
        ./computeclassification.py --smwenv DEV < $< > $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOMOTIFPROD_RDF): $(GLYCOMOTIFPROD_ALIGMENTS) $(GLYCOMOTIFPROD_CLASSIFICATIONS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOMOTIFSCR); \
	./alignments_tsv2rdf.py $^ | gzip -9 -c > $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

GNOMESCR=$(ROOT)/PyGly/scripts
GNOMERAW=$(GNOMESCR)/subsumption.$(DATESTAMP).txt

$(GNOMERAW):
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GNOMESCR); \
	./gnome_compute.py > $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

stage1: stage1start glycomotif_nocache $(GLYCOMOTIFDEV_RDF) $(GLYCOMOTIFPROD_RDF) $(GNOMERAW)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete" 

#
# STAGE2: LOAD THE GLYCOMOTIF TRIPLESTORES
#

SMWROOT=/data/projects/smw/docker

GLYCOMOTIFSMW=$(SMWROOT)/glycomotif
GLYCOMOTIFSMW_LOADED=$(GLYCOMOTIFSMW)/loadts.$(DATESTAMP).done

GLYCOMOTIFDEVSMW=$(SMWROOT)/glycomotifdev
GLYCOMOTIFDEVSMW_LOADED=$(GLYCOMOTIFDEVSMW)/loadts.$(DATESTAMP).done

$(GLYCOMOTIFDEVSMW_LOADED): $(GLYCOMOTIFDEV_RDF)
	@echo "[`$(TIMESTAMP_CMD)`]: Start load of GlycoMotif (DEV) TripleStore" 
	cp $< $(GLYCOMOTIFDEVSMW);
	ssh trypsin "cd $(GLYCOMOTIFDEVSMW); ../bin/dumprdf.sh";
	cd $(GLYCOMOTIFDEVSMW); \
	TRIP1=`ls glycomotifdev.*.rdf.gz | tail -n 1`; \
	TRIP2=$(<F); \
	ssh trypsin "cd $(GLYCOMOTIFDEVSMW); ../bin/loadts.sh $$TRIP1 $$TRIP2"
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: Load of GlycoMotif (DEV) TripleStore complete"

$(GLYCOMOTIFSMW_LOADED): $(GLYCOMOTIFPROD_RDF)
	@echo "[`$(TIMESTAMP_CMD)`]: Start load of GlycoMotif (PROD) TripleStore" 
	cp $< $(GLYCOMOTIFSMW);
	ssh trypsin "cd $(GLYCOMOTIFSMW); ../bin/dumprdf.sh";
	cd $(GLYCOMOTIFSMW); \
	TRIP1=`ls glycomotif.*.rdf.gz | tail -n 1`; \
	TRIP2=$(<F); \
	ssh trypsin "cd $(GLYCOMOTIFSMW); ../bin/loadts.sh $$TRIP1 $$TRIP2"
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: Load of GlycoMotif (PROD) TripleStore complete" 

stage2: stage1 stage2start $(GLYCOMOTIFDEVSMW_LOADED) $(GLYCOMOTIFSMW_LOADED)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete" 

#
# STAGE3: REBUILD GLYCOTREE 
#

GLYCOTREE=$(ROOT)/glycoTree
GLYCOTREEDEV=$(ROOT)/glycoTreeDev
GLYCOTREEDEV_DATA=$(GLYCOTREEDEV)/data

GLYCOTREE=$(ROOT)/glycoTree
GLYCOTREE_DATA=$(GLYCOTREE)/data

GLYCOTREEDEV_DATADL=$(GLYCOTREEDEV_DATA)/glycotree_data.$(DATESTAMP).done
GLYCOTREEDEV_BUILD=$(GLYCOTREEDEV)/glycotree_build.$(DATESTAMP).done
GLYCOTREEDEV_POPULATE=$(GLYCOTREEDEV)/glycotree_populate.$(DATESTAMP).done
GLYCOTREEDEV_EXPORTS=$(GLYCOTREEDEV)/glycotree_exports.$(DATESTAMP).done

GLYCOTREE_DATADL=$(GLYCOTREE_DATA)/glycotree_data.$(DATESTAMP).done
GLYCOTREE_BUILD=$(GLYCOTREE)/glycotree_build.$(DATESTAMP).done
GLYCOTREE_POPULATE=$(GLYCOTREE)/glycotree_populate.$(DATESTAMP).done
GLYCOTREE_EXPORTS=$(GLYCOTREE)/glycotree_exports.$(DATESTAMP).done

# FIXME, SHOULD PULL FROM GLYCOMOTIFSMW (PROD)
$(GLYCOTREEDEV_DATADL): $(GLYCOMOTIFDEVSMW_LOADED)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREEDEV_DATA) && \
	./dl.sh >& dl.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREEDEV_BUILD): $(GLYCOTREEDEV_DATADL)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREEDEV) && \
	./build_all.sh clear >& ./build_all.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREEDEV_POPULATE): $(GLYCOTREEDEV_BUILD)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREEDEV) && \
	ssh trypsin "cd $(GLYCOTREEDEV); ./scripts/populate.sh" && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREEDEV_EXPORTS): $(GLYCOTREEDEV_POPULATE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREEDEV) && \
	ssh trypsin "cd $(GLYCOTREEDEV); ./scripts/exports.sh" && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREE_DATADL): $(GLYCOMOTIFSMW_LOADED)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREE_DATA) && \
	./dl.sh >& dl.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREE_BUILD): $(GLYCOTREE_DATADL)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREE) && \
	./build_all.sh clear >& ./build_all.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREE_POPULATE): $(GLYCOTREE_BUILD)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREE) && \
	ssh trypsin "cd $(GLYCOTREE); ./scripts/populate.sh" && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

$(GLYCOTREE_EXPORTS): $(GLYCOTREE_POPULATE)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)" 
	cd $(GLYCOTREE) && \
	ssh trypsin "cd $(GLYCOTREE); ./scripts/exports.sh" && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete" 

stage3: stage2 stage3start $(GLYCOTREEDEV_EXPORTS) $(GLYCOTREE_EXPORTS)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete" 

APIFRAMEWORK=$(ROOT)/APIFramework
GLYMAGE_SRC=$(APIFRAMEWORK)/src/Application/Glymage
GLYMAGE_IMAGES=$(GLYMAGE_SRC)/glymage_images.$(DATESTAMP).done

$(GLYMAGE_IMAGES): $(GLYCOTREE_POPULATE) $(GLYCOTREEDEV_POPULATE) $(GLYCOMOTIFSMW_LOADED) $(GLYCOMOTIFDEVSMW_LOADED)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cd $(GLYMAGE_SRC) && \
	./makeimages.sh >& makeimages.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete"

stage4: stage3 stage4start $(GLYMAGE_IMAGES)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete" 

GLYCANDATA=$(ROOT)/PyGly/smw/glycandata
GLYCANDATA_SCR=$(GLYCANDATA)/scripts
GLYCANDATA_DATA=$(GLYCANDATA)/data
GLYCANDATA_IMAGES=$(GLYCANDATA_SCR)/glycandata_images.$(DATESTAMP).done
GLYCANDATA_ACCS=$(GLYCANDATA_SCR)/glycandata_accs.$(DATESTAMP).done
GLYCANDATA_BUILD=$(GLYCANDATA_SCR)/glycandata_build.$(DATESTAMP).done
GLYCANDATA_EXPORT=$(GLYCANDATA_SCR)/glycandata_export.$(DATESTAMP).done

$(GLYCANDATA_IMAGES): $(GLYMAGE_IMAGES)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cd $(GLYCANDATA_SCR) && \
	./makeimages1.sh >& makeimages1.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete"

$(GLYCANDATA_ACCS): $(GNOMERAW)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cp $< $(GLYCANDATA_DATA)/gnome_subsumption_raw.txt
	cd $(GLYCANDATA_SCR) && \
	./collaccs.sh --nolog -d $(DATESTAMP) -f >& collaccs.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete"

$(GLYCANDATA_BUILD): $(GLYCANDATA_ACCS)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cd $(GLYCANDATA_SCR) && \
	./build.sh --nolog -d $(DATESTAMP) -f >& build.log && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete"

$(GLYCANDATA_EXPORT): $(GLYCANDATA_BUILD) $(GLYCANDATA_IMAGES)
	@echo "[`$(TIMESTAMP_CMD)`]: Start $(@F)"
	cd $(GLYCANDATA_SCR) && \
	./exports.sh cache >& exports.log && \
	cp -f images.tsv images-*tbz* ../export && \
	touch $@
	@echo "[`$(TIMESTAMP_CMD)`]: $(@F) complete"

stage5: stage4 stage5start $(GLYCANDATA_EXPORT)
	@echo "[`$(TIMESTAMP_CMD)`]: $@ complete" 
