* Freeze HH, Chong JX, Bamshad MJ, Ng BG. Solving glycosylation disorders: fundamental approaches reveal complicated pathways. Am J Hum Genet. 2014 Feb 6;94(2):161-75. doi: 10.1016/j.ajhg.2013.10.024. PMID: 24507773; PMCID: PMC3928651.
  - Human CDG Genes (Supplementary Data, mmc2.xls)
* IMPC download (FTP site):
  - http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/phenotypeHitsPerGene.csv.gz
  - http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz
* Human Phenotype Ontology:
  - https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa
  - https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt
  - https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt
* Human Protein Atlas:
  - https://www.proteinatlas.org/download/proteinatlas.tsv.zip
* Human-Mouse Gene orthologs, per MGI
  - ./listenzymes.py > ../data/enzyme.gene.txt
  - wget -O - -q https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt | acut.py 5 2 1 | fgrep -w -f enzyme.gene.txt > human_mouse_orthologs.txt
