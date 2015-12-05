


all: CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo taxonomy/gi_taxid_nucl.dmp library
	Rscript generate-card-db.R
	@echo Now run 
	@echo      kraken-build --build --db .
	@echo to build the database, and 
	@echo      kraken --db . library/* | kraken-report --db . - 
	@echo to ensure everything went fine.

library:
	mkdir library

taxonomy:
	mkdir taxonomy

taxonomy/gi_taxid_nucl.dmp: taxonomy
	printf "1\t1\n" > taxonomy/gi_taxid_nucl.dmp

CARD-files:
	mkdir CARD-files

CARD-files/AT-genes.fa: CARD-files
	CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo

CARD-files/AROtags.txt: CARD-files
	CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo
	
CARD-files/aro.obo: CARD-files
	CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo

clean:
	rm -r library/
	rm -r taxonomy/
	rm -r CARD-files/
