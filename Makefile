

define result_msg
Now run 
     kraken-build --build --db .
to build the database, and 
     kraken --db . library/* | kraken-report --db . - 
to ensure everything went fine.
endef
export result_msg

all: CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo taxonomy/gi_taxid_nucl.dmp library
	Rscript generate-card-db.R
	@echo "$$result_msg"

library:
	mkdir -p library

taxonomy:
	mkdir -p taxonomy

taxonomy/gi_taxid_nucl.dmp: taxonomy
	printf "1\t1\n" > taxonomy/gi_taxid_nucl.dmp

CARD-files:
	mkdir -p CARD-files

CARD-files/AT-genes.fa: CARD-files
	wget -O - http://arpcard.mcmaster.ca/blast/db/nucleotide/AT-genes.fa.gz | gunzip -c > CARD-files/AT-genes.fa

CARD-files/AROtags.txt: CARD-files
	wget -O - http://arpcard.mcmaster.ca/blast/db/AROtags.txt.gz | gunzip -c > CARD-files/AROtags.txt
	
CARD-files/aro.obo: CARD-files
	wget -O - http://arpcard.mcmaster.ca/obo-download/aro.obo.gz | gunzip -c > CARD-files/aro.obo

clean:
	rm -r library/
	rm -r taxonomy/
	rm -r CARD-files/
