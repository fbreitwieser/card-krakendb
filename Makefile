

define result_msg
Congratulations, the database structure was sucessfully generated.

To build the kraken database, run:
     kraken-build --build --db .

To ensure everything went well, you can test it by assigning the library/ sequences:
    kraken --db . library/* > test.kraken
    kraken-report --db . > test.report
endef
export result_msg

all: CARD-files/AT-genes.fa CARD-files/AROtags.txt CARD-files/aro.obo taxonomy/gi_taxid_nucl.dmp library
	Rscript generate-card-db.R
	@echo "\n$$result_msg\n"

test_msg:
	@echo "\n$$result_msg\n"

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
