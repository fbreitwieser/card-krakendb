

define result_msg

Congratulations, the database structure was sucessfully generated.

To build the kraken database, run:
     kraken-build --build --db --kmer-len 24 --minimizer-len 12 .

To ensure everything went well, you can test it by assigning the library/ sequences:
    kraken --db . library/* > test.kraken
    kraken-report --db . test.kraken > test.report
endef
export result_msg

all: CARD-files/nucleotide_fasta_protein_homolog_model.fasta CARD-files/aro.csv CARD-files/aro.obo taxonomy/gi_taxid_nucl.dmp library
	Rscript generate-card-db.R
	@echo "$$result_msg"

test_msg:
	@echo "$$result_msg"

library:
	mkdir -p library

taxonomy:
	mkdir -p taxonomy

taxonomy/gi_taxid_nucl.dmp: taxonomy
	printf "1\t1\n" > taxonomy/gi_taxid_nucl.dmp

CARD-files:
	mkdir -p CARD-files
	cd CARD-files
	wget https://card.mcmaster.ca/download/0/broadstreet-v1.1.3.tar.gz
	tar xvvf broadstreet-v1.1.3.tar.gz 
	cd ..

CARD-files/nucleotide_fasta_protein_homolog_model.fasta: CARD-files

CARD-files/AROtags.txt: CARD-files
	
CARD-files/aro.obo: CARD-files

CARD-files/aro.csv: CARD-files

clean:
	rm -r library/
	rm -r taxonomy/
	rm -r CARD-files/
