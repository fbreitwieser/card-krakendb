
Scripts to generate a Kraken database from the 
The Comprehensive Antibiotic Resistance Database
-------------------------------------------------

These scripts 

1) Download the antibiotic resistance genes and ontology from the Comprehensive Antibiotic Resistance Database

2) Create a taxonomy based on the ontology file

3) Update the sequences such that the new taxonomy IDs 

Note that when either a ontology term has multiple parents, or a sequence is assigned to more than one term,
one is arbitrarily taken.

The repository already contains the processed library/ and taxonomy/ files - so you can go ahead and run 

    kraken-build --build --db .

If you would like to start from scratch, run

   make clean
   make
