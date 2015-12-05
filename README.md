
Kraken database from the CARD Antibiotic Resistance Database
-------------------------------------------------

These scripts 

1. Download the antibiotic resistance genes and ontology from the Comprehensive Antibiotic Resistance Database
2. Create a taxonomy based on the ontology file
3. Update the sequences such that the new taxonomy IDs 

Note that when either a ontology term has multiple parents, or a sequence is assigned to more than one term,
one is arbitrarily taken.

## Usage

The repository already contains the processed `library/` and `taxonomy/` files - so you can go ahead and run 

    kraken-build --build --db .

If you would like to start from scratch, run

    make clean
    make
    kraken-build --build --db .

## Testing

To test the new database, you can try re-classifying the library sequences

    kraken --db . library/* > test.kraken
    kraken-reprot --db . > test.report
    
(c) Florian Breitwieser, 2015
