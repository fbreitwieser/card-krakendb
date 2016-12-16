#!/bin/env Rscript

if (!"ontoCAT" %in% rownames(installed.packages())) {
    message("Installing ontoCAT from Bioconductor")
    source("http://bioconductor.org/biocLite.R")
    tryCatch ({
        biocLite("ontoCAT")
        },error = function(e) {
            message("Installation failed. Make sure you have Java installed, ",
                    " then type 'R CMD javareconf -e' and restart.\n",
                    "If it stills fails, make sure liblzma-dev is installed.\n")
        })
}

library(ontoCAT)

aro_obo_file <- "CARD-files/aro.obo"
aro_tags_file <- "CARD-files/aro.csv"

root_id <- "ARO_1000001"

assert <- function(TEST) {
  if (!isTRUE(TEST)) stop("Test failed.")
}

aro <- ontoCAT::getOntology(aro_obo_file)
root_term <- ontoCAT::getTermById(aro, root_id)
all_children <- ontoCAT::getAllTermChildren(aro, root_term)
all_acs <- c(ontoCAT::getAccession(root_term), sapply(all_children,ontoCAT::getAccession))
aro_ac_to_id <- setNames(seq_along(all_acs),all_acs)

##############################################################################################
## generate TAXONOMY
message("Generating taxonomy ...")
if (!dir.exists("taxonomy")) dir.create("taxonomy")

## write names.dmp
  all_labels <- c(ontoCAT::getLabel(root_term), sapply(all_children,ontoCAT::getLabel))
  assert(all(table(all_labels)==1))   ## ensure all labels are unique
  assert(length(all_labels)==length(all_acs)) 
  
  #label_to_id <- setNames(seq_along(all_labels),paste0(all_labels,"(",all_acs,")"))
  label_to_id <- setNames(seq_along(all_labels),all_labels)
  
  ####### write nodes.dmp

## difficulty: some terms have more than one parent
## for example

# [Term]
# id: ARO:3000857
# name: 16S ribosomal RNA methyltransferase
# namespace: antibiotic_resistance
# def: "Methyltransferases that modify the 16S rRNA of the 30S subunit of bacterial ribosomes, conferring # resistance to drugs that target 16S rRNA." [PMID:21177880, PMID:14667745, PMID:22547620]
# is_a: ARO:3000104 ! aminoglycoside resistance gene
# is_a: ARO:3000164 ! rRNA methyltransferase conferring antibiotic resistance

## solutions:
## 1) open aro.obo in vim before loading, find all the terms with more than one is_a (/is_a.*\nis_a)
##    and delete the ones you don't like
## 2) in the code here, just the first time a child appears it is taken

  ## matrix columns: child, parent, depth
  nodes_dmp <- matrix(NA,ncol=3,nrow=length(all_labels))
  nodes_dmp[1,] <- c(label_to_id[1],label_to_id[1],0)
  current_row <- 2
  
  get_id <- function(label) {
    if (!label %in% names(label_to_id)) {
      stop("Cannot find label ",label)
    }
    label_to_id[label]
  }
  
  recurse_to_leafs <- function(aro, parent_term, depth) {
    children <- ontoCAT::getTermChildren(aro, parent_term)
    if (length(children) > 0) {
      for (child in children) {
        child_id <- recurse_to_leafs(aro, child, depth + 1)
      
        # filter out children that have already occured
        if (!child_id %in% nodes_dmp[seq(from=1, to=current_row-1),1]) {
          nodes_dmp[current_row, ] <<- c(child_id, 
                                         get_id(ontoCAT::getLabel(parent_term)),
                                         depth)
          current_row <<- current_row + 1
        }
      }
    }
    return(get_id(ontoCAT::getLabel(parent_term)));
  }
  ret <- recurse_to_leafs(aro,root_term, 1)
  assert(as.numeric(ret) == 1)
  assert(all(table(nodes_dmp[,1]) == 1))
  

##############################################################################################
## update LIBRARY sequences

  
message("Updating library sequences ...")
if (!dir.exists("library")) dir.create("library")

#library(RJSONIO)  
#aro_data <- read.csv("CARD-files/aro.csv")
#card_data <- fromJSON("CARD-files/card.json")

tags <- read.delim(aro_tags_file, stringsAsFactors = FALSE, header = FALSE)
tags[,1] <- sub("\\.p0?[123]$","",tags[,1])

## Some tags are duplicated, remove the second ones
duplicated_tags <- duplicated(tags[,1])
message("Removing ",sum(duplicated_tags), " duplicated tags")
tags <- tags[!duplicated_tags,]

## Reads the whole FASTA at once - not advisable for big files
#fasta_files <- Sys.glob("CARD-files/nucleotide_fasta_protein_*.fasta")
fasta_files <- Sys.glob("CARD-files/nucleotide_fasta_protein_homolog_model.fasta")
AT_genes <- do.call(c, lapply(fasta_files,readLines))

header <- grepl("^>",AT_genes)
## get the names from the headers, which we'll save in names.dmp
header_sentences <- AT_genes[header]

## use gene_names for one level in the hierachy
gene_names <- sapply(header_sentences, function(x) x[1])

seq_names <- sub("\\|ARO:[0-9]*", "", header_sentences)
full_gene_names <- sub(".*\\|", "", header_sentences)
gene_names <- sub(" .*", "", full_gene_names)
aro_tags <- sub(".*\\|(ARO:[0-9]*).*", "\\1", header_sentences)

assert(all(!duplicated(seq_names)))

## Creating a seq-to-ac map. Note that only the first non-root ARO tag is taken!
seq_to_aro_ac <- aro_tags

## Mapping the sequemces to a parent taxId (from ARO)
assert(all(sub(":","_",seq_to_aro_ac) %in% names(aro_ac_to_id)))
seq_parent_taxids <- as.numeric(aro_ac_to_id[sub(":","_",seq_to_aro_ac)])

## Add a new tax-ID for all the individual sequences to names.dmp
## 10000 is arbitrarily chosen to be above all others
gene_taxids <- setNames(seq(from=10000,length.out = length(unique(gene_names))),unique(gene_names))

seq_taxids <- seq(from=20000,length.out = sum(header))

label_to_ac <- setNames(sub("_",":",names(aro_ac_to_id)),seq_names)

#TODO: update all_labels to include AC
all_labels_ext <- paste0(all_labels," (",sub("_",":",all_acs),")")

names_dmp <- cbind(as.character(c(label_to_id[all_labels],
								  gene_taxids,
								  seq_taxids)),
                   as.character(c(all_labels_ext,
								  names(gene_taxids),
								  seq_names)),
                   "","scientific name\t|")

write.table(names_dmp, 
            file="taxonomy/names.dmp",
            quote=FALSE, sep = "\t|\t",
            col.names = FALSE, row.names = FALSE)

AT_genes[header] <- paste0(">kraken:taxid|",seq_taxids," ",substring(AT_genes[header],2))
writeLines(AT_genes, "library/lib-sequences.fa")


## write connections from parent_taxid to sequence taxids into nodes_dmp
rownames(nodes_dmp) <- nodes_dmp[,1]

nodes_dmp_genes <- unique(cbind(gene_taxids[gene_names],
								seq_parent_taxids,
								nodes_dmp[seq_parent_taxids]+1))
nodes_dmp_genes <- nodes_dmp_genes[!duplicated(nodes_dmp_genes[,1]),]


nodes_dmp_seq <- cbind(seq_taxids,
					   gene_taxids[gene_names],
					   nodes_dmp[seq_parent_taxids]+2)

nodes_dmp <- rbind(nodes_dmp, nodes_dmp_genes, nodes_dmp_seq)

assert(sum(duplicated(nodes_dmp[,1]))==0)

nodes_dmp[,3] <- paste0("R",nodes_dmp[,3],"\t|")
assert(!any(seq_taxids %in% label_to_id))
new_leaf_ids <- seq(from=max(label_to_id) + 1,length.out = nrow(names_dmp) - 1)



## Append nodes.dmp and names.dmp with names and links from the file
write.table(nodes_dmp, 
            file="taxonomy/nodes.dmp",
            quote=FALSE, sep = "\t|\t",
            col.names = FALSE, row.names = FALSE)


