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
aro_tags_file <- "CARD-files/AROtags.txt"
at_genes_file <- "CARD-files/AT-genes.fa"

root_id <- "ARO_1000001"

assert <- function(TEST) {
  if (!isTRUE(TEST)) stop("Test failed.")
}

aro <- ontoCAT::getOntology(aro_obo_file)
root_term <- ontoCAT::getTermById(aro, root_id)

##############################################################################################
## generate TAXONOMY
message("Generating taxonomy ...")
dir.create("taxonomy")

## write names.dmp
all_children <- ontoCAT::getAllTermChildren(aro, root_term)
all_labels <- c(ontoCAT::getLabel(root_term), sapply(all_children,ontoCAT::getLabel))
all_acs <- c(ontoCAT::getAccession(root_term), sapply(all_children,ontoCAT::getAccession))
assert(all(table(all_labels)==1))   ## ensure all labels are unique
assert(length(all_labels)==length(all_acs)) 

label_to_id <- setNames(seq_along(all_labels),all_labels)
ac_to_id <- setNames(seq_along(all_acs),all_acs)
names_dmp <- cbind(label_to_id[all_labels],all_labels,"","scientific name\t|")

write.table(names_dmp,
            file="taxonomy/names.dmp",
            quote=FALSE, sep = "\t|\t",
            col.names = FALSE, row.names = FALSE)

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

nodes_dmp[,3] <- paste0("R",nodes_dmp[,3],"\t|")
write.table(nodes_dmp,
            file="taxonomy/nodes.dmp",
            quote=FALSE, sep = "\t|\t",
            col.names = FALSE, row.names = FALSE)

##############################################################################################
## update LIBRARY sequences
message("Updating library sequences ...")

dir.create("library")

tags <- read.delim(aro_tags_file, stringsAsFactors = FALSE, header = FALSE)
head(tags)

tags[,1] <- sub("\\.p0?[123]$","",tags[,1])

## Some tags are duplicated, remove the second ones
table(table(tags[,1]))
duplicated_tags <- duplicated(tags[,1])
message("Removing ",sum(duplicated_tags), "duplicated tags")
tags <- tags[!duplicated_tags,]

## Reads the whole FASTA at once - not advisable for big files
AT_genes <- readLines(at_genes_file)
header <- grepl("^>",AT_genes)
seq_ids <- sub("^>([^ ]*) .*","\\1",AT_genes[header])

## For some sequence ids, there's no entry in AROtags. let's take the tags from the headers instead
seq_ids[!seq_ids %in% tags[,1]]

## Creating a seq-to-ac map. Note that only the first non-root ARO tag is taken!
seq_to_ac <-sapply(strsplit(AT_genes[header]," "),
                   function(x) x[grepl("^ARO",x) & x != root_id][1] )
assert(is(seq_to_ac,"character"))

## Mapping AC to taxId
assert(all(sub(":","_",seq_to_ac) %in% names(ac_to_id)))
seq_taxids <- as.numeric(ac_to_id[sub(":","_",seq_to_ac)])

AT_genes[header] <- paste0(">kraken:taxid|",seq_taxids," ",substring(AT_genes[header],2))
writeLines(AT_genes, "library/AT-genes-updated.fa")


