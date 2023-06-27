setwd("C://Users/mtrw8/OneDrive - The University of Melbourne/workspace/species_assignment")

#read in imm

#create t-SNE
dimm <- readRDS("data/dimm.rds")
td <- tsne::tsne(dimm,perplexity = 30)
td <- readRDS("data/tsne_ibs.Rds") 
