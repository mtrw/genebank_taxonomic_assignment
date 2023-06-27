

setwd("C://Users/mtrw8/OneDrive - The University of Melbourne/workspace/species_assignment")
source("scripts/0_sp_assign_setup_params.R")

##############################################################################
######################### Curate distance matrix #############################
##############################################################################
imm <- readRDS("data/G2P_pepper_ibs_raw.rds") #read in raw data

#Kill columns/rows containing NA entries
has_na <- apply(imm,1,function(r){any(is.na(r))})
imm <- imm[!has_na,!has_na]

#Attach raw taxonomy information from passport data
spinfo <- fread("data/passport_data_incl_taxassignment.csv",select=c(2,4,13,14),col.names = c("g2p_id","id","genus","species"))
#Curate and remove taxonomic assignments with fewer than 10 representatives
spinfo[,species := paste0(genus,"_",species)][,genus:=NULL]
spinfo[,.N,by=species][grepl("^_$|sp|/",species)]$species -> unassign
spinfo[,.N,by=species][N<10]$species -> kill
spinfo[species %in% unassign,species:="Undetermined"]
spinfo <- spinfo[!species %in% kill]
rcnames <- colnames(imm) #row and column names
rcnames <- gsub("[^[:alnum:] ]", "", rcnames)
spinfo[,id:=gsub("[^[:alnum:] ]", "", id)]

#Records individuals marked as duplicates
dupdat <- fread("data/passport_data_incl_duplicationclustering.csv",select=c(2,8),col.names=c("dupclust","g2p_code"))[!is.na(dupclust)]

#Remove non-G2P samples and those flagged as duplicates
rcnames[rcnames %in% spinfo$g2p_id] #keep
rcnames[!rcnames %in% dupdat[duplicated(dupclust)]$g2p_code] #keep
rcnames[rcnames %in% spinfo$g2p_id & !rcnames %in% dupdat[duplicated(dupclust)]$g2p_code]
imm <- imm[rcnames %in% spinfo$g2p_id & !rcnames %in% dupdat[duplicated(dupclust)]$g2p_code,rcnames %in% spinfo$g2p_id & !rcnames %in% dupdat[duplicated(dupclust)]$g2p_code]
rcnames <- rcnames[rcnames %in% spinfo$g2p_id & !rcnames %in% dupdat[duplicated(dupclust)]$g2p_code]

#Sanity checks
all(colnames(imm) == rcnames) #true or bad
all(rownames(imm) == rcnames) #true or bad
imm[rcnames=="NA",rcnames=="NA"] #empty or it's a bad
all(!is.na(imm)) #true or bad

#Convert taxon names to ranks and create a conversion table for use later
spids_converter <- data.table(
  g2p_id = rcnames
)
spids_converter <- spinfo[spids_converter,on=.(g2p_id)]
spids_converter[,sp_id:=frank(species,ties.method = "dense")]
colnames(imm) <- rownames(imm) <- spids_converter$sp_id

spids <- spids_converter[,.(n_sp=.N),by=.(species,sp_id)]

dimm <- as.dist(imm)
saveRDS(dimm,"data/dimm.rds")

#do not let "undetermined" act like an actual classification
imm[,spids_converter$species=="Undetermined"] <- Inf
#don't let remaining duplicates or selves count
imm[imm==0] <- Inf
diff(dim(imm))==0 #TRUE or it is bad

saveRDS(imm,"data/imm.rds")
saveRDS(spids,"data/spids.rds")
saveRDS(spids_converter,"data/spids_converter.rds")






