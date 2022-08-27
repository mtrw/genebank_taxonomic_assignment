#
##############################################################################
######################### Set up environment #################################
##############################################################################

setwd("WORKING DIR GOES HERE; MUST INCLUDE A SUBDIR `DATA` THAT INCLUDES ALL THE DATA FILES ASSOCIATED WITH THE PAPER (DOI:10.5281/zenodo.7016070)")

source("https://raw.githubusercontent.com/mtrw/tim_r_functions/master/tim_functions.R")
require(pheatmap)
require(tsne)

#useful wrapper for kernel density plots
pd <- function(x,add=F,main=NA,...){
  if(!add){
    x %>% density(na.rm=TRUE,...) %>% plot(main=main)
  } else {
    x %>% density(na.rm=TRUE,...) %>% lines(main=main)
  }
}

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

#Convert taxon names to ranks and create a conversiton table for use later
spids_converter <- data.table(
  g2p_id = rcnames
)
spids_converter <- spinfo[spids_converter,on=.(g2p_id)]
spids_converter[,sp_id:=frank(species,ties.method = "dense")]
colnames(imm) <- rownames(imm) <- spids_converter$sp_id

spids <- spids_converter[,.(n_sp=.N),by=.(species,sp_id)]

#Original input matrix, before any manipulation of values, just useful to have around for debugging
imo <- copy(imm) #imm <- imo

#do not let "undetermined" act like an actual classification
imm[,spids_converter$species=="Undetermined"] <- Inf

#don't let remaining duplicates or selves count
imm[imm==0] <- Inf
dim(imm)

#create t-SNE
#require(tsne)
#imod <- as.dist(imo)
#td <- tsne::tsne(imod,perplexity = 30)
################################################################
td <- readRDS("data/tsne_ibs.Rds") #DELETE LINE BEFORE PUB



#Set up PDF for output
dev.off()
pdf(file="plots/kernelclassifier_param_optimisation.pdf",w=12,h=21)
par(mfrow=c(5,3))


##############################################################################
######### Run analysis in loop over parameters to compare outputs ############
##############################################################################
weight_by_freq <- T


distthres=0.2; decay=320; unassignthres=0.66; knn=Inf

#Reset and set up score matrix d
gpids <- spids$sp_id
ngps <- length(gpids)
nsamps <- nrow(imm)
d <- exp(-(decay*(imm+1))) #score based on distances

d[imm>distthres] <- 0 #no score contribution for non-neighbours
#k nearest neighbours
l_ply(1:nrow(imm),function(i) { d[i,][order(imm[i,],decreasing=T)[1:(length(imm[i,])-min(knn,length(imm[i,])))]] <<- 0 } )

#plot distribution of ibs scores and weighting scheme
plotx <- seq(0,1,l=5000)
ploty <- exp(-(decay*(plotx+1)))
ploty[plotx>distthres] <- 0
ploty <- ploty %>% scale_between(0,max(density(imm)$y))
pd(imm,main=paste0("Distribution of IBS values with relative\nscaling factor for weighting (blue)\nlambda=",decay," d=",distthres," r=",unassignthres," wbf=",weight_by_freq," knn=",knn))
lines(plotx,ploty,col="#1111FFCC")

#design matrix: which inds are in what tax group
designmat <- matrix(0,nrow=nsamps,ncol=ngps)
l_ply(1:nsamps,function(i){
  designmat[i,spids_converter$sp_id[i]] <<- 1
})

#count non-NA inds within range of each ind from each group, so we can normalise for it
if(weight_by_freq==T){ #we're using my method, weight the neighbours by (the inverse of) their frequency in the neighbour pool
  weightsPerIndGp <- apply(designmat==1,2,function(s){
    #selector for species s --- run on each col of design
    apply(imm,1,function(i){
      as.numeric(sum(i[s]<distthres))
    })
  })
  #convert to normalisation factors
  weightsPerIndGp[weightsPerIndGp!=0.0] <- 1.0/weightsPerIndGp[weightsPerIndGp!=0.0]
  #sum the scores and perform the normalisation
  scores <- (d %*% designmat) * weightsPerIndGp
} else {
  scores <- (d %*% designmat)
}
#add tax id info in nice printable format
rownames(scores) <- paste0("ind:",1:nsamps," gp:",spids_converter$sp_id)
colnames(scores) <- paste0("dist_to_gp_",1:ncol(scores))
#convert scores to proportions
scores <- t(apply(scores,1,function(r){r/sum(r)}))

#convert to assignment table and make assignments
assign <- as.data.table(scores)
assign[,g2p_code:=spids_converter$g2p_id]
assign[,sp_id:=spids_converter$sp_id]
assign[,description:=rownames(scores)]
assign[,sp_id_new := apply(scores,1,function(r) {which(r==max(r))[1]})]
assign[,max_prop := apply(scores,1,function(r) {max(r)})]

#plot to show effect of cutoff param r
pd(assign$max_prop,main=paste0("Distribution of dominant proportional\nassignments with undetermined\ncutoff (blue; ",unassignthres,")"))
abline(v=unassignthres,col="blue")

#join taxonomic info to assignment calls
assign <- spids[,.(species,sp_id)][assign,on=.(sp_id)]
assign <- spids[,.(species_new=species,sp_id_new=sp_id)][assign,on=.(sp_id_new)]

#apply cutoff param
assign[max_prop < unassignthres,species_new := "Undetermined"]

#add useful columns for printing and plotting
assign[,n_sp:=.N,,by=.(species)]
assign[,n_sp_new:=.N,,by=.(species_new)]

#wrangle taxon names
spnamer <- assign[,.(n_sp=.N),,by=.(species)]
spnamerfill <- spnamer[,.(species_new=species)]
spnamernew <- assign[,.(n_sp_new=.N),by=.(species_new)][spnamerfill,on="species_new"]
spnamernew[is.na(n_sp_new),n_sp_new:=0]
spnamer[,sp_name:=paste0(species," (",n_sp,")")]
spnamernew[,sp_name_new:=paste0(species_new," (",n_sp_new,")")]

#create summary of changes (best for a heatmap)
changesumm <- assign[,.(n_sp_sp_new=.N),by=.(species,species_new,n_sp,n_sp_new)]
changesumm <- changesumm[,pr_sp_sp_new:=n_sp_sp_new/n_sp][]
filler <- setDT(expand.grid(species=unique(changesumm$species),species_new=unique(changesumm$species)))
changesumm <- changesumm[,.(species,species_new,pr_sp_sp_new,n_sp_sp_new)][filler,on=.(species,species_new)]
changesumm[is.na(pr_sp_sp_new),pr_sp_sp_new:=0.0]
changesumm[is.na(n_sp_sp_new),n_sp_sp_new:=0]
rm(filler)
t <- dcast(changesumm,species~species_new,value.var="pr_sp_sp_new")
changemat <- as.matrix(t[,2:ncol(t)]) * 100
rm(t)
nameconv <- data.table(species=colnames(changemat),species_new=colnames(changemat))
rownames(changemat) <- spnamer[nameconv,on="species"]$sp_name
colnames(changemat) <- spnamernew[nameconv,on="species_new"]$sp_name_new

#heatmap
colnames(changemat) <- paste0("... to ",sub("_"," ",colnames(changemat)))
rownames(changemat) <- paste0("From ",sub("_"," ",rownames(changemat))," ...")
#pheatmap(mat = changemat,cluster_rows = F,cluster_cols = F)

#accumulate plot data
tpdat <- cbind(assign,data.table(X1=td[,1],X2=td[,2]))
tpdat[,cex:=0.8]
tpdat[,col:=replace_levels_with_colours(species,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
tpdat[,col_new:=replace_levels_with_colours(species_new,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
tpdat[species=="Undetermined",col:="#0000007D"]
tpdat[species_new=="Undetermined",col:="#0000007D"]

#Cluster data from the original Tripodi et al paper, for plotting in comparison if required
clusts <- fread("data/passport_data_incl_assignmentclustering.csv",select=2:4,col.names=c("g2p_code","species","clust"))
clusts[,species:=swap(
  species,
  c(
    "C.baccatum_pend",
    "C.baccatum",
    "C.pubescens",
    "C.chacoense",
    "C.baccatum_bac",
    "C.annuum",
    "C.praetermissum",
    "C.eximium",
    "C.frutescens",
    "C.chinense",
    "Undefined",
    "C.annuum_glab",
    "C.cardenasii",
    "C.tovarii",
    "C.galapagoense"
  ),
  c(
    "Capsicum_baccatum",
    "Capsicum_baccatum",
    "Capsicum_pubescens",
    "Capsicum_chacoense",
    "Capsicum_baccatum",
    "Capsicum_annuum",
    "Kill",
    "Capsicum_eximium",
    "Capsicum_frutescens",
    "Capsicum_chinense",
    "Undetermined",
    "Capsicum_annuum",
    "Kill",
    "Kill",
    "Capsicum_galapagoense"
  )
)]
clusts <- clusts[species!="Kill"]

clustsumm <- clusts[,n_sp_in_clust:=.N,by=.(clust,species)]
clustsumm[,n_clust:=.N,by=.(clust)]
clustsumm[,pr_clust:=n_sp_in_clust/n_clust]
setkey(clustsumm,clust,n_sp_in_clust)
setorder(clustsumm,clust,-n_sp_in_clust)

clustsumm[,species_new:={
  if(pr_clust[1]>0.8){
    species[1]
  } else {
    species
  }
},by=clust]


clustsumm[,n_sp:=.N,,by=.(species)]
clustsumm[,n_sp_new:=.N,,by=.(species_new)]

spnamer <- clustsumm[,.(n_sp=.N),,by=.(species)]
spnamernew <- clustsumm[,.(n_sp_new=.N),,by=.(species_new)]
spnamer[,sp_name:=paste0(species," (",n_sp,")")]
spnamernew[,sp_name_new:=paste0(species_new," (",n_sp_new,")")]

changesumm <- clustsumm[,.(n_sp_sp_new=.N),by=.(species,species_new,n_sp,n_sp_new)]
changesumm <- changesumm[,pr_sp_sp_new:=n_sp_sp_new/n_sp][]

filler <- setDT(expand.grid(species=unique(changesumm$species),species_new=unique(changesumm$species)))
changesumm <- changesumm[,.(species,species_new,pr_sp_sp_new,n_sp_sp_new)][filler,on=.(species,species_new)]
changesumm[is.na(pr_sp_sp_new),pr_sp_sp_new:=0.0]
changesumm[is.na(n_sp_sp_new),n_sp_sp_new:=0]
rm(filler)

t <- dcast(changesumm,species~species_new,value.var="pr_sp_sp_new")
changemat2 <- as.matrix(t[,2:ncol(t)]) * 100
rm(t)

nameconv <- data.table(species=colnames(changemat2),species_new=colnames(changemat2)) #name conversion table


rownames(changemat2) <- spnamer[nameconv,on="species"]$sp_name
colnames(changemat2) <- spnamernew[nameconv,on="species_new"]$sp_name_new

pdat2 <- clustsumm[,.(g2p_code,species,species_new,species_changed=(species!=species_new))]
pdat2[,col:=replace_levels_with_colours(species_new,fun="qualitative_hcl",palette="Dark 3",alpha=.5)]
pdat2[,cex:=0.8]
#pdat2[species_changed==TRUE,cex:=3.0]
tpdat1 <- tpdat[,.(g2p_code,X1,X2)][pdat2,on=.(g2p_code)]
tpdat1$col <- pdat2$col
tpdat1$cex <- pdat2$cex
tpdat1$g2p_code <- pdat2$g2p_code
tpdat1$species <- pdat2$species
tpdat1$species_new <- pdat2$species_new
tpdat1$species_changed <- pdat2$species_changed
pdat2[,col:=replace_levels_with_colours(species_new,fun="qualitative_hcl",palette="Dark 3",alpha=.5)]
colourtable <- pdat2[,.N,by=.(col,sp=species_new)][,N:=NULL][]
colourtable[sp=="Undetermined",col:="#0000007D"]

tpdat[,col:=NULL]
tpdat1[,col:=NULL]

tpdat <- colourtable[,.(col,species=sp)][tpdat,on="species"]
tpdat1 <- colourtable[,.(col_new=col,species_new=sp)][tpdat1,on="species_new"]
tpdat <- colourtable[,.(col_new=col,species_new=sp)][tpdat,on="species_new"]

shuffler <- data.table(g2p_code=intersect(tpdat$g2p_code,tpdat1$g2p_code))
shuffler <- shuffler[sample(1:.N)]
tpdat <- shuffler[tpdat,on=.(g2p_code)]
tpdat1 <- shuffler[tpdat1,on=.(g2p_code)]

# heatmap for Tripodi 2021 reclassification
colnames(changemat2) <- paste0("... to ",sub("_"," ",colnames(changemat2)))
rownames(changemat2) <- paste0("From ",sub("_"," ",rownames(changemat2))," ...")
#pheatmap(mat = changemat2,cluster_rows = F,cluster_cols = F)


  #### Uncomment these depending on which plots you wish to see
#original passport assignments
plot(tpdat$X1,tpdat$X2,col=tpdat$col,pch=20,cex=.8,main="tSNE (perplexity=30) on first three PCs: Original assignment",xlab="X_1",ylab="X_2")
  #uncomment to overlay the legend
#legend(x=10,y=30,legend = colourtable$sp,fill = colourtable$col,cex=.7)

#h-clustering method (tripodi et al)
plot(tpdat1$X1,tpdat1$X2,col=tpdat1$col,pch=20,cex=tpdat1$cex,main="tSNE (perplexity=30) on first three PCs: Cluster-based assignment",xlab="X_1",ylab="X_2")
#legend(x=10,y=30,legend = colourtable$sp,fill = colourtable$col,cex=.7)

#Kernel-based method
#plot(tpdat$X1,tpdat$X2,col=tpdat$col_new,pch=20,cex=tpdat$cex,main="tSNE (perplexity=30) on first three PCs: Kernel-based assignment",xlab="X_1",ylab="X_2")
plot(tpdat$X1,tpdat$X2,col=tpdat$col_new,pch=20,cex=tpdat$cex,main=NA,xlab="X_1",ylab="X_2")
#key
#legend(x=10,y=30,legend = colourtable$sp,fill = stringi::stri_sub(colourtable$col,1,7),cex=.7)

