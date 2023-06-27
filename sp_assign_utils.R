
#useful wrapper for kernel density plots
pd <- function(x,add=F,main=NA,...){
  if(!add){
    x %>% density(na.rm=TRUE,...) %>% plot(main=main)
  } else {
    x %>% density(na.rm=TRUE,...) %>% lines(main=main)
  }
}



run <- function(knn = Inf, distthres = 0.2, decay = 320, unassignthres = 0.66, weight_by_freq = T,imm,spids_converter,spids,tsneDat,colourtable,plot=F){
  
  ce(paste0("lambda=",decay," d=",distthres," r=",unassignthres," knn=",knn))
  
  #Reset and set up score matrix d
  gpids <- spids$sp_id
  ngps <- length(gpids)
  nsamps <- nrow(imm)
  d <- exp(-(decay*(imm+1))) #score based on distances
  
  d[imm>distthres] <- 0 #no score contribution for any pair beyond this far apart
  #k nearest neighbours
  if(knn < Inf){
    l_ply( 1:nrow(imm),function(i) { d[i,][order(imm[i,],decreasing=T)[1:(ncol(imm)-knn)]] <<- 0 } )
  }
  d[1:10,1:10]
  imm[1:10,1:10]
  
  
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
        #browser()
        as.numeric(sum(i[s]<distthres))
      })
    })
    
    weightsPerIndGp <- t(apply(weightsPerIndGp,1,function(r){
      r/sum(r)
    }))
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
  assign[,sp_id_new := apply(scores,1,function(r){
    #browser()
    s <- sort(r,d=T)
    if(any(is.nan(r))){
      NA
    } else if(s[1]!=s[2]){ # if((max(c)/sum(c)) > 0.5){
      which(r==max(r))[1]
    } else {
      NA
    }
  })]
  assign[,max_prop := apply(scores,1,function(r) {max(r)})]
  
  #join taxonomic info to assignment calls
  assign <- spids[,.(species,sp_id)][assign,on=.(sp_id)]
  assign <- spids[,.(species_new=species,sp_id_new=sp_id)][assign,on=.(sp_id_new)]
  
  #apply cutoff param
  assign[max_prop < unassignthres,species_new := "Undetermined"]
  
  # NAs occur when no classified points were within the dist thres or best guesses were ties
  assign[is.na(species_new),species_new := "Undetermined" ]
  assign[is.na(species_new),sp_id_new := spids[species=="Undetermined"]$sp_id ]
  
  tpdat <- cbind(assign,data.table(X1=tsneDat[,1],X2=tsneDat[,2]))
  tpdat[,cex:=0.8]
  tpdat[,col:=replace_levels_with_colours(species,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
  tpdat[,col_new:=replace_levels_with_colours(species_new,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
  tpdat[species=="Undetermined",col:="#0000007D"]
  tpdat[species_new=="Undetermined",col:="#0000007D"]
  
  tpdat[,col:=NULL]
  #tpdat1[,col:=NULL]
  
  tpdat <- colourtable[tpdat,on="species"]
  #tpdat1 <- colourtable[,.(col_new=col,species_new=sp)][tpdat1,on="species_new"]
  tpdat <- colourtable[,.(col_new=col,species_new=species)][tpdat,on="species_new"]
  
  
  
  shuffler <- data.table(g2p_code=tpdat$g2p_code)
  shuffler <- shuffler[sample(1:.N)]
  tpdat <- shuffler[tpdat,on=.(g2p_code)]
  #tpdat1 <- shuffler[tpdat1,on=.(g2p_code)]
  
  if(plot==T) {plot(tpdat$X1,tpdat$X2,col=tpdat$col_new,pch=20,cex=tpdat$cex,main=paste0("lambda=",decay," d=",distthres," r=",unassignthres," knn=",knn),xlab="X_1",ylab="X_2")}
  
  #plot(assign)
  
  
  notNaPos <- (!is.na(optimalGuesses$sp_id_optimal)) & (!assign[,species_new=="Undetermined"])
  
  data.table(
    knn,
    distthres,
    decay,
    unassignthres,
    weight_by_freq,
    pr_wrong=sum(optimalGuesses[notNaPos]$sp_id_optimal!=assign[notNaPos]$sp_id_new)/sum(notNaPos),
    n_additional_na=sum(!notNaPos)-sum(is.na(optimalGuesses$sp_id_optimal))
  ) %>% return()
}








run_heatmap <- function(knn = Inf, distthres = 0.2, decay = 320, unassignthres = 0.66, weight_by_freq = T,imm,spids_converter,spids,tsneDat,colourtable){
  
  ce(paste0("lambda=",decay," d=",distthres," r=",unassignthres," wbf=",weight_by_freq," knn=",knn))
  
  #Reset and set up score matrix d
  gpids <- spids$sp_id
  ngps <- length(gpids)
  nsamps <- nrow(imm)
  d <- exp(-(decay*(imm+1))) #score based on distances
  
  d[imm>distthres] <- 0 #no score contribution for any pair beyond this far apart
  #k nearest neighbours
  if(knn < Inf){
    l_ply( 1:nrow(imm),function(i) { d[i,][order(imm[i,],decreasing=T)[1:(ncol(imm)-knn)]] <<- 0 } )
  }
  d[1:10,1:10]
  imm[1:10,1:10]
  
  
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
        #browser()
        as.numeric(sum(i[s]<distthres))
      })
    })
    
    weightsPerIndGp <- t(apply(weightsPerIndGp,1,function(r){
      r/sum(r)
    }))
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
  
  #join taxonomic info to assignment calls
  assign <- spids[,.(species,sp_id)][assign,on=.(sp_id)]
  assign <- spids[,.(species_new=species,sp_id_new=sp_id)][assign,on=.(sp_id_new)]
  
  #apply cutoff param
  assign[max_prop < unassignthres,species_new := "Undetermined"]
  
  # NAs occur when no classified points were within the dist thres
  assign[is.na(species_new),species_new := "Undetermined" ]
  assign[is.na(species_new),sp_id_new := spids[species=="Undetermined"]$sp_id ]
  
  tpdat <- cbind(assign,data.table(X1=tsneDat[,1],X2=tsneDat[,2]))
  tpdat[,cex:=0.8]
  tpdat[,col:=replace_levels_with_colours(species,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
  tpdat[,col_new:=replace_levels_with_colours(species_new,fun="qualitative_hcl",palette="Dark 3",alpha=.5,plot=F,newplot = T)]
  tpdat[species=="Undetermined",col:="#0000007D"]
  tpdat[species_new=="Undetermined",col:="#0000007D"]
  
  tpdat[,col:=NULL]
  #tpdat1[,col:=NULL]
  
  tpdat <- colourtable[tpdat,on="species"]
  #tpdat1 <- colourtable[,.(col_new=col,species_new=sp)][tpdat1,on="species_new"]
  tpdat <- colourtable[,.(col_new=col,species_new=species)][tpdat,on="species_new"]
  
  
  
  shuffler <- data.table(g2p_code=tpdat$g2p_code)
  shuffler <- shuffler[sample(1:.N)]
  tpdat <- shuffler[tpdat,on=.(g2p_code)]
  #tpdat1 <- shuffler[tpdat1,on=.(g2p_code)]
  
  plot(tpdat$X1,tpdat$X2,col=tpdat$col,pch=20,cex=tpdat$cex,main=paste0("lambda=",decay," d=",distthres," r=",unassignthres," wbf=",weight_by_freq," knn=",knn),xlab="X_1",ylab="X_2")
  
  plot(tpdat$X1,tpdat$X2,col=tpdat$col_new,pch=20,cex=tpdat$cex,main=paste0("lambda=",decay," d=",distthres," r=",unassignthres," wbf=",weight_by_freq," knn=",knn),xlab="X_1",ylab="X_2")
  
  #plot(assign)
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
  
  browser()
  #heatmap
  colnames(changemat) <- paste0("... to ",sub("_"," ",colnames(changemat)))
  rownames(changemat) <- paste0("From ",sub("_"," ",rownames(changemat))," ...")
  pheatmap(mat = changemat,cluster_rows = F,cluster_cols = F,display_numbers = apply(changemat,2,function(c) paste0(round(c,digits = 2),"%")))
  
  
}
