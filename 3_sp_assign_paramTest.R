setwd("C://Users/mtrw8/OneDrive - The University of Melbourne/workspace/species_assignment")
source("scripts/0_sp_assign_setup_params.R")

##############################################################################
######################### Curate distance matrix #############################
##############################################################################

imm <-  readRDS("data/imm.rds") #read in raw data
# Species ID converter
spids <- readRDS("data/spids.rds")
# Sample list also used to convert species
spids_converter <- readRDS("data/spids_converter.rds")
spids_converter[,sp_id_orig:=sp_id]

#do not let "undetermined" act like an actual classification
imm[,spids_converter$species=="Undetermined"] <- Inf

# For dev purposes
imo <- copy(imm)

# Get tSNE coords
td <- readRDS("data/tsne_ibs.Rds")

# Colour tables
colourtable <- readRDS("data/colourTable.rds")
colourtable <- spids[colourtable[,.(col,species=sp)],on="species"]

optimalGuesses <- readRDS("data/optimal_assignments.rds")

#par(mfrow=c(5,3))

##############################################################################
######### Run analysis in loop over parameters to compare outputs ############
##############################################################################


# # KNN settings in publication
# knn <- 6
# distthres <- Inf
# decay <- 0
# unassignthres <- 4/6
# weight_by_freq = F
# 
# # In publication
knn <- Inf
distthres <- 0.2
decay <- 320
unassignthres <- 0.666 #0
weight_by_freq = F

decays <- seq(0,600,l=10)#0
distthress <- seq(0.01,0.03,l=9)#Inf
knns <- Inf#c(3,5,6,7,10,20,100)
unassignthress <- c(0,1/3,2/3)
weight_by_freqs <- F#F

decays <- 0
distthress <- Inf
knns <- c(2,3,5,6,7,10,20,100)
unassignthress <- c(0,1/3,2/3)
weight_by_freqs <- F

(nPanels <- (length(knns) * length(unassignthress) * length(weight_by_freqs) * length(distthress)) + 1)/2

# Set up PDF for output
fn <- 0
while(file.exists(fName<-paste0("plots/kernelclassifier_param_optimisation_knn",fn<-fn+1,".pdf"))){}
fName
dev.off()
pdf(file=fName,w=8,h=(nPanels/2)*4)
par(mfrow=c(nPanels/2,2))


for( decay in decays ){ # parameter for score decay w/ distance (lambda)
  for( distthres in  distthress ){ # neighbour distance cutoff (d)
  # plotx <- seq(0,1,l=5000)
  # ploty <- exp(-(decay*(plotx+1)))
  # ploty[plotx>distthres] <- 0
  # ploty <- ploty %>% scale_between(0,max(density(imm)$y))
  # pd(imm,main=paste0("Distribution of IBS values with relative\nscaling factor for weighting (blue)\nlambda=",decay," d=",distthres))
  # lines(plotx,ploty,col="#1111FFCC")
    for(weight_by_freq in weight_by_freqs){
      for(knn in knns){ # k nearest neighbours, set k
        for(unassignthres in unassignthress){ # param r; max score must make up proportionally this much or individual will be marked unassigned
          {
            time <- Sys.time()
            #debugonce(run)
            run(
              knn = knn,
              distthres = distthres,
              decay = decay,
              unassignthres = unassignthres,
              weight_by_freq = weight_by_freq,
              imm=copy(imm),
              spids_converter=copy(spids_converter),
              spids=copy(spids),
              tsneDat = copy(td),
              colourtable=copy(colourtable),
              plot = T
            )
          }
          ce("That took ",Sys.time()-time,"s")
          
        #   run_heatmap(
        #     knn = knn,
        #     distthres = distthres,
        #     decay = decay,
        #     unassignthres = unassignthres,
        #     weight_by_freq = weight_by_freq,
        #     imm=copy(imm),
        #     spids_converter=copy(spids_converter),
        #     spids=copy(spids),
        #     tsneDat = copy(td),
        #     colourtable=copy(colourtable))
        # }
        }
      }
    }
  }
}
#
dev.off()

#saveRDS(assign[,.(g2p_code,sp_id_orig=sp_id,sp_id_optimal=sp_id_new)],"data/optimal_assignments.rds")

