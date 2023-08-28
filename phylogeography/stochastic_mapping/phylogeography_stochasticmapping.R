## Loading required libraries
library(foreach)
library(doParallel)
library(phytools)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggtree)
library(plyr)
library(gridExtra)
library(circlize)

## Read in the maximum likelihood global phylogeny of ST34, as well as the metadata
tree <- read.newick(file="./RAxML_bipartitions.raxml_ST34rmblk_454str.japanroot.newick")
metadata <- read.csv(file="./ST34_454str_short_meta.csv", header=TRUE, stringsAsFactors = FALSE)
tree$tip.label <- str_replace_all(tree$tip.label,"'", "")

intree <- match(tree$tip.label, metadata$ID)
geo <- as.vector(metadata[intree, "Region"])
table(geo) ## minium of 32 genomes per geography
names(geo) <- tree$tip.label

## Create a function to subsample from a phylogeny
resample_trees <- function(dat, ML.tree, state, K, iti) { # stands for metadata, tree, the geographical state file, number to retain in each state, number of subsampled tree. 
      require(ape)
      # First subsampling the names to include exact K per character state, to avoid sampling bias
      sampling_func <- function(state, dat, K) {
            sam <- c()
            for (i in levels(factor(state))){sam <- c(sam,as.vector(sample(dat[which(dat$Region==i),"ID"], K, replace=FALSE)))}
            return(sam)
      }
      tak <- replicate(iti, sampling_func(state,dat,K))
      tak <- t(tak)
      ntax <- ncol(tak) # no. of isolates per tree
      fit.trees <- list(as.phylo(rtree(ntax)))
      for (i in 1:iti) {
            inc.tip <- ML.tree$tip.label[which(!(ML.tree$tip.label %in% tak[i,]))]
            ### Subsample the ML tree to create iti new trees with subsampled data
            fit.trees[[i]] <- drop.tip(ML.tree, tip=inc.tip, trim.internal=T)
      }
      return(fit.trees)
}

## creating a set of 1000 trees, each is a subsampled tree of the original ML tree with equal sampling from each geographical region
st34.trees1000 <- resample_trees(dat=metadata, ML.tree=tree, state=geo, K=30, iti=1000) ## select K=30

class(st34.trees1000) <- 'multiPhylo'

### Stochastic mapping using simmap
setCores <- round(detectCores()*0.85) ## use 85% of available logical processors
cl <- parallel::makeCluster(getOption("cl.cores", setCores))
doParallel::registerDoParallel(cl = cl)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

## Function for stochastic mapping with parallel computing
simmap_resample_v4 <- function(trees, states){
      ntax <- length(trees[[1]]$tip.label) # Get the number of taxa for each tree
      fak <- c(rep("ka", ntax/2), rep("ak", ntax/2)) # create a fake matrix first
      names(fak) <- paste("t", seq(1,ntax),sep ="")
      M <- nlevels(factor(states))
      res <- foreach (i = 1:length(trees), .combine='rbind') %dopar% {
            require(phytools)
            require(stringr)
            cat("Running simmap for iteration ",i, "\n")
            trait <- states[which(names(states) %in% trees[[i]]$tip.label)]
            cen.sim <- make.simmap(trees[[i]], trait, model='ARD', nsim=100, Q='mcmc',tol=1e-8, prior="estimated")
            pds <- describe.simmap(cen.sim)
            ## Summarize the state transition from the pds files
            kk <- summary(pds$count)
            colnames(kk) <- colnames(pds$count)
            ll <- str_replace(string=kk[3,], pattern="Median.*:", "")
            ll <- as.numeric(str_replace(string = ll, pattern="  $",""))
            mm <- summary(pds$times)
            colnames(mm) <- colnames(pds$times)
            nn <- str_replace(string=mm[4,], pattern="Mean.*:", "")
            nn <- as.numeric(str_replace(string = nn, pattern="  $",""))
            nn <- nn/nn[length(nn)]
            c(ll,nn)
      }
      return(res)
}

## Running a test simmap on one subtree to get column names of result file
test.geo <- geo[which(names(geo) %in% st34.trees1000[[1]]$tip.label)]
test.simmap <- make.simmap(st34.trees1000[[1]], test.geo, model='ARD', nsim=100, Q='mcmc',tol=1e-8, prior="estimated")
pds <- describe.simmap(test.simmap)
simmap.colnames <- c(colnames(pds$count), colnames(pds$times))

#############################################################
# Main command, stochastic mapping is performed on each subsampled phylogeny
## This will take very longgg 
############################################################

simmap_results <- simmap_resample_v4(st34.trees1000, geo)

# Then stop paralellization
parallel::stopCluster(cl = cl)
unregister_dopar <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
}
unregister_dopar()

#### Result visualization of stochastic mapping on 1,000 subsampled phylogeny
simmap_results <- as.data.frame(simmap_results)
colnames(simmap_results) <- simmap.colnames
dim(simmap_results)
summary(simmap_results$N) ## all 1000 iterations completed successfully

## Summary of estimations from 1,000 subsampled trees
simmap_summary <- apply(simmap_results, 2, mean)
simmap_summary 
## Transition events between geographies
simmap_summary[c(2:21)]
## values under 1 are considered insignificant and would not be counted
simmap_summary[c(2:21)][simmap_summary[c(2:21)] > 1]
## significant transitions are from Europe and from Southeast Asia.
names(simmap_summary)
## output the events that are significant
st34.trans <- simmap_results[,c(10:13,15,18:21)]
colnames(st34.trans) <- c("EU to AM", "EU to EA", "EU to OC", "EU to SEA",
                          "OC to EA", "SEA to AM", "SEA to EA", "SEA to EU", "SEA to OC")
st34.trans.long <- melt(st34.trans)
colnames(st34.trans.long) <- c("geography", 'transitions')

################################
## Draw circlize plot to visualize transitions between geography
##########################

# First create a dataframe showing directionality between regions
summary(st34.trans)
trans.cir.df <- data.frame(from=c(rep("Europe",4), "Oceania", rep("SoutheastAsia",4)),
                           to=c("America", "EastAsia", "Oceania", "SoutheastAsia", "EastAsia", "America","EastAsia", "Europe", "Oceania"),
                           value=c(3,4,2,3,1,2,6,5,8))
trans.cir.df

#chordDiagram(trans.cir.df, order = c("SEA", "Europe", "America", "EA", "Oceania"))

circos.par(gap.after=rep(5,5))
#c("#89092A", "#15AE55", "#F073DE", "#F08519", "#0E28EA")
grid.col <- c(America="#97b3d0ff", EastAsia="#002147", Europe="#f0e1b9ff", Oceania="#bf87b3", SoutheastAsia="#ed254eff" )
r1 <- chordDiagram(trans.cir.df, order = c("SoutheastAsia", "Europe", "America", "EastAsia", "Oceania"),
                   grid.col = grid.col, transparency = 0.4)
circos.clear()
## This becomes Figure 5a.

####################################
## Draw boxplots to visualize time spent in each geography
####################################
## Time spent in each geography
st34.time <- simmap_results[,c(22:26)]
st34.time.long <- melt(st34.time)
colnames(st34.time.long) <- c("geography", "time")

x1 <- ggplot(st34.time.long, aes(geography, time, fill=geography)) + geom_boxplot(color="sienna4",alpha=0.9) + theme_bw() + 
      scale_fill_manual(values = c("#97b3d0ff", "#002147", "#f0e1b9ff","#bf87b3","#ed254eff")) + 
      ylab("Proportion of time\nspent at each geography") + theme(axis.text = element_text(colour="black")) + 
      theme(legend.position = "none")
x1
## This becomes figure 5b

#########################################################3
## Randomization for robustness

# The downstream analysis aims to randomize tip-geography information to generate randomization dataset
## These will be subjected to stochastic mapping again, with 500 subtrees (iterations) to save computing power
# 10 randomization dataset will be created

## create function to randomize tip location, keeping the number of each location unchanged
randomize_location <- function(geo) {
      require(gtools)
      ran <- permute(geo)
      names(ran) <- names(geo)
      return(ran)
}

### Function to create randomization set and perform stochastic mapping 
simmap_randomize <- function(ori.dat, ori.tree, ori.state) {
      rand.geo <- randomize_location(ori.state)
      rand.dat <- ori.dat
      rand.geo <- rand.geo[match(rand.dat$ID, names(rand.geo))]
      rand.dat$Region <- rand.geo
      rand.trees <- resample_trees(dat=rand.dat, ML.tree = ori.tree, state=rand.geo, K=30, iti=500) ## 500 trees
      class(rand.trees) <- 'multiPhylo'
      rand.simm <- simmap_resample_v4(rand.trees, rand.geo) ## simmap on randomization set 
      return(rand.simm)
}

## Conduct 10 sets of randomization stochastic mapping
rdomsim1 <- simmap_randomize(metadata, tree, geo)
rdomsim2 <- simmap_randomize(metadata, tree, geo)
rdomsim3 <- simmap_randomize(metadata, tree, geo)
rdomsim4 <- simmap_randomize(metadata, tree, geo)
rdomsim5 <- simmap_randomize(metadata, tree, geo)
rdomsim6 <- simmap_randomize(metadata, tree, geo)
rdomsim7 <- simmap_randomize(metadata, tree, geo)
rdomsim8 <- simmap_randomize(metadata, tree, geo)
rdomsim9 <- simmap_randomize(metadata, tree, geo)
rdomsim10 <- simmap_randomize(metadata, tree, geo)

## Compare time spent in the geographies between the original results and those from randomization
## 
geo.time <- rbind(simmap_results[,22:26], rdomsim1[,22:26], rdomsim2[,22:26], rdomsim3[,22:26],
                  rdomsim4[,22:26], rdomsim5[,22:26], rdomsim6[,22:26], rdomsim7[,22:26],
                  rdomsim8[,22:26], rdomsim9[,22:26], rdomsim10[,22:26])
dim(geo.time)
## create another column to identify original results and randomizations
geo.time$run <- c(rep('original',nrow(simmap_results)),rep('random1',nrow(rdomsim1)), rep('random2',nrow(rdomsim2)),
                  rep('random3',nrow(rdomsim3)),rep('random4',nrow(rdomsim4)),rep('random5',nrow(rdomsim5)),
                  rep('random6',nrow(rdomsim6)),rep('random7',nrow(rdomsim7)),rep('random8',nrow(rdomsim8)),
                  rep('random9',nrow(rdomsim9)),rep('random10',nrow(rdomsim10)))

##############################################
## Proportion of time spent in each geography
p1 <- ggplot(geo.time, aes(run, America, color=run)) + geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Time spent in America') +
      scale_color_manual(values=c('red1', rep('black',10)))

p2 <- ggplot(geo.time, aes(run, `East Asia`, color=run)) + geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Time spent in East Asia') +
      scale_color_manual(values=c('red1', rep('black',10)))

p3 <- ggplot(geo.time, aes(run, Europe, color=run)) + geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Time spent in Europe') +
      scale_color_manual(values=c('red1', rep('black',10))) # significant

p4 <- ggplot(geo.time, aes(run, Oceania, color=run)) + geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Time spent in Oceania') +
      scale_color_manual(values=c('red1', rep('black',10)))

p5 <- ggplot(geo.time, aes(run, `Southeast Asia`, color=run)) + geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Time spent in Southeast Asia') +
      scale_color_manual(values=c('red1', rep('black',10))) ## significant

##########################################
### Transitions from Southeast Asia
colnames(simmap_results)
K <- c(18:21)
sea.trans <- rbind(simmap_results[,K], rdomsim1[,K], rdomsim2[,K], rdomsim3[,K],
                   rdomsim4[,K], rdomsim5[,K], rdomsim6[,K], rdomsim7[,K],
                   rdomsim8[,K], rdomsim9[,K], rdomsim10[,K])
head(sea.trans)

sea.trans$run <- c(rep('original',nrow(simmap_results)),rep('random1',nrow(rdomsim1)), rep('random2',nrow(rdomsim2)),
                   rep('random3',nrow(rdomsim3)),rep('random4',nrow(rdomsim4)),rep('random5',nrow(rdomsim5)),
                   rep('random6',nrow(rdomsim6)),rep('random7',nrow(rdomsim7)),rep('random8',nrow(rdomsim8)),
                   rep('random9',nrow(rdomsim9)),rep('random10',nrow(rdomsim10)))
colnames(sea.trans) <- c("SEA.AM", "SEA.EA", "SEA.EU", "SEA.OC", "run")

g1 <- ggplot(sea.trans, aes(run, SEA.AM, color=run))+ geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Southeast Asia to America') +
      scale_color_manual(values=c('red1', rep('black',10)))

g2 <- ggplot(sea.trans, aes(run, SEA.EA, color=run))+ geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Southeast Asia to East Asia') +
      scale_color_manual(values=c('red1', rep('black',10))) ## significant

g3 <- ggplot(sea.trans, aes(run, SEA.EU, color=run))+ geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Southeast Asia to Europe') +
      scale_color_manual(values=c('red1', rep('black',10))) ## significant

g4 <- ggplot(sea.trans, aes(run, SEA.OC, color=run))+ geom_boxplot() + theme(axis.text.x=element_blank(),  axis.title.x=element_blank(), legend.position = 'none') + 
      ylab('Southeast Asia to Oceania') +
      scale_color_manual(values=c('red1', rep('black',10))) ## significant

## Combine graphs
#### Southeast Asia 
grid.arrange(g1, g2, g3, g4, nrow=2) ## transitions from SEA
## This becomes supplementary Figure 9

#####################################
## Statistical testing
# use analysis of variance (AOV) with Tukey posthoc test to compare the estimtation from true dataset
# and those of randomization datasets
#################################

# Time spent in SEA
TukeyHSD(aov(with(geo.time, lm(`Southeast Asia` ~ run))),ordered = FALSE)$run[1:10,] ## all significant
## Time spent in other locations
TukeyHSD(aov(with(geo.time, lm(Europe ~ run))))$run[1:10, 'p adj'] <= 0.05## 7/10 are significant
TukeyHSD(aov(with(geo.time, lm(America ~ run))))$run[1:10, 'p adj'] < 0.01 # 8/10 significant
TukeyHSD(aov(with(geo.time, lm(Oceania ~ run))))$run[1:10, 'p adj'] < 0.01 # 8/10 significant

## Transitions from SEA 
TukeyHSD(aov(with(sea.trans, lm(SEA.EU ~ run))))$run[1:10, 'p adj'] < 0.05 ## 2/10 significant
TukeyHSD(aov(with(sea.trans, lm(SEA.OC ~ run))))$run[1:10, 'p adj'] < 0.05 # 9/10 significant
TukeyHSD(aov(with(sea.trans, lm(SEA.EA ~ run))))$run[1:10, 'p adj'] < 0.05 # 3/10 significant
TukeyHSD(aov(with(sea.trans, lm(SEA.AM ~ run))))$run[1:10, 'p adj'] < 0.05 # 7/10 significant

## It seems like transitions from Southeast Asia to Oceania have the most reliable
## transitions from SEA to EU is the least reliable.

##############################################
### Re-do simmap on a bigger dataset, after the addition of 91 genomes positive with blaCTX-M-55 and 16 mcr-3.1 (N=561)
##########################################
tree2 <- read.newick(file="./RAxML_bipartitions.raxml2_ST34_core_561str_postgubbins.jprooted.newick")
meta2 <- read.csv(file="./ST34_561str_short_meta.csv", header=TRUE, stringsAsFactors = FALSE)
tree2$tip.label <- str_replace_all(tree2$tip.label,"'", "")

intree2 <- match(tree2$tip.label, meta2$ID)
geo2 <- as.vector(meta2[intree2, "Region"])
table(geo2) ## minium of 43 genomes per geography
names(geo2) <- tree2$tip.label

## Create 1,000 subsampled trees with at least 40 genomes from each geography
glob.trees1000 <- resample_trees(dat=meta2, ML.tree=tree2, state=geo2, K=40, iti=1000)
class(glob.trees1000) <- 'multiPhylo'

## Redo simmap for this. This will take very long. 
### Stochastic mapping using simmap
setCores <- round(detectCores()*0.85) ## use 85% of available logical processors
cl <- parallel::makeCluster(getOption("cl.cores", setCores))
doParallel::registerDoParallel(cl = cl)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
# Main command
simmap_results_big <- simmap_resample_v4(glob.trees1000, geo2)

# Then stop paralellization
parallel::stopCluster(cl = cl)
unregister_dopar <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
}
unregister_dopar()

## summarizing for the simmap from bigger data
colnames(simmap_results_big) <- simmap.colnames
simmap_results_big <- as.data.frame(simmap_results_big)
dim(simmap_results_big)
summary(simmap_results_big$N) ## all 1,000 iterations completed successfully 

simmap_summary_big <- apply(simmap_results_big, 2, mean)

## Transition events between geographies 
simmap_summary_big[c(2:21)]
## values under 1 are considered trivial and would not be counted 
simmap_summary_big[c(2:21)][simmap_summary_big[c(2:21)] > 1]
## output events that are > 1
st34new.trans <- simmap_results_big[,c(2:4,6:8,10:13,15,18:21)]
names(st34new.trans)
st34new.trans

colnames(st34new.trans) <- c("AM to EA", "AM to EU", "AM to OC", "EA to AM", "EA to EU",
                             "EA to OC", "EU to AM", "EU to EA", "EU to OC", "EU to SEA",
                             "OC to EA", "SEA to AM", "SEA to EA", "SEA to EU", "SEA to OC")
st34new.trans.long <- melt(st34new.trans)
colnames(st34new.trans.long) <- c("geography", "transition")
## Output the visualization 
ggplot(st34new.trans.long, aes(geography, transition)) + geom_boxplot() + 
      scale_y_continuous(name="total transition events between regions", 
                         breaks = c(0,2,4,6,8,10,12,14,16))
## Evolutionary time spent in each geography
simmap_results_big[c(22:26)]
simmap_summary_big[c(22:26)] ## mean estimations
