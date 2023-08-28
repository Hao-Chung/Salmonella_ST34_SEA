## code for visualization of phylogenies and data
library(ape)
library(phytools)
library(stringr)
library(coda)
library(lubridate)
library(ggplot2)
library(ggtree)
library(phytools)
library(treeio)
library(RColorBrewer)
library(ggnewscale)
library(reshape2)
library(dplyr)
library(ggpubr)

## Import BEAST tree
beast.tre <- read.beast(file="./ST34_222str_constr.rlxlg.exp.combined.best.tre")
tipnames <- beast.tre@phylo$tip.label
tipnames <- str_replace(tipnames, "@.*", "")
beast.tre@phylo$tip.label <- tipnames
tipnames

meta2 <- read.csv(file="./ST34_454str_abricate_summary_edit2.txt", header=TRUE,
                  sep="\t", stringsAsFactors = FALSE)
which(!(tipnames %in% meta2$ID)) ## all true, names matched
colnames(meta2)

## investigate the AMR profile of all ST34
dim(meta2)
colnames(meta2)
with(meta2, table(ARR.3))
with(meta2, table(tet.AB.)) ## 421/454
with(meta2, table(sul2)) ## 416/454
with(meta2, table(strA)) ## 366/454
with(meta2, table(strB)) ## 383/454
with(meta2, table(blaTEM)) ## 365/454
with(meta2, table(IncQ1_1)) # 330/454
with(meta2, table(floR))

table(meta2$fastbaps1)
meta2$fastbaps1 <- factor(meta2$fastbaps1)
## Get a table aggregating counts of AMR genes in each BAPS population
A <- meta2 %>% 
      group_by(fastbaps1) %>%
      summarize(tet=sum(tet.AB.), sul2=sum(sul2), blaTEM =sum(blaTEM),strAB = sum(strB),
                floR=sum(floR), catA2 = sum(catA2), catB3 = sum(catB3), cmlA = sum(cmlA1), arr3 = sum(ARR.3),
                aac3iid = sum(aac.3..Iid), aac3iva = sum(aac.3..IVa), aac6Ibcr=sum(aac.6...Ib.cr),
                aph3Ia= sum(aph.3...Ia), aph4Ia=sum(aph.4..Ia))
A <- as.data.frame(A)

B <- meta2 %>% 
      group_by(fastbaps1) %>%
      summarize(blaCMY2=sum(blaCMY.2), blaCTXM14=sum(blaCTX.M.14), blaCTXM15 =sum(blaCTX.M.15), blaCTXM55=sum(blaCTX.M.55),
                blaCTXM65=sum(blaCTX.M.65),mcr1.1=sum(mcr.1.1), mcr3.11=sum(mcr.3.11), mcr3.1=sum(mcr.3.1),
                mcr3.2=sum(mcr.3.2),mcr4.1=sum(mcr.4.1), mcr5.1=sum(mcr.5.1), qnrS=sum(qnrS),mphA=sum(mph))
B <- as.data.frame(B)
bapstotal <- as.numeric(table(meta2$fastbaps1))

B[,c(-1)] <- apply(B[,c(-1)], 2, function(x) x/bapstotal*100)
Bm <- melt(B,id.vars = "fastbaps1")
colnames(Bm) <- c("BAPS", "AMRgene", "percentage")
#baps_col <- c("#ffdde2ff", "#faa094ff", "#9ed9ccff", "#ccee95")

## Down a shade of these BAPS colour using color-hex website, to increase visibility
baps_col2 <- c("#e5c6cb", "#e19085", "#8ec3b7", "#b7d686")
## Create barplot of AMR pattern between the BAPS lineages
m2 <- ggplot(Bm, aes(AMRgene, percentage, fill=BAPS)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~BAPS, ncol=4) + 
      scale_fill_manual(values = baps_col2)  +    
      scale_y_continuous(breaks=c(0,20,40)) + theme_bw() + theme(panel.grid.minor=element_blank()) 


A[,c(-1)] <- apply(A[,c(-1)],2,function(x) x/bapstotal*100)
Am <- melt(A, id.vars = "fastbaps1")
colnames(Am) <- c("BAPS", "AMRgene", "percentage")
m1 <- ggplot(Am, aes(AMRgene, percentage, fill=BAPS)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~BAPS, ncol=4) + 
      scale_fill_manual(values = baps_col2)  +    
      scale_y_continuous(breaks=c(0,25,50,75,100)) + theme_bw() + theme(panel.grid.minor=element_blank())

#####
# visualize the carriage of plasmids 
#####
plasmid.meta <- read.table(file="./ST34_plasmidfinder_summary.csv", header = TRUE,sep = ",",stringsAsFactors = FALSE)
dim(plasmid.meta)

plasmid.meta$fastbaps1 <- as.factor(plasmid.meta$fastbaps1)
table(plasmid.meta$fastbaps1)
plasmid.meta$IncFIA.total <- plasmid.meta$IncFIA.HI1._1_HI1 + plasmid.meta$IncFIA_1
plasmid.meta$IncFIB.total <- plasmid.meta$IncFIB.AP001918._1 + plasmid.meta$IncFIB.K._1_Kpn3
plasmid.meta$IncFII.total <- plasmid.meta$IncFII.p96A._1_p96A + plasmid.meta$IncFII.pCoo._1_pCoo + plasmid.meta$IncFII_1
plasmid.meta$IncI.total <- plasmid.meta$IncI1_1_Alpha + plasmid.meta$IncI2_1 + plasmid.meta$IncI2_1_Delta
plasmid.meta$IncX.total <- plasmid.meta$IncX1_1 + plasmid.meta$IncX1_4 + plasmid.meta$IncX4_1

Q <- plasmid.meta %>%
      group_by(fastbaps1) %>%
      summarize(IncFIA = sum(IncFIA.total), IncFIB=sum(IncFIB.total), IncFIC=sum(IncFIC.FII._1),
                IncFII = sum(IncFII.total), IncAC2 = sum(IncA.C2_1), IncI = sum(IncI.total), IncX=sum(IncX.total),
                IncY = sum(IncY_1), IncHI1 = sum(IncHI1A_1), IncHI2 = sum(IncHI2_1), IncN=sum(IncN_1),
                IncR = sum(IncR_1), p0111 = sum(p0111_1))
Q
Q[,c(-1)] <- apply(Q[,c(-1)], 2, function(x) x/bapstotal*100)

dim(Q)
Qm <- melt(Q)
colnames(Qm) <- c("BAPS", "plasmid", "percentage")

## Visualizing the carriage of plasmids.
m3 <- ggplot(Qm, aes(plasmid, percentage, fill=BAPS)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~BAPS, ncol=4) + 
      scale_fill_manual(values = baps_col2)  +    
      scale_y_continuous(breaks=c(0,25,50,75,100)) + theme_bw() + theme(panel.grid.minor=element_blank())
library(gridExtra)

ggarrange(m1, m2, m3, nrow=3)
## This creates the main Figure 3 in the paper

############################################
## Repeat the same procedures with the four major VN clones
#################################
table(meta2$clonal_clade)
metaVNclones <- meta2[which(meta2$clonal_clade == "VN1" | meta2$clonal_clade == "VN2" | meta2$clonal_clade == "VN3" | meta2$clonal_clade == 'VN4'),]
dim(metaVNclones)
metaVNclones$clonal_clade
table(metaVNclones$Region)
table(metaVNclones$Country) #122/130 originated from Vietnam

Avn <- metaVNclones %>% 
      group_by(clonal_clade) %>%
      summarize(tet=sum(tet.AB.), sul2=sum(sul2),blaTEM =sum(blaTEM),strAB = sum(strB),dfrA14 = sum(dfrA14),
                floR=sum(floR), catA2 = sum(catA2), cmlA = sum(cmlA1), arr3 = sum(ARR.3),
                aac3iid = sum(aac.3..Iid),aph3Ia= sum(aph.3...Ia),
                blaCTXM55=sum(blaCTX.M.55),
                mcr1.1=sum(mcr.1.1), mcr3.1=sum(mcr.3.1),
                qnrS=sum(qnrS),mphA=sum(mph))
Avn <- as.data.frame(Avn)
Avn.m <- melt(Avn)
colnames(Avn.m) <- c("clone", "AMRgene", "count")
clonecol <- brewer.pal(5, "Dark2")[2:5]
v1 <- ggplot(Avn.m, aes(AMRgene, count, fill=clone)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~clone, ncol=4) + 
      scale_fill_manual(values = clonecol)  +    
      scale_y_continuous(breaks=c(0,20,40)) + theme_bw() + theme(panel.grid.minor=element_blank()) 

table(plasmid.meta$clonal_clade)
pmetaVN <- plasmid.meta[which(plasmid.meta$clonal_clade == "VN1" | plasmid.meta$clonal_clade == "VN2" | plasmid.meta$clonal_clade == "VN3" | plasmid.meta$clonal_clade == 'VN4'),]
dim(pmetaVN)

Qvn <- pmetaVN %>%
      group_by(clonal_clade) %>%
      summarize(IncFIB=sum(IncFIB.total), IncFIC=sum(IncFIC.FII._1),
                IncFII = sum(IncFII.total), IncAC2 = sum(IncA.C2_1), IncI = sum(IncI.total), IncX=sum(IncX.total),
                IncY = sum(IncY_1), IncHI2 = sum(IncHI2_1), IncN=sum(IncN_1), p0111 = sum(p0111_1))
Qvn
Qvn.m <- melt(Qvn)
colnames(Qvn.m) <- c("clone", "plasmid", "count")
v2 <- ggplot(Qvn.m, aes(plasmid, count, fill=clone)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~clone, ncol=4) + 
      scale_fill_manual(values = clonecol)  +    
      scale_y_continuous(breaks=c(0,20,40)) + theme_bw() + theme(panel.grid.minor=element_blank()) 
ggarrange(v1, v2, ncol=1)
## This becomes Supplementary Figure 2

##########
# Visualizing BEAST tree
###
meta2.ed <- meta2[,c("ID", "Region", "serovar", "number_acquired_AMR_genes","clonal_clade",
                     "blaCTX.M.55", "qnrS", "mcr.3.1", "mph", "IncA.C2_1", "IncHI2_1")]
meta2.ed$IncA.C2_1
meta2.ed$IncHI2_1

colnames(meta2.ed) <- c("ID", "Region", "serovar", "NoAMRgenes", "clade", 
                        "blaCTX-M-55", "qnrS", "mcr3-1", "mphA", "IncA/C", "IncHI2")

rownames(meta2.ed) <- meta2.ed$ID
meta2.ed <- meta2.ed[match(tipnames, meta2.ed$ID),]
meta2.ed$Region <- as.factor(meta2.ed$Region)

dfheat <- meta2.ed[,c(3,6:11)]
dfheat$`IncA/C` <- str_replace(dfheat$`IncA/C`, "\\.", "0")
dfheat$`IncA/C` <- as.numeric(dfheat$`IncA/C`)
dfheat$IncHI2 <- str_replace(dfheat$IncHI2, "\\.", "0")
dfheat$IncHI2 <- as.numeric(dfheat$IncHI2)

head(dfheat)

## setting margin
par(mar=c(0.5, 0.5, 1.5, 0.5))

# set posterior at tip of tree as 1 for accurate visualization
naposterior <- which(is.na(beast.tre@data$posterior))
beast.tre@data$posterior[naposterior] <- 1

p <- ggtree(beast.tre, mrsd="2019-11-01", size=0.6, aes(color=posterior)) + 
      theme_tree2() +
      scale_color_continuous(low='skyblue1', high='black') + theme(legend.position="left") +
      scale_x_continuous(breaks=c(1996,2000,2010,2020), minor_breaks=seq(1996,2020,2)) +
      theme(panel.grid.major = element_line(color='grey50', size=.2),
            panel.grid.minor = element_line(color='grey', size=.2, linetype=3),
            panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())

p1 <- p %<+% meta2.ed + new_scale_color() + geom_tippoint(aes(color=Region), size=2.5, alpha=1, shape=20) + 
      scale_colour_manual(name=levels(factor(meta2.ed$Region)), values = c("#97b3d0ff", "#002147", "#f0e1b9ff", "#bf87b3", "#ed254eff"))

p1

## find nodes that represent VN clades
tre1 <- beast.tre@phylo
library(caper)
trenodelist <- clade.members.list(tre1, tip.labels = TRUE, include.nodes = FALSE)
VN1 <- meta2.ed[which(meta2.ed$clade == "VN1"), "ID"]
VN2 <- meta2.ed[which(meta2.ed$clade == "VN2"), "ID"]
VN3 <- meta2.ed[which(meta2.ed$clade == "VN3"), "ID"]
VN4 <- meta2.ed[which(meta2.ed$clade == "VN4"), "ID"]
Vbsi <- meta2.ed[which(meta2.ed$clade == "V-BSI"), "ID"]

findnode <- function(tree, tips) {
      require(caper)
      trenodelist <- clade.members.list(tree, tip.labels = TRUE, include.nodes = FALSE)
      kk <- c()
      for (i in 1:length(trenodelist)) {
            kk[i] <- all(all(trenodelist[[i]] %in% tips) & (length(trenodelist[[i]]) == length(tips)))
      }
      val <- names(trenodelist)[which(kk)]
      return(val)
}

ann_node <- c(findnode(tre1, VN1),findnode(tre1, VN2), findnode(tre1, VN3),
              findnode(tre1, VN4), findnode(tre1, Vbsi))
dt <- data.frame(node=as.numeric(ann_node), name=c("VN1", "VN2", "VN3", "VN4", "V-BSI"))
###
p1 <- p1 + new_scale_color() + geom_hilight(data=dt, mapping=aes(node=node, fill=name), alpha=0.5, type="rect")

p2 <- gheatmap(p1, dfheat, offset=0.07, width=0.3, colnames_position = "top", font.size=3, color="grey80",
               colnames_angle=30, colnames_offset_y=0.09, hjust=0.1) + theme(legend.title = element_blank(), legend.position = "left") + 
      scale_fill_manual(values=c("honeydew", "grey0", "#eea47fff", "#00539cff", brewer.pal(5, "Dark2")))
p2
## This becomes figure 2a

### Visualize original big tree of 454 strains STm
bigtre <- read.newick(file="./RAxML_bipartitions.raxml_ST34rmblk_454str")
bigtre$tip.label
ggtree(bigtre)

### reroot the tree based on most probable root 
findnode(bigtre, c("SAMD00234571", "SAMD00097511", "SAMD00097512")) # 689
roottre <- reroot(bigtre, node.number = 689,position=0.0011)
ggtree(roottre)

all(roottre$tip.label %in% meta2$ID) ## all ids are correct
which(roottre$node.label == "") # 36
roottre$node.label <- c(100,100,roottre$node.label[2:35], roottre$node.label[37:453])
roottre$node.label <- as.numeric(roottre$node.label)

g <- ggtree(roottre, size=0.4,layout='fan', open.angle=10,
            mapping=aes(color=c(rep(100,454), as.numeric(label)[455:907]))) + 
      scale_color_gradient(low='cyan', high='black') + theme(legend.position="left") + labs(color="bootstrap")

### 
meta3 <- meta2
rownames(meta3) <- meta3$ID
colnames(meta3)
head(meta3)

meta3.ed <- meta3[,c(6,3,9,5,7)]
colnames(meta3.ed) <- c("serovar", "Region", "BAPS", "source", "disease")

#test <- meta3[,c(9,8)]
table(meta3.ed$serovar)
table(meta3.ed$Region)
table(meta3.ed$source)
meta3.ed$source <- str_replace(meta3.ed$source,"0"," ")
table(meta3.ed$disease)
table(meta3.ed$BAPS)

###
## annotate Vietnamese clades
table(meta3$clonal_clade)
table(meta3$BAPS)
Vbsi_b <- meta3[which(meta3$clonal_clade == "V-BSI"),"ID"]
VN1_b <- meta3[which(meta3$clonal_clade == "VN1"),"ID"]
VN2_b <- meta3[which(meta3$clonal_clade == "VN2"),"ID"]
VN3_b <- meta3[which(meta3$clonal_clade == "VN3"),"ID"]
VN4_b <- meta3[which(meta3$clonal_clade == "VN4"),"ID"]
mc5 <- meta3[which(meta3$clonal_clade == "mclade5"),"ID"]
mc7 <- meta3[which(meta3$clonal_clade == "mclade7"),"ID"]
mc8 <- meta3[which(meta3$clonal_clade == "mclade8"),"ID"]
mc9 <- meta3[which(meta3$clonal_clade == "mclade9"),"ID"]

ann_node2 <- c(findnode(roottre, Vbsi_b), findnode(roottre, VN1_b), findnode(roottre, VN2_b),
               findnode(roottre, VN3_b), findnode(roottre, VN4_b), findnode(roottre, mc5), 
               findnode(roottre, mc7), findnode(roottre, mc8), 
               findnode(roottre, mc9))

dt2 <- data.frame(node=as.numeric(ann_node2), name=c("V-BSI","VN1", "VN2", "VN3", "VN4", "mc5", "mc7", 
                                                     "mc8", "mc9"))

g <- g + new_scale_color() + geom_hilight(data=dt2, mapping=aes(node=node, fill=name), alpha=0.5, type="rect")

brks <- c(" ", levels(factor(meta3.ed$serovar)), levels(factor(meta3.ed$Region)),
          levels(factor(meta3.ed$source))[2:6], levels(factor(meta3.ed$disease))[2:4],
          levels(factor(meta3.ed$BAPS)),dt2$name)

col_val2 <- c("white", "#eea47fff", "#00539cff","#97b3d0ff", "#002147", "#f0e1b9ff", "#bf87b3", "#ed254eff",
              "#d7a9e3ff", "#dfdce5ff", "#97b3d0ff", "#e7ebe0ff",  "#dbb04aff",
              "#f0e1b9ff", "#a13941ff","#b3c7d6ff", 
              "#ffdde2ff", "#faa094ff", "#9ed9ccff", "#ccee95",
              brewer.pal(5, "Dark2"), "grey20", "grey40", "grey60", "grey80")

g1 <- gheatmap(g, meta3.ed, offset=0, width=0.4, colnames_position = "top", font.size=2, color=NULL,
               colnames_angle=80, colnames_offset_y=0.2, hjust=0) + theme(legend.title = element_blank(), legend.position = "left") + 
      scale_fill_manual(breaks=brks, values=col_val2)
g1 
## This becomes main Figure 1 in the paper
### examining the over-presence of AMR genes in SEA isolates.
meta3$SEA <- meta3$Region
meta3[which(meta3$SEA != 'Southeast Asia'), 'SEA'] <- 'other'
chisq.test(with(meta3, table(SEA, qnrS))) ## p value < 2.2e-16
chisq.test(with(meta3, table(SEA, blaCTX.M.55))) # p value < 1.27e-06
chisq.test(with(meta3, table(SEA, mph))) # p value = 0.0002866
chisq.test(with(meta3, table(SEA, mcr.3.1))) # p value = 9.359e-7

## Number of AMR genes by time
## Excluding the constiutive gene acc(6)-Iaa
table(meta3$aac.6...Iaa) ## almost universal present
meta3$number_acquired_AMR_genes <- meta3$number_acquired_AMR_genes - 1
VNdat <- meta3[which(meta3$Country == "Vietnam"),]
dim(VNdat)
ggplot(VNdat, aes(year, number_acquired_AMR_genes)) + geom_point() ## no obvious trend in time
ggplot(VNdat, aes(clonal_clade, number_acquired_AMR_genes)) + geom_boxplot()

## number of AMR genes by clone
## only include the major clones
VNdat2 <- VNdat[str_detect(VNdat$clonal_clade, "V",),]
VNdat2$number_acquired_AMR_genes
k1 <- ggplot(VNdat2, aes(clonal_clade, number_acquired_AMR_genes, fill=clonal_clade,alpha=0.7)) + geom_boxplot(color="grey20", width=0.5) +
      theme_bw() + coord_flip() + 
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none") + 
      ylab("Number of acquired AMR genes") + 
      scale_fill_manual(values=brewer.pal(5,"Dark2"))

TukeyHSD(aov(with(VNdat2, lm(number_acquired_AMR_genes ~ clonal_clade))))
k1
## This becomes figure 2c

#########
## Read in beast trees file to summarize the divergence dates of major clones
###########
library(rBt)
library(beastio)
st34.trees <- readTreeLog(file="./ST34_222str_constr.rlxlg.exp.combined.trees", burnin = 0)
beastclone <- read.csv(file="./beast_tree_annotation.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
vn1 <- beastclone[beastclone$clone == "VN1",'beast_name']
vn2 <- beastclone[beastclone$clone == "VN2",'beast_name']
vn3 <- beastclone[beastclone$clone == "VN3",'beast_name']
vn4 <- beastclone[beastclone$clone == "VN4",'beast_name']
vbsi <- beastclone[beastclone$clone == "VBSI",'beast_name']

## the max year of sampling is 2019.833 
himax <- 2019.833

length(st34.trees)
vn1.tmrca <- rep(0,length(st34.trees))
vn2.tmrca <- rep(0, length(st34.trees))
vn3.tmrca <- rep(0, length(st34.trees))
vn4.tmrca <- rep(0, length(st34.trees))
vbsi.tmrca <- rep(0, length(st34.trees))

## output all the data points for each clone.
for (x in 1:length(st34.trees)) {
      vn1.tmrca[x] <- getCladeHeight(st34.trees[[x]], vn1)
}
for (x in 1:length(st34.trees)) {
      vn2.tmrca[x] <- getCladeHeight(st34.trees[[x]], vn2)
}
for (x in 1:length(st34.trees)) {
      vn3.tmrca[x] <- getCladeHeight(st34.trees[[x]], vn3)
}
for (x in 1:length(st34.trees)) {
      vn4.tmrca[x] <- getCladeHeight(st34.trees[[x]], vn4)
}
for (x in 1:length(st34.trees)) {
      vbsi.tmrca[x] <- getCladeHeight(st34.trees[[x]], vbsi)
}
mrcadat <- cbind(vbsi.tmrca, vn1.tmrca, vn2.tmrca, vn3.tmrca, vn4.tmrca)
colnames(mrcadat) <- c("VBSI", "VN1", "VN2", "VN3", "VN4")
mrcadat <- himax - mrcadat 
mrcamelt <- melt(mrcadat)
colnames(mrcamelt) <- c("var1", "clone", "year")

ggplot(mrcamelt, aes(clone, year, color=clone)) + geom_boxplot(width=0.6,size=0.8,alpha=0.6) +
      scale_y_continuous(breaks=c(1996,2000,2010,2020),limits = c(1996,2020),minor_breaks = seq(from=1996,to=2020,by=2)) +
      coord_flip() + theme_bw() + theme(axis.text.y = element_blank(), axis.title = element_blank(),legend.position ='none', axis.ticks.y=element_blank()) + 
      scale_color_manual(values = brewer.pal(5, "Dark2"))
## This becomes figure 2b 

############################
### After including more global isolates, mostly of ones carrying blaCTX-M-55 and mcr-3.1
############################
glob.tre <- read.newick(file="./RAxML_bipartitions.raxml_ST34_panaroocore_329str")
glob.tre$tip.label
ggtree(glob.tre)
findnode(glob.tre, c("SAMD00234571", "SAMD00097511", "SAMD00097512")) # node 410
rootglob.tre <- reroot(glob.tre, node.number = 410)

ggtree(rootglob.tre)
which(rootglob.tre$node.label == "") #11 
rootglob.tre$node.label <- c(100,100,rootglob.tre$node.label[2:10],rootglob.tre$node.label[12:328])
## Read in metadata used for global tre
metaglob <- read.csv(file="./Supplementary_Data_1.csv", header=TRUE, sep=",", 
                     stringsAsFactors = FALSE)
is.na(match(rootglob.tre$tip.label, metaglob$Sample_ID)) ## all taxa are included
rootglob.tre$node.label <- as.numeric(rootglob.tre$node.label)

w <- ggtree(rootglob.tre, size=0.6,
            mapping=aes(color=c(rep(100,329), as.numeric(label)[330:657]))) + 
      scale_color_gradient(low='cyan', high='black') + theme(legend.position="left") + labs(color="bootstrap")
w

colnames(metaglob)
metaglob$newdata <- NA
metaglob[456:nrow(metaglob),"newdata"] <- 1
colnames(metaglob)
meta5 <- metaglob[,c("Sample_ID", "serovar","Region", "clonal_clade","IncA.C2_plasmid",
                     "blaCTX.M.55", "newdata")]
rownames(meta5) <- meta5$Sample_ID
meta5$newdata <- as.factor(meta5$newdata)

w1 <- w %<+% meta5 + new_scale_color() + geom_tippoint(aes(color=Region), size=2.5, alpha=1, shape=20) + 
      scale_colour_manual(name=levels(factor(meta2.ed$Region)), values = c("#97b3d0ff", "#002147", "#f0e1b9ff", "#bf87b3", "#ed254eff"))
w1

meta6 <- meta5[which(meta5$Sample_ID %in% rootglob.tre$tip.label),]

## annotate Vietnamese clades
table(meta6$clonal_clade)

VN2_g <- meta6[which(meta6$clonal_clade == "VN2"),"Sample_ID"]
VN4_g <- meta6[which(meta6$clonal_clade == "VN4"),"Sample_ID"]
Aus1 <- meta6[which(meta6$clonal_clade == "Australia_1"),"Sample_ID"]
Cnchrom <- meta6[which(meta6$clonal_clade == "chrombla_China"),"Sample_ID"]
Ausm <- meta6[which(meta6$clonal_clade =="SEA_min"), "Sample_ID"]

ann_node3 <- c(findnode(rootglob.tre, VN2_g), findnode(rootglob.tre, VN4_g), findnode(rootglob.tre, Aus1),
               findnode(rootglob.tre, Cnchrom), findnode(rootglob.tre, Ausm))
#w %<+% meta5 + new_scale_color() + geom_tippoint(aes(color=newdata), size=2.5, alpha=1, shape=20, na.rm=TRUE)
dt3 <- data.frame(node=as.numeric(ann_node3), name=c("VN2", "VN4", "Australia1", "China", "SEAm"))
w1 <- w1 + new_scale_color() + geom_hilight(data=dt3, mapping=aes(node=node, fill=name), alpha=0.3, type="rect")
w1

meta7 <- meta6[,c(7,2,5,6)]
meta7$newdata <- as.numeric(meta7$newdata)
meta7[is.na(meta7$newdata) , "newdata"] <- 0
meta7$IncA.C2_plasmid <- as.numeric(str_replace(meta7$IncA.C2_plasmid, "\\.", "0"))
meta7$blaCTX.M.55 <- as.numeric(str_replace(meta7$blaCTX.M.55, "\\.", "0"))

brks3 <- c(0,1,levels(factor(meta7$serovar)),dt3$name)
col_val3 <- c("honeydew", "grey0", "#eea47fff", "#00539cff",
              brewer.pal(5,"Dark2")[c(3,5)], "cadetblue4","brown1","darkgoldenrod4")
w2 <- gheatmap(w1, meta7, offset=0, width=0.2, colnames_position = "top", font.size=3, color=NULL,
               colnames_angle=30, colnames_offset_y=0.2, hjust=0.1) + theme(legend.title = element_blank(), legend.position = "left") + 
      scale_fill_manual(breaks=brks3, values=col_val3)
w2
## This then becomes Supplementary Figure 8