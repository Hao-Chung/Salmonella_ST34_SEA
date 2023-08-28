library(stringr)
library(foreach)
library(rhierbaps)
library(geiger)
library(seqinr)
library(genoPlotR)
library(phytools)
library(ggplot2)

ann_path <- "../annotation_files/"
filelist <- list.files(ann_path, pattern ="*.txt",full.names = TRUE)

ann_txt <- lapply(filelist, function(x) read.table(x, header=T, stringsAsFactors = FALSE))
ann_txt[[1]]

filenames <- str_replace(str_replace(filelist,pattern = "../annotation_files//",""),
            pattern="_bwafrbTWStm6.frb.fil.snpSift_extract.txt","")
names(ann_txt) <- filenames

KK <- foreach (i=1:length(filenames), .combine="rbind") %do% {
      ann_txt[[i]][,c("POS", "REF","ALT","ANN.0..EFFECT","ANN.0..GENE")]
}
sort(unique(KK$POS),decreasing = FALSE)
## Exploring 
KK[KK$POS==852,]

snp_fasta <- load_fasta(msa="../ST34_only_454str_bwafrbTWStm6_cons_rmblk.snp.fa", keep.singletons = TRUE)
snp_fasta <- t(snp_fasta)
colnames(snp_fasta) <- str_replace(colnames(snp_fasta),"_trimmed", "")
dim(snp_fasta)
## read in the original position of the snp alignment
pos <- read.csv(file="../snp_alignment_TWStm6_position.csv", header=FALSE,stringsAsFactors = FALSE)
pos <- pos$V2

which(!(pos %in% unique(KK$POS))) ## all pos value are within the vcf annotation positions.
A <- KK[match(pos, KK$POS),]
rownames(A) <- NULL
dim(A)
colnames(A) <- c("chrom_position","ref", "snp", "ann_effect", "gene")
table(A$ann_effect)
A[which(A$ann_effect == "downstream_gene_variant" | A$ann_effect == "upstream_gene_variant"),"coding"] <- "intergenic"
A[which(is.na(A$coding)),'coding'] <- "CDS"

## read in tRNA list file
tRNA_list <- read.csv(file="../TWStm6_tRNA-list.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
A[match(tRNA_list$gene, A$gene,nomatch = 0),'coding'][c(1:3,5:10,12,13)] <- 'tRNA'
A[A$coding == 'intergenic',]
A[A$ann_effect == "synonymous_variant","syno_nonsyn"] <- "S"
A[A$ann_effect == "missense_variant", "syno_nonsyn"] <- "N"
A[A$ann_effect == "stop_gained", "syno_nonsyn"] <- "STOP"
A[A$ann_effect == "start_lost" | A$ann_effect == "splice_region_variant&stop_retained_variant" |
      A$ann_effect == "stop_lost&splice_region_variant", "syno_nonsyn"] <- 'other'
A[(which(is.na(A$syno_nonsyn))),"syno_nonsyn"] <- "-" 
head(A)

## Read in annotation for genes
gene_list <- read.csv(file="../TWStm6_all_CDS_list.csv", header=TRUE, stringsAsFactors = FALSE)
gene_list$strand <- str_replace(gene_list$strand,"c","-1")
gene_list[which(gene_list$strand != "-1"),"strand"] <- "1"
gene_list$strand <- as.numeric(gene_list$strand)

A$strand <- gene_list[match(A$gene, gene_list$gene),"strand"]
A$product <- gene_list[match(A$gene, gene_list$gene),'product']
dim(A)
dim(snp_fasta)

snp_fasta <- toupper(snp_fasta)
## Find ambiguous character in the alignment
find_gap <- function(x) {
      return(length(which(x == "-")))
}
no_ambiguous <- apply(snp_fasta, 1, find_gap)
summary(no_ambiguous)
## Find the number of SNPs at a position
#find_snp <- function(x) {
 #     table_snp <- sort(table(x),decreasing=TRUE)
  #    y <- sum(table_snp[which(names(table_snp) != "-")][-c(1)])
   #   return(y)
#}

ss <- c()
for(i in 1:nrow(snp_fasta)) {
      table_snp <- table(snp_fasta[i,])
      ref_base <- A[i, "ref"]
      ss[i] <- paste(as.vector(table_snp[which(names(table_snp) != "-" & names(table_snp) != ref_base)]), collapse=",")
}

## checking how many allele per position
calc_allele <- function(x) {
      table_snp <- table(x)
      y <- length(which(names(table_snp)!= "-"))
      return(y)
}
no_allele <- apply(snp_fasta, 1, calc_allele)      
summary(no_allele)
which(no_allele==3)
table(snp_fasta[which(no_allele == 3)[2],])

#no_snp <- apply(snp_fasta, 1, find_snp)
#no_snp
cbind(no_ambiguous, ss)
head(A)
B <- cbind(A[,c(1,6,8,5,9,7,2,3)], no_ambiguous, ss, snp_fasta)

write.table(B, file="./ST34_454str_rmblk_snp_metadata.txt", quote = FALSE,sep = "\t", col.names = TRUE, row.names = FALSE)

##############
# Read in edited snp metadata file
##########

variants <- read.csv(file="../ST34_454str_rmblk_snp_metadata_edit.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")

### Read in fasta file of the reference
readFastaRef = function(refFile) {
      row = scan(refFile,what=character(0),sep="\n")
      # Convert each base to an individual character
      chars = substr(row,1,1)
      base = chars!=">"
      seq = paste(row[base],collapse="")
      return(toupper(unlist(strsplit(seq,""))))
}
TWStm6 <- readFastaRef("../TWStm6.fa")
table(TWStm6)
chromLen <- length(TWStm6)

colnames(variants)
colnames(variants) <- str_replace(colnames(variants), "^X", "")
## columns from 12 to 465 refer to 454 samples
varHeader <- colnames(variants)

nSeq <- length(varHeader[12:465])
newSeqLabels <- varHeader[12:465]

## Extract base calls and transpose
vSites <- t(variants[,c(12:465)])
vPos <- variants$chrom_position
table(variants$ref, TWStm6[vPos]) ## all matched

##########
## processing the output from PAML 

ml.tree <- read.newick(file="./RAxML_bipartitions.raxml_ST34rmblk_454str.japanroot.newick")
node.layout <- readNexus(file="./ST34_454str_ancrecon_tree_internalnodelayout.nexus",format = "raxml")

ml.tree$tip.label
ml.tree$tip.label <- str_replace_all(ml.tree$tip.label, pattern = "'", "")

node.layout$tip.label
node.layout$tip.label <- str_replace_all(node.layout$tip.label, pattern="'", "")
node.layout$tip.label <- str_replace_all(node.layout$tip.label,pattern = "^[0-9]*_","")

ml.tree$tip.label == node.layout$tip.label
identical(ml.tree$tip.label, node.layout$tip.label) ## all true

node.layout$node.label
ml.tree$node.label <- paste0("N",node.layout$node.label)

##############
## classification of substitutions

gbk <- read_dna_seg_from_file("../TWStm6.gbk", fileType="genbank", extra_fields = "*matching_locus_tag")
## read in the imputed sequence of ancestral reconstruction by PAML 
MLSeqFile <- scan("./ST34_454str_paml_ancrecon_joint_reconstruction.fasta", what=character(0))
MLSeqs <- MLSeqFile[seq(2,length(MLSeqFile), by =2)]

#In order to identify the reconstructed sequences corresponding to each node in the tree, we can store the sequences in a matrix 
# Create an empty matrix of reconstructed nucleotide of 4962 SNP position, for 454 taxa plus 453 internal nodes
M <- matrix("", length(MLSeqs), nchar(MLSeqs[1]))
# Fill each row of M, where each column is a variant site and each row is a sequence
for (i in 1:length(MLSeqs)) {
      seq <- as.matrix(unlist(strsplit(MLSeqs[i], "")))
      M[i,] <- seq
}
## Attach names to the M matrix
seqLabels <- gsub(">", "", MLSeqFile[seq(1,length(MLSeqFile),by=2)])
rownames(M) <- seqLabels

## Get all internal node and tip labels
nodeLabels <- c(ml.tree$tip.label, ml.tree$node.label)
# Reorder M so that sequences are in the same order as nodeLabels
M <- M[match(nodeLabels, rownames(M)),]
# For each node, get the row index of its ancestor
M_ancestorNodeIndex <- sapply(1:nrow(M), function(x) {
      edgeRowNumber <- match(x, ml.tree$edge[,2])
      return(ml.tree$edge[edgeRowNumber, 1])
})
M_ancestorNodeIndex

## Provide the genetic code for translation
geneticCode = list(
      "TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
      "TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
      "TAT"="Tyr","TAC"="Tyr","TAA"="STO","TAG"="STO",
      "TGT"="Cys","TGC"="Cys","TGA"="STO","TGG"="Trp",
      "CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
      "CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
      "CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
      "CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
      "ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
      "ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
      "AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
      "AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
      "GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
      "GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
      "GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
      "GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")
complement <- c("A"="T","C"="G","G"="C","T"="A","N"="N")

## Calculate the number of alleles at each variable site
nNuc <- apply(vSites, 2, function(x) table(factor(x, levels=c("A", "C", "G", "T"))))
nAlleles <- colSums(nNuc>0)

#Some specialized functions are needed. The first function reads in the index (i) of each variant site in the vector vPos and the number of alleles at site i.
#It returns a vector of information about substitutions at that site.
subType = function(i, numberOfAlleles){
      # Set all counts of substitution types to zero
      nonSynCount = 0
      synCount = 0
      readThroughCount = 0
      nonsenseCount = 0
      position = as.numeric(vPos[i])
      # We will fill in the entries in codonInfo below
      codonInfo = c("Type"="", "NonSynCount"=0,"SynCount"=0,
                    "ReadThroughCount"=0, "NonsenseCount"=0, "Refcodon"="",
                    "Nonrefcodon"="-","Frame"="-","Codonposition"="-",
                    "Refaa"="-","Nonrefaa"="-")
      refAllele = TWStm6[position]
      if(is.na(position)) {
            codonInfo[1] = "Unknown"
            return(codonInfo)
      }
      geneIndex = which((gbk$feature=="CDS" & gbk$start<=position & gbk$end>=position))
      if(length(geneIndex)==0) {
            codonInfo[1] = "Intergenic"
            return(codonInfo)
      }
      if(numberOfAlleles==1){
            codonInfo[1] = "Invariant"
            return(codonInfo)
      }
      # Find reading frame of gene and indices of all 3 sites in codon
      if(gbk$strand[geneIndex]==1) { # Codon is on the positive strand
            frame = (position-gbk$start[geneIndex]) %% 3
            codonPositions = position-frame+(0:2)
            refCodon = paste(TWStm6[codonPositions],collapse="")
      }else{ # Codon is on the negative strand and requires reverse complementation
            frame = (gbk$end[geneIndex]-position) %% 3
            codonPositions = position+frame-(0:2)
            refCodon = paste(complement[TWStm6[codonPositions]],collapse="")
      }
      for(j in 1:nrow(M)){ # For each branch in tree
            # Get node and ancestral bases
            nodeVarBase = M[j, i]
            ancestralNodeNumber = M_ancestorNodeIndex[j]
            ancestorVarBase = M[ancestralNodeNumber, i]
            # If not root and bases are different
            if(!is.na(ancestorVarBase) && nodeVarBase != ancestorVarBase){
                  nodeCodon = TWStm6[codonPositions] # Set node codon to reference
                  ancestorCodon = TWStm6[codonPositions] # Set ancestral codon to reference
                  for(x in c(1, 2, 3)){
                        if(codonPositions[x] %in% vPos){
                              # Get column index in matrix M 
                              M_columnIndex = match(codonPositions[x], vPos)
                              # Substitute codon position in node and ancestral codons for correct bases
                              nodeBase = M[j, M_columnIndex]
                              nodeCodon[x] = nodeBase 
                              ancestralBase = M[ancestralNodeNumber, M_columnIndex]
                              ancestorCodon[x] = ancestralBase
                        }
                  }
                  if(gbk$strand[geneIndex]==-1) { # If negative strand, complement bases
                        nodeCodon = as.vector(complement[nodeCodon])
                        ancestorCodon = as.vector(complement[ancestorCodon])
                  }
                  # If "N" in either codon, we can't count any substitutions
                  if(!"N" %in% nodeCodon && !"N" %in% ancestorCodon){      
                        nodeCodonCol = paste0(nodeCodon, collapse="")
                        ancestorCodonCol = paste0(ancestorCodon, collapse="")
                        nodeAA = as.character(geneticCode[nodeCodonCol])
                        ancestralAA = as.character(geneticCode[ancestorCodonCol])
                        if(nodeAA == ancestralAA){
                              synCount = synCount + 1
                        }else if(ancestralAA == "STO"){
                              readThroughCount = readThroughCount + 1
                        }else if(nodeAA == "STO") {
                              nonsenseCount = nonsenseCount + 1
                        }else{
                              nonSynCount = nonSynCount + 1
                        }
                  }
            }
      }
      codonInfo = c("Type"="Coding", "NonSynCount"=nonSynCount,
                    "SynCount"=synCount,"ReadThroughCount"=readThroughCount,
                    "NonsenseCount"=nonsenseCount, "Refcodon"=refCodon,
                    "Nonrefcodon"="-","Frame"=frame+1,
                    "Codonposition"=((position-gbk$start[geneIndex]) %/% 3)+1,
                    "Refaa"="-","Nonrefaa"="-")
      return(codonInfo)
}

## a function is also needed to return the gene name(s) within which a site lies
gbkLocate = function(gbk,position) {
      position = as.numeric(position)
      geneIndex = which(gbk$start<=position & position<=gbk$end)
      start = gbk$start[geneIndex]
      end = gbk$end[geneIndex]
      noMatch = length(start)==0 | length(end)==0
      if(!noMatch){
            noMatch = noMatch | is.na(start) | is.na(end)
      }
      if(noMatch) {
            tp = gbk[1,]
            tp[1,1:ncol(tp)] = "-"
            return(tp)
      }
      as.data.frame(t(apply(gbk[geneIndex,],2,paste0,collapse=":")))
}
### Finally, using the functions above, construct a table called 'subAnnotations', 
#which contains counts of non-synonymous and synonymous substitutions for each variant site, 
#in addition to information about the site.

## But first, we need to remove a few sites
# 
ori_vPos <- vPos

vPos <- ori_vPos[-c(1682, 2138,4627,4863)] ##  delete positions in vPos which SNPs fall into DNA regions overlapping two CDS
nNuc <- nNuc[,-c(1682, 2138,4627,4863)]
nAlleles <- nAlleles[-c(1682, 2138,4627,4863)]
ori_M <- M
M <- ori_M[,-c(1682, 2138,4627,4863)]

subAnnotation = list()
for(i in 1:length(vPos)){
      alleles = c("A","C","G","T")[order(nNuc[,i])][c(4,3)]
      siteInfo = data.frame("Position"=vPos[i],"MajorAllele"=alleles[1],
                            "MinorAllele"=alleles[2], "NumberOfAlleles"=nAlleles[i],
                            stringsAsFactors = FALSE)
      subInfo = t(subType(i, siteInfo$NumberOfAlleles))
      geneInfo = gbkLocate(gbk,siteInfo$Position[1])[,c("name","start","end","strand",
                                                        "length","pid","gene","synonym",
                                                        "product","proteinid","feature",
                                                        "*matching_locus_tag")]
      # Take the first gene if site in multiple overlapping genes
      subAnnotation = rbind(subAnnotation,
                            cbind(siteInfo[1,], subInfo, geneInfo[1,]))
}
tail(subAnnotation)

subAnnotation = transform(subAnnotation,
                          NonSynCount = as.numeric(as.character(NonSynCount)),
                          SynCount = as.numeric(as.character(SynCount)),
                          ReadThroughCount = as.numeric(as.character(ReadThroughCount)),
                          NonsenseCount = as.numeric(as.character(NonsenseCount)))
table(subAnnotation$Type)
head(subAnnotation)

variants[variants$gene == "ftsK",c(1:10)] ## large number of missing at position 1017873
# mutations in ftsK at this position is not reliable 
# So delete the annotations at this position
subAnnotation <- subAnnotation[which(!(subAnnotation$Position == 1017873)),]

### Now that we have counted and classified all subsitutions in the tree, 
#we need to calculate empirical codon usage from the reference genome before we can quantify the selection acting on the data

## Counting the frequency of codon appearance
concatTranscripts = c()
for(i in 1:nrow(gbk)) {
      if(gbk$strand[i]==1){
            concatTranscripts = c(concatTranscripts,apply(matrix(
                  toupper(TWStm6[gbk$start[i]:gbk$end[i]]),nrow=3),2,paste,collapse=""))
      }
      if(gbk$strand[i]==-1){
            concatTranscripts = c(concatTranscripts,apply(matrix(complement[
                  toupper(TWStm6[gbk$end[i]:gbk$start[i]])],nrow=3),2,paste,collapse=""))
      }
}
codonUsage = table(factor(concatTranscripts,levels=names(geneticCode)))
#Using the list of information about each substitution created in the subAnnotation variable, 
#count the number of synonymous and non-synonymous substitutions that occur in each gene.

geneNames = gbk$name
# Rename duplicate gene names
dup_genes <- geneNames[duplicated(geneNames)]
unique(variants[which(variants$gene %in% dup_genes),"gene"]) # 69 genes among the duplicated names
#for (i in 1:length(dup_genes)) {
#     geneNames[which(geneNames %in% dup_genes[i])[2]] <- str_c(dup_genes[i], "_2")
#}
#geneNames[duplicated(geneNames)]

## extract using the synonym field rather than names
which(table(gbk$synonym) > 1) # only one duplicate
geneNames <- gbk$synonym
geneNames[duplicated(geneNames)]
geneNames <- geneNames[!duplicated(geneNames)]

## Lengths of genes
geneLengths = xtabs((gbk$end-gbk$start+1)~factor(gbk$synonym,levels=geneNames))

rawNS <- data.frame(matrix(0, nrow=length(geneNames), ncol=2))
rownames(rawNS) <- geneNames
colnames(rawNS) <- c("N", "S")
## Finding all coding substitutions
codingSubs <- subAnnotation[which(subAnnotation$Type =='Coding'),]
# convert counts of non-synonymous and synonymous substitutions to numeric objects
aggCounts <- aggregate(codingSubs[,6:9], list(as.character(codingSubs$synonym)), sum)
## Count synonymous, nonsynonymous, nonsense and readthrough count 
aggCounts$sumNonsynCount <- rowSums(aggCounts[,c(2,4:5)])
rawNS[match(aggCounts[,1], geneNames),] <- aggCounts[,c(6,3)]
#rawNS[match(aggCounts[,1], geneNames),] <- aggCounts[,c(2,3)]

## (A) elevated substitution rates signal positive selection
# The first approach tests for differences in selection across genes by testing the presence of
# significant different number of substitutions than expected in a gene, given the substitution rate
# estimate for the whole genome. This approach does not distinguish between different types of substitution

## total number of nonsyn and syn substitutions for each gene
nSubs <- rawNS[,1] + rawNS[,2]

names(nSubs) <- rownames(rawNS)
dat <- as.data.frame(cbind(nSubs, geneLengths))

ggplot(dat, aes(x=geneLengths, y=nSubs)) + geom_point() + geom_smooth(method="lm")

## detect outlier
#which(nSubs > 20)
#gbk[which(gbk$synonym == "B0X74_06575"),] # ftsK DNA translocase

which(nSubs == 21) ## this is not too significant as gene length is big
gbk[which(gbk$synonym == "B0X74_24880"),] # siiE large repetive protein

which(nSubs == 19)
gbk[which(gbk$synonym == "B0X74_05735"),] ## oadA oxaloacetate decarboxylase subunit A

##
#nSubs.cut <- nSubs[!names(nSubs) %in% c("B0X74_06575", "B0X74_05735")]
#geneLengths.cut <- geneLengths[!names(geneLengths) %in% c("B0X74_06575", "B0X74_05735")]
cor(nSubs, geneLengths,method = 'pearson') ## high correlation at 0.4886
#cor(nSubs.cut, geneLengths.cut, method='pearson') ## correlation at 0.499
# showing the number of mutations is dependent on gene lengths

lnfit <- lm(nSubs ~ geneLengths)
summary(lnfit) ## adjusted R squared of 0.2386

#lnfit.cut <- lm(nSubs.cut ~ geneLengths.cut)
#summary(lnfit.cut) ## adjusted R-squared of 0.2497

#subRate <- sum(nSubs.cut)/sum(geneLengths.cut) 
subRate <- sum(nSubs)/sum(geneLengths)
subRate ## rate of 0.0009371538 ## rate of substitution mutation.

## calculate the p-value for each gene, this is for the total number of mutations
subRatePVal = sapply(1:length(nSubs),function(i){
      t = poisson.test(nSubs[i],r=(subRate*geneLengths[i]))
      return(t$p.value)
})
subRatePadj <- p.adjust(subRatePVal, method="BH")
subRateSignif = nSubs[which(subRatePadj<0.1)]
subRateSignif

gbk[gbk$synonym %in% names(subRateSignif),]
## The genes having elevated mutation rate include:
# ramR, glycosyltransferase 1, GT2, oadA, phoQ, cspC, fliC,transcriptional regulator, rpoS, btuB
##
# interestingly bapA has 0 mutations despite of its large size. It is possibly that bapA is 
# negatively selected. BapA plays role in biofilm formation.
## genes that are prone to disruptions/deactivation include GT1, cspC, rpoS, btuB

## interpretation:
# rpoS mutation is adaptive to storage condition
# btuB disruption is beneficial in term of colicin desensitization 

#### (B) Estimation of dN/dS 

## load matrix defining types of substitutions
chTypes <- read.table("../Codon_change_types.txt", h=T, sep="\t")
# remove stop codons
notStop = (chTypes[3:66,1])!="STOP"
ct61 = matrix(as.integer(as.matrix(chTypes[3:66,4:67][notStop,notStop])),nrow=61)
ct61[lower.tri(ct61)] = t(ct61)[lower.tri(ct61)]
# Calculate the NY98 substitution rate matrix
makeQ = function(kappa,omega,pi=rep(1/61,61),chTypes) {
      aa = chTypes[3:66,1]
      notStop = aa!="STOP"
      C = matrix(as.integer(as.matrix(chTypes[3:66,4:67][notStop,notStop])),nrow=61)
      Q = matrix(0,nrow=61,ncol=61)
      Q[C==2] = 1
      Q[C==3] = kappa
      Q[C==4] = omega
      Q[C==5] = kappa*omega
      Q[lower.tri(Q)] = t(Q)[lower.tri(Q)]
      Q = t(Q * pi)
      diag(Q) = -apply(Q,1,sum)
      Q = Q/sum(-pi*diag(Q))
      return(Q)
}
kappa = 5.2
dndsNull <- 1.0
codonCount = as.vector(codonUsage[unlist(geneticCode)!="STO"]) # Exclude stop codons
codonFreq = codonCount/sum(codonCount) # Normalize
# Relative frequency of each substitution
ny98 = makeQ(kappa,dndsNull,codonFreq,chTypes)*codonFreq
ny98
# Define the expected rate of various types of substitutions
dsTv = sum(ny98[!is.na(ct61) & ct61==2]) # Synonymous transversions
dsTs = sum(ny98[!is.na(ct61) & ct61==3]) # Synonymous transitions
dnTv = sum(ny98[!is.na(ct61) & ct61==4]) # Non-synonymous transversions
dnTs = sum(ny98[!is.na(ct61) & ct61==5]) # Non-synonymous transitions

# Compute the relative frequency of non-synonymous vs synonymous
# Substitutions under the null hypothesis
(r0 = (dnTv+dnTs)/(dsTv+dsTs))

## This show that under the null hypothesis, we expect to see many more non-synonymous than synonymous substitutions. 
#Since many of the counts of non-synonymous and synonymous substitutions for each genes are zero, there is not sufficient power to estimate dN/dS for each gene.

rawNS2 <- data.frame(matrix(0, nrow=length(geneNames), ncol=2))
rownames(rawNS2) <- geneNames
colnames(rawNS2) <- c("N", "S")
rawNS2[match(aggCounts[,1], geneNames),] <- aggCounts[,c(2,3)]
## output only nonsynonymous (excluding stop and start codon changes)

subsPergene <- cbind(rawNS2, rowSums(rawNS2))
table(t(subsPergene[,3]))
colnames(subsPergene) <- c("N", "S", "nSubs")

#Group genes by the total number of substitutions
aggCountsSubs <- aggregate(subsPergene, list(subsPergene$nSubs), sum)
head(aggCountsSubs)

## For ones that have nonzero no. of syn mutations
rawNS2.cut <- rawNS2[rawNS2[,"S"]!=0,]
ow <- (rawNS2.cut[,"N"]/rawNS2.cut[,"S"])/(r0)
sel1 <- rawNS2.cut[which(ow > 1.5),]
sel1
gbk[gbk$synonym %in% rownames(sel1),]

### Don't use dN/dS on the genes with 0 synonymous mutations, 
# extract out genes with high nonsyn mutations (N>4) and 0 syn mutations, 
sel2 <- rawNS2[which(rawNS2$S==0 & rawNS2$N>4),]
gbk[gbk$synonym %in% rownames(sel2),]

negsel1 <- rawNS[which(rawNS$N==0 & rawNS$S > 4),]
gbk[gbk$synonym %in% rownames(negsel1),]
## output the genes which are identified as positively selected 

positive.sel <- rbind(sel1, sel2)
P <- gbk[gbk$synonym %in% rownames(positive.sel),c("gene", "start", "end", "strand", "length",
                                              "product", "synonym")]
positive.sel <- positive.sel[match(P$synonym, rownames(positive.sel)),]
P <- cbind(positive.sel, P)
write.table(P, file="./ST34_454str_positive_selection_genes.txt", quote=FALSE, sep="\t",
            row.names = TRUE, col.names = TRUE)
##
## Genes with significant number of deactivation mutations
table(aggCounts$NonsenseCount)
pseudogenes <- aggCounts[which(aggCounts$NonsenseCount > 2),1]
L <- gbk[gbk$synonym %in% pseudogenes,c("gene", "start", "end", "strand", "length", "product", "synonym")]
aggCounts[which(aggCounts$NonsenseCount > 2),]
L <- cbind(L,aggCounts[which(aggCounts$NonsenseCount > 2), "NonsenseCount"] )

## The results from this new tree is not different compared to the one from midpoint/beast-rooted tree.
## The same genes returned to be under positive selection. 