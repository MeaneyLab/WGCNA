##########################################################################################
#WGCNA CODE
##########################################################################################

library('WGCNA')
options(stringsAsFactors = FALSE)

femData = read.csv("pDG_Sing_transcript_rpkms_log2.csv") 
datExprSing=as.data.frame(t(femData[, -c(1)]))
names(datExprSing)=femData$gene_id
rownames(datExprSing)=names(femData)[-c(1)]

gsg = goodSamplesGenes(datExprSing, verbose = 3);
gsg$allOK

##########################################################################################

traitData = read.csv("trait_pDG.csv")
allTraits=traitData[,c(1, 3)] 
names(allTraits) 
Rat=rownames(datExprSing)
traitRows = match(Rat, allTraits$Sample)
datTraits = allTraits[traitRows, ]
rownames(datTraits) = allTraits$Sample 

##########################################################################################

A=adjacency(t(datExprSing),type="distance")
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k) 
thresholdZ.k=-5 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
traitColors= data.frame(numbers2colors(datTraits$Condition,signed=FALSE))
dimnames(traitColors)[[2]]= 'Condition'
datColors=data.frame(outlierC=outlierColor,traitColors)
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap") 

##########################################################################################

powers=c(1:20) 
sft=pickSoftThreshold(datExprSing,powerVector=powers, networkType = "signed")
par(mfrow=c(1,1))
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.90,col="red") 
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red") 

##########################################################################################

options(stringsAsFactors = FALSE);
ADJ1=abs(cor(datExprSing,use="p"))^6
k=softConnectivity(datE=datExprSing,power=20)
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

##########################################################################################

bwnet = blockwiseModules(datExprSing,corType="pearson",
                         maxBlockSize=7000,networkType="signed",power=20,minModuleSize=30,
                         mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                         pamRespectsDendro=FALSE,saveTOMFileBase="femaleSingMTOM-blockwise",verbose=3) 

moduleColorsBlockwise=labels2colors(moduleLabelsBlockwise) 
moduleColorsFemale=moduleColorsAutomatic
nGenes = ncol(datExprSing)
nSamples = nrow(datExprSing)
MEs0 = moduleEigengenes(datExprSing,moduleColorsFemale)$eigengenes
MEsFemale = orderMEs(MEs0)
Condition<- datTraits$Condition
modTraitCor = cor(MEsFemale,Condition, use = "p")
modTraitP = round(corPvalueStudent(modTraitCor, nSamples),digits = 2)

##########################################################################################

datKME=signedKME(datExprSing, MEsFemale) 

GeneAnnotation=read.delim(file="annotated_pDG_C_vs_B_S_de_data.txt")
probes = names(datExprSing)
probes2annot = match(probes,GeneAnnotation$EnsemblID)
datGS.Traits=data.frame(cor(datExprSing,Condition,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datOutput=data.frame(ProbeID=names(datExprSing),
                     GeneAnnotation[probes2annot,],moduleColorsFemale,datKME,datGS.Traits)
write.table(datOutput,"Results.csv",row.names=F,sep=",") 

##########################################################################################
# Module preservation
##########################################################################################

VietnamData = read.csv("pDG_V_transcript_rpkms_log2.csv")
names(VietnamData)
datExprViet=data.frame(t(VietnamData[,-c(1)]))
names(datExprViet)=VietnamData$gene_id
gsg = goodSamplesGenes(datExprViet, verbose = 3);
gsg$allOK
setLabels = c("Singapore", "Vietnam")
multiExpr=list(Singapore=list(data=datExprSing),
               Vietnam=list(data=datExprViet))
moduleColorsFemale=moduleColorsAutomatic
multiColor=list(Singapore=moduleColorsFemale) 
nPermutations1=200
set.seed(1)
system.time({mp = modulePreservation(multiExpr, multiColor,
                  referenceNetworks = 1, nPermutations = nPermutations1,
                  randomSeed = 1, quickCor = 0, verbose = 3)})
save(mp, file = "modulePreservation_pDG.RData")
ref=1; test = 2 
Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
Obs.PreservationStats 
Z.PreservationStats 


