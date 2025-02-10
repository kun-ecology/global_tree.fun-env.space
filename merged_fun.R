# functions  from Carmona, Carlos P. et al. "Erosion of global functional diversity across the tree of life." Science Advances 7, no. 13 (2021): eabf2675.

 # start (00_Aux_graphics.R)
library(RColorBrewer)
 
alphaCol<-60
colorsPick <- "RdYlGn"
colThreat <- brewer.pal(11,colorsPick)[c(9,3)]
colThreatLines <- brewer.pal(11,colorsPick)[c(11,1)]
colWorld <- "black" 
colGradient <- c("white", "white", "yellow", "red")

for (i in 1:length(colThreat)){
  colThreat[i] <- rgb(t(col2rgb(colThreat)/255), alpha=alphaCol/255)[i]
}
### CONTOUR LEVELS:
conLines<-c(50, 99)
lineTypeWorld <- 1
lineTypeContour <- 1
lineTypeWorldExt <- 1
thickCountour <- 2+(length(conLines):1)
thickWorldExt <- max(thickCountour)

ptSize<-30
cexMain <- 1.25

### PCA Arrows:
multArrow <- c(0.85, 1.5, 2.5, 1.7, 0.85, 1)
multiplierTextGr <- c(1.2, 1.2, 1.25, 1.15, 1.3, 1)
### legend position:
legPos <- c("bottomleft", "bottomright", "topright", "bottomleft", "bottomleft", "topleft")
linesLegPos <- c("bottomright", "bottomleft", "topleft", "bottomright", "bottomright", "topright")


limXlist <- list( c(-5.7, 5.5), #Plants
                  c(-4.6, 8), #Mammals
                  c(-4.7, 10.1), #Aves
                  c(-5, 7), #reptiles
                  c(-5, 5.9),#Amphibians
                  c(-7, 7)) #FWFish
#
limYlist <- list( c(-8.5, 6.5), #Plants
                  c(-3.5, 4), #Mammals
                  c(-3, 6.1), #Aves
                  c(-4, 3.6), #reptiles
                  c(-5.5,5),
                  c(-7, 7)) #FWFish


 # end (00_Aux_graphics.R) 


 # start (00_PCA_functions.R)

 
sorted.data.for.pca<-function(TraitsImputed = read.table("dataPrepared/Plants/plantTraitsImputedIUCN.txt"),
                              TraitsMissing = read.table("dataPrepared/Plants/plantTraitsMissingIUCN.txt")){
  TraitsImputed[, 1:(ncol(TraitsImputed)-1)] <- log10(TraitsImputed[, 1:(ncol(TraitsImputed)-1)]+1)
  TraitsImputed[, 1:(ncol(TraitsImputed)-1)] <- scale(TraitsImputed[, 1:(ncol(TraitsImputed)-1)])
  TraitsMissing[, 1:(ncol(TraitsImputed)-1)] <- replace(TraitsMissing[, 1:(ncol(TraitsImputed)-1)],TraitsMissing[, 1:(ncol(TraitsImputed)-1)]> -1000, 1) * TraitsImputed[, 1:(ncol(TraitsImputed)-1)]
  removeSP <- which(rowSums(!is.na(TraitsMissing[, 1:(ncol(TraitsImputed)-1)]))==0)
  TraitsMissing <- TraitsMissing[- removeSP, ]
  TraitsImputed <- TraitsImputed[- removeSP, ]
  return(list(TraitsImputed=TraitsImputed,TraitsMissing=TraitsMissing))
}
make.PCA<-function(traitsImputed,traitsMissing,group=c("Plants","Birds","Mammals","Reptiles","Amphibians","Fishes")){
  require(paran)
  paranAux <- paran(traitsImputed[, 1:(ncol(traitsImputed) - 1)])
  if(paranAux$Retained<2){paranAux$Retained=2}
  dimensionsAux <- paranAux$Retained
  PCAImpute <- list()
  PCAImpute$traits <- traitsImputed[, 1:(ncol(traitsImputed) - 1)]
  PCAImpute$traitsMissing <- traitsMissing[, 1:(ncol(traitsImputed) - 1)]
  PCAImpute$nImputed <- rowSums(is.na(PCAImpute$traitsMissing))
  PCAImpute$PCA <- princomp(traitsImputed[, 1:(ncol(traitsImputed) - 1)])
  PCAImpute$dimensions <- dimensionsAux
  PCAImpute$variance <- (summary(PCAImpute$PCA)[1][[1]]^2)[1:dimensionsAux] /
    sum(summary(PCAImpute$PCA)[1][[1]]^2)
  if(group %in% c("Birds","Mammals") ) { #Mammals & Birds
    if(PCAImpute$PCA$loadings["bm", "Comp.1"]>0){
      multPC1 <- 1
    } else{
      multPC1 <- -1
    }
    PCAImpute$PCA$scores[,1] <- multPC1 * PCAImpute$PCA$scores[,1]
    PCAImpute$PCA$loadings[, "Comp.1"] <- multPC1 * PCAImpute$PCA$loadings[, "Comp.1"]
  }
  if(group == "Reptiles") { #Reptiles
    if(PCAImpute$PCA$loadings["inc", "Comp.2"]>0){
      multPC2 <- -1 
    } else{
      multPC2 <- 1
    }
    PCAImpute$PCA$scores[,2] <- multPC2 * PCAImpute$PCA$scores[,2]
    PCAImpute$PCA$loadings[, "Comp.2"] <- multPC2 * PCAImpute$PCA$loadings[, "Comp.2"]
  }
  
  if(group == "Fishes") { #Fishes
    if(PCAImpute$PCA$loadings["svl", "Comp.1"]>0){
      multPC1 <- 1
    } else{
      multPC1 <- -1
    }
    if(PCAImpute$PCA$loadings["svl", "Comp.2"]<0){
      multPC2 <- 1
    } else{
      multPC2 <- -1
    }
    PCAImpute$PCA$scores[,1] <- multPC1 * PCAImpute$PCA$scores[,1]
    PCAImpute$PCA$loadings[, "Comp.1"] <- multPC1 * PCAImpute$PCA$loadings[, "Comp.1"]
    PCAImpute$PCA$scores[,2] <- multPC2 * PCAImpute$PCA$scores[,2]
    PCAImpute$PCA$loadings[, "Comp.2"] <- multPC2 * PCAImpute$PCA$loadings[, "Comp.2"]
  }
  PCAImpute$loadings <- PCAImpute$PCA$loadings
  PCAImpute$traitsUse <- data.frame(PCAImpute$PCA$scores[, 1:PCAImpute$dimensions]) 
  PCAImpute$IUCN <- factor(plantTraitsImputed$IUCN)
  names(PCAImpute$IUCN)<-rownames(plantTraitsImputed)
  treathCat <- c("EX", "EW", "CR", "EN", "VU")
  noThreatCat <- c("LC", "NT")
  PCAImpute$Threat <- with(traitsImputed, ifelse(IUCN %in% treathCat, 1, 
                                                                  ifelse(IUCN %in% noThreatCat, 0, NA)))
  return(PCAImpute)
  
}

make.TPD.trad<-function(targetGroupName,PCAImpute,alphaUse=0.95,
                        saveFile= paste0(getwd(),"/TPDsTradi.rds")){
  require(TPD)
  cat(paste("\n ESTIMATING TPDs for: ", targetGroupName, "\n"))
  data <- PCAImpute
  dimensions <- PCAImpute$dimensions
  traitsUSE <- PCAImpute$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:2)
  gridSize <- ifelse(dimensions == 4, 30, 100)
  sdTraits <- sqrt(diag(Hpi.diag(traitsUSE)))
  TPDsAux <- TPDsMean(species = rownames(traitsUSE), 
                      means = traitsUSE, 
                      sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                      alpha = alphaUse,
                      n_divisions = gridSize)
  saveRDS(TPDsAux,saveFile)
  TPDsAux <-NULL
  gc()
}

make.TPD.2D.high.def<-function(targetGroupName,PCAImpute,alphaUse=0.95,
                        saveFile= paste0(getwd(),"/TPDs2DHighDef.rds")){
  require(TPD)
  cat(paste("\n ESTIMATING TPDs for: ", targetGroupName, "\n"))
  data <- PCAImpute
  dimensions <- 2
  traitsUSE <- PCAImpute$traitsUse[, 1:dimensions]
  colnames(traitsUSE)<-paste0("Comp.",1:2)
  gridSize <- ifelse(dimensions == 4, 30, 200)
  sdTraits <- sqrt(diag(Hpi.diag(traitsUSE)))
  TPDsAux <- TPDsMean(species = rownames(traitsUSE), 
                      means = traitsUSE, 
                      sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                      alpha = alphaUse,
                      n_divisions = gridSize)
  saveRDS(TPDsAux, saveFile)
  TPDsAux <-NULL
  gc()
}

make.TPD.2D.low.def<-function(targetGroupName,PCAImpute,alphaUse=0.95,
                           saveFile= paste0(getwd(),"/TPDs2DLowDef.rds")){
  require(TPD)
  cat(paste("\n ESTIMATING TPDs for: ", targetGroupName, "\n"))
  data <- PCAImpute
  dimensions <- 2
  traitsUSE <- PCAImpute$traitsUse[, 1:dimensions]
  colnames(traitsUSE)<-paste0("Comp.",1:2)
  gridSize <- ifelse(dimensions == 4, 30, 100)
  sdTraits <- sqrt(diag(Hpi.diag(traitsUSE)))
  TPDsAux <- TPDsMean(species = rownames(traitsUSE), 
                      means = traitsUSE, 
                      sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                      alpha = alphaUse,
                      n_divisions = gridSize)
  saveRDS(TPDsAux,saveFile)
  TPDsAux <-NULL
  gc()  
}

make.TPDLarge<-function(targetGroupName,PCAImpute,alphaUse=0.95,TPDsMean_large,
                              saveFile= paste0(getwd(),"/TPDs2DLowDef.rds")){
  cat(paste("\n ESTIMATING TPDs for: ", targetGroupName, "\n"))
  data <-PCAImpute
  dimensions <- PCAImpute$dimensions
  traitsUSE <-  PCAImpute$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  gridSize <- ifelse(dimensions == 4, 30, 100)
  sdTraits <- sqrt(diag(Hpi.diag(traitsUSE)))
  TPDsAux <- TPDsMean_large(species = rownames(traitsUSE), 
                            means = traitsUSE, 
                            sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                            alpha = alphaUse,
                            n_divisions = gridSize)
  saveRDS(TPDsAux, saveFile)
  TPDsAux <-NULL
  gc()
}

make.functional.space<-function(species,traitsUSE,TPDsAux,PCAImpute,Rank,taxo,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  var<-NULL
  var<-PCAImpute$variance
  if(is.null(var)){var<-PCAImpute$Variance}
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * var[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * var[2], 2), "%)")
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]


  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1,cex.lab=1)
  graphics::box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  
  
  if (Rank == "all"){
    contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
            drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
    
    contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
            drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
  }

  if(Rank != "all"){
    contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.998,
            drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
    
    nbGr<-length(namesOrd)
    ## for the world
    colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
    gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
    ColorordR <- rev(gradientColorsF(nbGr))
    subGBIF<-taxo
    j<-1
    for (ordR in namesOrd){
      commNames <- c("All")
      
      eval(parse(text=paste0("commMatrix <- matrix(0, nrow = length(commNames), ncol = length(rownames(subGBIF[subGBIF$",Rank,"==ordR,])),
                         dimnames = list(commNames, rownames(subGBIF[subGBIF$",Rank,"==ordR,])))")))
      
      
      
      commMatrix["All", ] <- 1
      commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop=F]
      
      TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
      imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
      
      trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
      trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
      
      contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
              drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=subGBIF[j], add=T)
      j=j+1
    }
    
  }
  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
}

make.functional.space.arrows<-function(species,traitsUSE,TPDsAux,PCAImpute,Rank,taxo,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA,dim=1){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  colAbove <- "darkolivegreen"
  colBelow <- "darkorange4"
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  var<-NULL
  var<-PCAImpute$variance
  if(is.null(var)){var<-PCAImpute$Variance}
  if(dim==1){
    var<-var[c(1,2)]
    }
  if(dim==6){
    var<-var[c(3,4)]
    }
  
  
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * var[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * var[2], 2), "%)")
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
  
  
  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1,cex.lab=1)
  graphics::box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  
    contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
            drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
    
    contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
            drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
  
  
  if(dim == 1) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 2)]
  if(dim == 2) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 3)]
  if(dim == 3) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 4)]
  if(dim == 4) loadingsUse <- PCAImpute$PCA$loadings[, c(2, 3)]
  if(dim == 5) loadingsUse <- PCAImpute$PCA$loadings[, c(2, 4)]
  if(dim == 6) loadingsUse <- PCAImpute$PCA$loadings[, c(3, 4)]
  
  multLoadings <- numeric()
  for(j in 1:nrow(loadingsUse)){
    multLoadings[j] <- dist(rbind(c(0,0), loadingsUse[j,1:2]))
  }
  multArrow <- 0.8
  multiplierTextGr <- 1.2
  maxDist <- min(abs(c(limY, limX)))
  offChar <- c(0.025 * maxDist, 0.05 * maxDist, 0.04 * maxDist, 
               0.025 * maxDist, 0.02 * maxDist, 0.035 *maxDist)
  multiplierArrow <- multArrow * maxDist
  multiplierText <- multiplierTextGr * multiplierArrow
  multOther <- -1.3
  minLoadPlot <- 0.35
  colAllTraits <- c(rep(colAbove, 6), rep(colBelow, 4))
  
  for(j in 1:nrow(loadingsUse)){
    if(any(abs(loadingsUse[j,]) > minLoadPlot)){
      radianAng <- base::atan2(loadingsUse[j,2], loadingsUse[j,1])
      txtAng <-  radianAng * (180 / pi)
      txtAng <- ifelse(txtAng < 0, txtAng + 360, txtAng)
      offsetText <- offChar[3] * nchar(rownames(loadingsUse)[j])
      if(txtAng > 92 & txtAng < 270){
        txtAng <- txtAng + 180
      }
      Arrows(x0=0, y0= 0, x1= multiplierArrow*loadingsUse[j,1], 
             y1= multiplierArrow*loadingsUse[j,2], 
             col=colAllTraits[j], lwd=2, arr.type="triangle", arr.length = 30/70)
      adjVal <- c(0.5, 0.5)
      posX <- multiplierArrow*loadingsUse[j,1] + (0.1*maxDist + offsetText) * cos(radianAng)
      posY <- multiplierArrow*loadingsUse[j,2] + (0.1*maxDist + offsetText) * sin(radianAng)
      text(x= posX, y= posY, labels=rownames(loadingsUse)[j], col=colAllTraits[j], cex=1.25,
           offset = 0, adj = adjVal, srt = txtAng) 
    }
  }
  
  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
}

make.functional.space.trait<-function(species,traitsUSE,TPDsA,PCAImpute,taxo,taxon,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMatALL <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1EdgesALL <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2EdgesALL <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMatALL[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMatALL[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]

  if (taxo %in% rownames(traitsUSE)){
    image(x=trait1EdgesALL, y=trait2EdgesALL, z=imageMatALL[, , "All"], col=ColorRamp_ex,
          xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
          main = "", asp = 1,cex.lab=1)
    graphics::box(which="plot")
    axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
    axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
    contour(x=trait1EdgesALL, y=trait2EdgesALL, z=imageMatALL[, , "All"], levels=0.99,
            drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
    
    contour(x = trait1EdgesALL, y = trait2EdgesALL, z = imageMatALL[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
            drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
    
    points(traitsUSE[taxo,1],traitsUSE[taxo,1],cex=3,col=1:length(taxo),pch=16)
    mtext(text=paste0("Species: ",taxo),side=3,  cex = 2, line=1.5) 
  }
  
  if(taxo %in% taxon$order){
    nbGr<-length(taxo)
    ## for the world
    colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
    gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
    ColorordR <- rev(gradientColorsF(nbGr))
    j<-1
    for (ordR in taxo){
      commNames <- c("All")
      
      commMatrix <- matrix(0, nrow = length(commNames), ncol = length(rownames(taxon[taxon$order==ordR,])),
                         dimnames = list(commNames, rownames(taxon[taxon$order==ordR,])))
      
      
      
      commMatrix["All", ] <- 1
      commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop=F]
      
      TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
      imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
      
      trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
      trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
      
      image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
            xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
            main = "", asp = 1,cex.lab=1)
      graphics::box(which="plot")
      axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      
      contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
              drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=ColorordR[j], add=T)
      
      contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
              drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
      
      j=j+1
      mtext(text=paste0("Order: ",taxo),side=3,  cex = 2, line=1.5) 
    }
    
  }
  
  if(taxo %in% taxon$family){
    nbGr<-length(taxo)
    ## for the world
    colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
    gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
    ColorordR <- rev(gradientColorsF(nbGr))
    j<-1
    for (ordR in taxo){
      commNames <- c("All")
      
      commMatrix <- matrix(0, nrow = length(commNames), ncol = length(rownames(taxon[taxon$family==ordR,])),
                           dimnames = list(commNames, rownames(taxon[taxon$family==ordR,])))
      
      
      
      commMatrix["All", ] <- 1
      commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop=F]
      
      TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
      imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
      
      trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
      trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
      
      image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
            xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
            main = "", asp = 1,cex.lab=1)
      graphics::box(which="plot")
      axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      
      contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
              drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=ColorordR[j], add=T)
      
      contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
              drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
      
      j=j+1
      mtext(text=paste0("Family: ",taxo),side=3,  cex = 2, line=1.5) 
    }
    
  }
  
  if(taxo %in% taxon$genus){
    nbGr<-length(taxo)
    ## for the world
    colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
    gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
    ColorordR <- rev(gradientColorsF(nbGr))
    j<-1
    for (ordR in taxo){
      commNames <- c("All")
      
      commMatrix <- matrix(0, nrow = length(commNames), ncol = length(rownames(taxon[taxon$genus==ordR,])),
                           dimnames = list(commNames, rownames(taxon[taxon$genus==ordR,])))
      
      
      
      commMatrix["All", ] <- 1
      commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop=F]
      
      TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
      imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
      
      trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
      trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
      
      image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
            xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
            main = "", asp = 1,cex.lab=1)
      graphics::box(which="plot")
      axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
      
      contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
              drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=ColorordR[j], add=T)
      
      contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
              drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
      
      j=j+1
    }
      mtext(text=paste0("Genus: ",taxo),side=3,  cex = 2, line=1.5)  
  }

}

make.functional.space.trait.up<-function(traitsUSE,TPDsA,PCAImpute,taxo,imageTPD,
                                          ncolors,limX,limY,ColorRamp,title=NA){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  dataAnalysis<-traitsUSE#[rownames(traitsUSE)%in%species,]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMatALL <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1EdgesALL <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2EdgesALL <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMatALL[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMatALL[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
  commNames <- c("All")
  
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(taxo),
                       dimnames = list(commNames, rownames(taxo)))

  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop=F]
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  
  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1,cex.lab=1)
  graphics::box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  contour(x=trait1EdgesALL, y=trait2EdgesALL, z=imageMatALL[, , "All"], levels=0.99,
          drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
  contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
          drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=ColorRamp, add=T)
  contour(x = trait1Edges, y = trait2Edges, z = imageMat[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
          drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
}

make.functional.space.realms<-function(species,traitsUSE,TPDsA,PCAImpute,Rank,taxo,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA,taxRealm,Realm="Afrotropical"){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop = F]

  
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMatALL <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1EdgesAll <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2EdgesAll <- unique(TPDcAux$data$evaluation_grid[,2])
      
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%gsub(" ","_",taxRealm[[Realm]]$speciesGbifTraitRealms),]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop = F]
  commMatrixALL<-commMatrix
  
  TPDcAuxRealm <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  imageMatRealm <-imageTPD(TPDcAuxRealm, thresholdPlot = 1)
  
  trait1EdgesRealm <- unique(TPDcAuxRealm$data$evaluation_grid[,1])
  trait2EdgesRealm <- unique(TPDcAuxRealm$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMatRealm[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMatRealm[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
  ylabY<-ylab
  image(x=trait1EdgesRealm, y=trait2EdgesRealm, z=imageMatRealm[, , "All"], col=ColorRamp_ex,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylabY, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1,cex.lab=1)
  graphics::box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  
  contour(x=trait1EdgesAll, y=trait2EdgesAll, z=imageMatALL[, , "All"], levels=0.998,
          drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
  # 
  if (Rank == "all"){
    contour(x = trait1EdgesRealm, y = trait2EdgesRealm, z = imageMatRealm[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
            drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
    contour(x=trait1EdgesRealm, y=trait2EdgesRealm, z=imageMatRealm[, , "All"], levels=0.99,
            drawlabels = F, labcex = 0.8, lwd=2.3, lty=1, col="red", add=T)
  }
  
  if(Rank != "all"){
    nbGr<-length(namesOrd)
    ## for the world
    colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
    gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
    ColorordR <- rev(gradientColorsF(nbGr))
    subGBIF<-taxo
    j<-1
    for (ordR in namesOrd){
      commNames <- c("All")
      commMatrix<-commMatrixALL
      eval(parse(text=paste0("commMatrix['All',colnames(commMatrix)%in%rownames(subGBIF[subGBIF$",Rank,"==ordR,])] <- 1")))
      commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop=F]
      
      TPDcAuxRealm <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
      imageMat <-imageTPD(TPDcAuxRealm, thresholdPlot = 1)
      
      trait1Edges <- unique(TPDcAuxRealm$data$evaluation_grid[,1])
      trait2Edges <- unique(TPDcAuxRealm$data$evaluation_grid[,2])
      
      contour(x = trait1Edges, y = trait2Edges, z = imageMatRealm[, , "All"], levels=c(seq(0.5, 0.99, by=0.2)),
              drawlabels = F, labcex = 0.8, lwd=0.01, lty=1, col="grey75", add=T)
      contour(x=trait1Edges, y=trait2Edges, z=imageMatRealm[, , "All"], levels=0.99,
              drawlabels = F, labcex = 0.8, lwd=2.3, lty=1, col="red", add=T)
      j=j+1
    }
    
  }
  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
}

make.functional.space.loss<-function(species,traitsUSE,TPDsA,PCAImpute,Rank,taxo,namesOrd,imageTPD,
                                      ncolors,limX,limY,ColorRamp,title=NA,IUCNStatus){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  
  species<-species[species%in%TPDsA$data$species]

  names(IUCNStatus)<-names(PCAImpute$Threat)<-rownames(PCAImpute$traitsUse)
  IUCNStatus<-names(na.omit(IUCNStatus[names(IUCNStatus)%in%species]))
  dataAnalysis <- data.frame(traitsUSE[IUCNStatus,],
                             Threat = PCAImpute$Threat[IUCNStatus])
  commNames <- c("ALL", "IUCN", "WithoutThreatened")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["ALL", ] <- 1
  commMatrix["IUCN", rownames(dataAnalysis)] <- 1
  commMatrix["WithoutThreatened", 
             rownames(subset(dataAnalysis, dataAnalysis$Threat == 0))] <- 1
  
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 *PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMatALL <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"ALL"]
  imageMatIUCN <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"IUCN"]
  imageMatRemain <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"WithoutThreatened"]
  imageMatIUCN1 <-imageTPD(TPDcAux, thresholdPlot = 1)[,,"IUCN"]
  imageMatRemain1 <-imageTPD(TPDcAux, thresholdPlot = 1)[,,"WithoutThreatened"]
  
  imageMatChange <- imageMatIUCN - imageMatRemain
  imageMatLostSpace <- replace(imageMatIUCN, imageMatIUCN>0, NA)
  imageMatLostSpace <- replace(imageMatLostSpace, 
                               imageMatIUCN>0 & is.na(imageMatRemain), 1)
  Comp.1Vec <- unique(TPDcAux$data$evaluation_grid[,1])
  Comp.2Vec <- unique(TPDcAux$data$evaluation_grid[,2])

  ncol <- 1000
  colorsPick <- "vik"
  ColorRamp <-  rev(scico(n=ncol, palette=colorsPick))
  nHalf <- sum(!is.na(imageMatChange))/2
  Min <- -0.36
  Max <- 0.29
  Thresh <- 0
  ## Make vector of colors for values below threshold
  rc1 <- colorRampPalette(colors = c(ColorRamp[1:500]), space="Lab")(nHalf)    
  ## Make vector of colors for values above threshold
  rc2 <- colorRampPalette(colors = c(ColorRamp[501:1000]), space="Lab")(nHalf)
  rampcols <- c(rc1, rc2)
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  image(x=Comp.1Vec, y=Comp.2Vec, z=imageMatChange, xlim=limX, ylim=limY, 
        col=rampcols, breaks=rampbreaks, xaxs="r", yaxs="r", axes=F, asp=1, 
        xlab=paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)"), 
        ylab=paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)"))
  contIUCN99<-contourLines(x = Comp.1Vec, y = Comp.2Vec, z = imageMatIUCN1, levels= c(0.99))
  contRemain99<-contourLines(x = Comp.1Vec, y = Comp.2Vec, z = imageMatRemain1, levels= c(0.99))
  for(fig in 1:length(contIUCN99)){
    polygon(x=c(contIUCN99[[fig]]$x), y=c(contIUCN99[[fig]]$y), col=1)
  } 
  image(x=Comp.1Vec, y=Comp.2Vec, z=imageMatChange, xlim=limX, ylim=limY, 
        col=rampcols, breaks=rampbreaks, add=T)
  
  graphics::box(which="plot")
  axis(1,tcl=0.3,lwd=0.8)
  axis(2, las=1, tcl=0.3,lwd=0.8)
  
  
  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
  
}

make.functional.space.realm.loss<-function(species,traitsUSE,TPDsA,PCAImpute,Rank,taxo,namesOrd,imageTPD,
                                     ncolors,limX,limY,ColorRamp,title=NA,IUCNStatus,taxRealm,Realm="Afrotropical"){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  
  species<-species[species%in%gsub(" ","_",taxRealm[[Realm]]$speciesGbifTraitRealms)]
  species<-species[species%in%TPDsA$data$species]
  
  names(IUCNStatus)<-names(PCAImpute$Threat)<-rownames(PCAImpute$traitsUse)
  IUCNStatus<-names(na.omit(IUCNStatus[names(IUCNStatus)%in%species]))
  dataAnalysis <- data.frame(traitsUSE[IUCNStatus,],
                             Threat = PCAImpute$Threat[IUCNStatus])
  commNames <- c("ALL", "IUCN", "WithoutThreatened")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["ALL", ] <- 1
  commMatrix["IUCN", rownames(dataAnalysis)] <- 1
  commMatrix["WithoutThreatened", 
             rownames(subset(dataAnalysis, dataAnalysis$Threat == 0))] <- 1
  
  TPDcAux <- TPDc(TPDs = TPDsA, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 *PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMatALL <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"ALL"]
  imageMatIUCN <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"IUCN"]
  imageMatRemain <-imageTPD(TPDcAux, thresholdPlot = 0.99)[,,"WithoutThreatened"]
  imageMatIUCN1 <-imageTPD(TPDcAux, thresholdPlot = 1)[,,"IUCN"]
  imageMatRemain1 <-imageTPD(TPDcAux, thresholdPlot = 1)[,,"WithoutThreatened"]
  
  imageMatChange <- imageMatIUCN - imageMatRemain
  imageMatLostSpace <- replace(imageMatIUCN, imageMatIUCN>0, NA)
  imageMatLostSpace <- replace(imageMatLostSpace, 
                               imageMatIUCN>0 & is.na(imageMatRemain), 1)
  Comp.1Vec <- unique(TPDcAux$data$evaluation_grid[,1])
  Comp.2Vec <- unique(TPDcAux$data$evaluation_grid[,2])
  
  ncol <- 1000
  colorsPick <- "vik"
  ColorRamp <-  rev(scico(n=ncol, palette=colorsPick))
  nHalf <- sum(!is.na(imageMatChange))/2
  Min <- -0.36
  Max <- 0.29
  Thresh <- 0
  ## Make vector of colors for values below threshold
  rc1 <- colorRampPalette(colors = c(ColorRamp[1:500]), space="Lab")(nHalf)    
  ## Make vector of colors for values above threshold
  rc2 <- colorRampPalette(colors = c(ColorRamp[501:1000]), space="Lab")(nHalf)
  rampcols <- c(rc1, rc2)
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  image(x=Comp.1Vec, y=Comp.2Vec, z=imageMatChange, xlim=limX, ylim=limY, 
        col=rampcols, breaks=rampbreaks, xaxs="r", yaxs="r", axes=F, asp=1, 
        xlab=paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)"), 
        ylab=paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)"))
  contIUCN99<-contourLines(x = Comp.1Vec, y = Comp.2Vec, z = imageMatIUCN1, levels= c(0.99))
  contRemain99<-contourLines(x = Comp.1Vec, y = Comp.2Vec, z = imageMatRemain1, levels= c(0.99))
  for(fig in 1:length(contIUCN99)){
    polygon(x=c(contIUCN99[[fig]]$x), y=c(contIUCN99[[fig]]$y), col=1)
  } 
  image(x=Comp.1Vec, y=Comp.2Vec, z=imageMatChange, xlim=limX, ylim=limY, 
        col=rampcols, breaks=rampbreaks, add=T)
  
  graphics::box(which="plot")
  axis(1,tcl=0.3,lwd=0.8)
  axis(2, las=1, tcl=0.3,lwd=0.8)
  
  
  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
  
}

make.functional.space.orders<-function(species,traitsUSE,Rank,taxo,namesOrd,col,TPDsAux,PCAImpute,imageTPD,
                                       ncolors,limX,limY,ColorRamp,title=NA){
  
  nbGr<-length(namesOrd)
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
  
  
  image(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], col=ColorRamp_ex,
        xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
        main = "", asp = 1,cex.lab=1)
  box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  
  
  contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.998,
          drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
  
  ## for the world
  colGradient <- c("royalblue3","darkgreen", "lightsalmon3")
  gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
  ColorordR <- rev(gradientColorsF(nbGr))
  subGBIF<-taxo
  j<-1
  for (ordR in namesOrd){
    commNames <- c("All")
    
    eval(parse(text=paste0("commMatrix <- matrix(0, nrow = length(commNames), ncol = length(rownames(subGBIF[subGBIF$",Rank,"==ordR,])),
                         dimnames = list(commNames, rownames(subGBIF[subGBIF$",Rank,"==ordR,])))")))
    
    
    
    commMatrix["All", ] <- 1
    commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop=F]
    
    TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
    imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
    
    trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
    trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
    
    contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
            drawlabels = F, labcex = 0.8, lwd=2, lty=1, col=col[j], add=T)
    j=j+1
  }
  
  legend("topleft", legend=namesOrd, bty="n",cex=0.9,fill=col, ncol = 2)
  if(!is.na(title)){
    mtext(text=paste0(title,":",Rank),side=3,  cex = 2, line=1.5)  
  }
}

make.functional.diversity.all<-function(species,traitsUSE,TPDsA,target,taxRealm){
  
    
  #if(target==c("Fishes","Plants")){
    # functionTPDRichness <- TPDRichness_large
    # functionTPDc <- TPDc_large
    # functionredundancy <- redundancy_large
  # } else{
  #   functionTPDc <- TPDc
  #   functionTPDRichness <- TPDRichness
  #   functionredundancy <- redundancy
  # }
  biogeo=c("Afrotropical","Australian","Nearctic","Neotropical","Oriental","Palearctic")
  
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  commNames <- c("World",biogeo)
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["World", ] <- 1
  for(i in 1:length(biogeo)){
    commMatrix[biogeo[i], colnames(commMatrix)%in%gsub(" ","_",taxRealm[[i]]$speciesGbifTraitRealms) ] <- 1
  }
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsA$data$species,drop = F]
  
  TPDsC <- TPDc_large(TPDs = TPDsA, sampUnit = commMatrix)
  TPDc_realms <- REND_large(TPDc = TPDsC)
  FRed_realms <- redundancy_large(TPDsC)
  #
  #
  alphaFunc<- t(cbind.data.frame(TD=apply(commMatrix,1,sum),FRic=round((TPDc_realms$communities$FRichness/TPDc_realms$communities$FRichness[1])*100,2),
                                 FEve=round(TPDc_realms$communities$FEvenness,2),
                                 FDiv=round(TPDc_realms$communities$FDivergence,2),
                                 FRed=round(FRed_realms$redundancy),
                                 FRedRelative=round(FRed_realms$redundancyRelative,2)))
  
  alphaFunc
}


make.PCA.phylo<-function(x,y,n,biogeo,PCAImpute,TPDsAux,species,traitsUSE,limX,limY,colorsA,colors,title=NA){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  xlab <- paste0("PC1 (", round(100 * PCAImpute$variance[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * PCAImpute$variance[2], 2), "%)")
  xmin <- 0
  xmax <- 1
  
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  commNames <- c("All")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)

  plot(NA,xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=limX, ylim=limY,
       main = "", asp = 1,cex.lab=1)
  box(which="plot")
  axis(1,tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  axis(2, las=1, tcl=0.5,lwd=0.8,labels=T,cex=0.9)
  
  points(x[1:n],y[1:n],cex=1,pch=16,col=addTrans(colors,80))
  
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , "All"], levels=0.99,
          drawlabels = F, labcex = 0.8, lwd=1.8, lty=2, col="grey50", add=T)
  
  legend("topleft",legend = biogeo,col=colorsA,pch=16,border = "white")

  if(!is.na(title)){
    mtext(text=title,side=3,  cex = 2, line=1.5)  
  }
}

addTrans <- function(color,trans){
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


make.functional.diversity<-function(data,spList,targetGroupName,functionTPDRichness,functionTPDc,
                                    functionredundancy,functiondissim,functiondensityProfileTPD,boot=99,
                                    TPDsAux,namesOrd,Rank,saveFile=paste0(getwd(),"/FD.RDS")){
  
  t0 <- Sys.time()
  dimensions <- data$dimensions
  traitsUSE <- data$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  cat(paste("\n MODEL FOR: ", targetGroupName, " for ",namesOrd,Rank,"\n"))
  alphaFunc<-list()
  
  spList<-spList[rownames(spList)%in%TPDsAux$data$species,]
  eval(parse(text=paste0("spOccurenceObsOrder <- matrix(0,nc=nrow(spList),nr=length(levels(spList$",Rank,")),
                         dimnames = list(levels(spList$",Rank,"),rownames(spList)))")))
  
  eval(parse(text=paste0("lvLRank<-levels(spList$",Rank,")")))
  
  for(or in lvLRank){
    eval(parse(text=paste0("spOccurenceObsOrder['",or,"',colnames(spOccurenceObsOrder)%in%rownames(spList[spList$",Rank," == '",or,"',])]<-1")))
  }
  
  TPDc_realms <- functionTPDc(TPDsAux,spOccurenceObsOrder)
  alphaFunc[[1]]<- functionTPDRichness(TPDc_realms)$communities$FRichness
  rm(TPDc_realms)
  gc()
  
  
  for (bo in 2:boot){
    spOccurenceBoot<-spOccurenceObsOrder
    # for(or in 1:nrow(spOccurenceBoot)){
    #   spOccurenceBoot[or,]<-sample(spOccurenceBoot[or,])
    # }
    colnames(spOccurenceBoot)<-sample(colnames(spOccurenceObsOrder))
    TPDc_realms <- functionTPDc(TPDsAux,spOccurenceBoot)
    alphaFunc[[bo]]<- functionTPDRichness(TPDc_realms)$communities$FRichness
    rm(TPDc_realms)
    gc()
    t1 <- Sys.time()
    cat(paste0("\r GROUP:",targetGroupName," :",bo,"/",boot," .... Progress: ",round(bo/boot,4)*100,"% in ",hms_span(t0, t1),"\r"))
  }
  alphaFunc<-matrix(unlist(alphaFunc),nr=length(lvLRank),ncol=boot,byrow = F,
                    dimnames = list(lvLRank,c("Obs",paste0("Boot ",1:(boot-1)))))
  
  alphaFunc<-cbind(Species=apply(spOccurenceObsOrder,1,sum),alphaFunc)
  
  if(namesOrd=="all"){namesOrd<-lvLRank}
  alphaFunc<-alphaFunc[namesOrd,]
  saveRDS(alphaFunc,saveFile) 
  return(alphaFunc)
}

make.functional.diversity.loss<-function(data,spList,focal,targetGroupName,functionTPDRichness,functionTPDc,
                                         functionredundancy,functiondissim,functiondensityProfileTPD,boot=99,
                                         TPDsAux,namesOrd,saveFile=paste0(getwd(),"/FDLoss.RDS")){
  
  dimensions <- data$dimensions
  traitsUSE <- data$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  cat(paste("\n MODEL FOR: ", targetGroupName, "\n"))
  alphaFunc<-list()
  
  spList<-spList[rownames(spList)%in%TPDsAux$data$species,]
  spOccurenceObsOrder <- matrix(0,nc=nrow(spList),nr=length(levels(spList$order)),dimnames = list(levels(spList$order),rownames(spList)))
  
  for(or in levels(spList$order)){
    spOccurenceObsOrder[or,colnames(spOccurenceObsOrder)%in%rownames(spList[spList$order==or,])]<-1
  }
  
  spOccurenceObsOrderLoss<-spOccurenceObsOrder
  spOccurenceObsOrderLoss[,focal]<-0
  TPDc_realms <- functionTPDc(TPDsAux,spOccurenceObsOrderLoss)
  alphaFunc[[1]]<- functionTPDRichness(TPDc_realms)
  
  t0 <- Sys.time()
  for (bo in 2:boot){
    sam<-numeric()
    for(or in names(table(spList[focal,]$order))){
      boOrd<-sample(rownames(spList[which(spList$order==or),]),length(which(spList[focal,]$order==or)))
      sam<-c(sam,boOrd)
    }
    spOccurenceBoot<-spOccurenceObsOrder
    spOccurenceBoot[,sam]<-0
    TPDc_realms <- functionTPDc(TPDsAux,spOccurenceBoot)
    alphaFunc[[bo]]<- functionTPDRichness(TPDc_realms)
    rm(TPDc_realms)
    t1 <- Sys.time()
    cat(paste0("\r GROUP:",targetGroupName," :",bo,"/",boot," .... Progress: ",round(bo/boot,4)*100,"% in ",hms_span(t0, t1),"\r"))
  }
  alphaFunc<-matrix(unlist(alphaFunc),nr=length(levels(spList$order)),ncol=boot,byrow = F,
                    dimnames = list(levels(spList$order),c("Obs",paste0("Boot ",1:(boot-1)))))
  
  alphaFunc<-cbind(Species=apply(spOccurenceObsOrder,1,sum),Sprm=table(spList[focal,]$order),alphaFunc)
  
  if(namesOrd=="all"){namesOrd<-levels(spList$order)}
  alphaFunc<-alphaFunc[namesOrd,]
  saveRDS(alphaFunc,saveFile) 
  return(alphaFunc)
}

make.redundancy<-function(data,spList,targetGroupName,functionTPDRichness,functionTPDc,
                          functionredundancy,functiondissim,functiondensityProfileTPD,boot=99,
                          TPDsAux,namesOrd,Rank,saveFile=paste0("Analyses/dataPrepared/AlphaDiv_RAW_",targetGroupName,"_Family.RDS")){
  
  dimensions <- data$dimensions
  traitsUSE <- data$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  cat(paste("\n MODEL FOR: ", targetGroupName, " for ",namesOrd,Rank,"\n"))
  alphaFunc<-list()
  
  spList<-spList[rownames(spList)%in%TPDsAux$data$species,]
  eval(parse(text=paste0("spOccurenceObsOrder <- matrix(0,nc=nrow(spList),nr=length(levels(spList$",Rank,")),
                         dimnames = list(levels(spList$",Rank,"),rownames(spList)))")))
  
  eval(parse(text=paste0("lvLRank<-levels(spList$",Rank,")")))
  
  for(or in lvLRank){
    eval(parse(text=paste0("spOccurenceObsOrder[or,colnames(spOccurenceObsOrder)%in%rownames(spList[spList$",Rank,"==or,])]<-1")))
  }
  
  TPDc_realms <- functionTPDc(TPDsAux,spOccurenceObsOrder)
  alphaFunc[[1]]<- functiondissim(TPDc_realms)
  
  t0 <- Sys.time()
  for (bo in 2:boot){
    spOccurenceBoot<-spOccurenceObsOrder
    for(or in 1:nrow(spOccurenceBoot)){
      spOccurenceBoot[or,]<-sample(spOccurenceBoot[or,])
    }
    TPDc_realms <- functionTPDc(TPDsAux,spOccurenceBoot)
    alphaFunc[[bo]]<- functiondissim(TPDc_realms)
    rm(TPDc_realms)
    t1 <- Sys.time()
    cat(paste0("\r GROUP:",targetGroupName," :",bo,"/",boot," .... Progress: ",round(bo/boot,4)*100,"% in ",hms_span(t0, t1),"\r"))
  }
  
  # if(namesOrd=="all"){namesOrd<-lvLRank}
  # if(!"all"%in%namesOrd){
  #   for (i in 1:length(alphaFunc)){
  #     alphaFunc[[i]][[1]]<-lapply(alphaFunc[[i]][[1]],function(x){x[namesOrd,namesOrd]})
  #   }
  # }
  saveRDS(alphaFunc,saveFile)

  return(alphaFunc)
}

make.redundancy.loss<-function(data,spList,focal,targetGroupName,functionTPDRichness,functionTPDc,
                               functionredundancy,functiondissim,functiondensityProfileTPD,boot=99,
                               TPDsAux,namesOrd,saveFile=paste0("Analyses/dataPrepared/AlphaDiv_RAW_",targetGroupName,"_Family.RDS")){
  set.seed(1)
  dimensions <- data$dimensions
  traitsUSE <- data$traitsUse
  colnames(traitsUSE)<-paste0("Comp.",1:dimensions)
  cat(paste("\n MODEL FOR: ", targetGroupName, "\n"))
  alphaFunc<-list()
  
  spList<-spList[rownames(spList)%in%TPDsAux$data$species,]
  spOccurenceObsOrder <- matrix(0,nc=nrow(spList),nr=length(levels(spList$order)),dimnames = list(levels(spList$order),rownames(spList)))
  
  for(or in levels(spList$order)){
    spOccurenceObsOrder[or,colnames(spOccurenceObsOrder)%in%rownames(spList[spList$order==or,])]<-1
  }
  
  spOccurenceObsOrderLoss<-spOccurenceObsOrder
  spOccurenceObsOrderLoss[,focal]<-0
  TPDc_realms <- functionTPDc(TPDsAux,spOccurenceObsOrderLoss)
  alphaFunc[[1]]<- functiondissim(TPDc_realms)
  
  t0 <- Sys.time()
  for (bo in 2:boot){
    
    sam<-numeric()
    for(or in names(table(spList[focal,]$order))){
      boOrd<-sample(rownames(spList[which(spList$order==or),]),length(which(spList[focal,]$order==or)))
      sam<-c(sam,boOrd)
    }
    spOccurenceBoot<-spOccurenceObsOrder
    spOccurenceBoot[,sam]<-0
    TPDc_realms <- functionTPDc(TPDsAux,spOccurenceBoot)
    alphaFunc[[bo]]<- functiondissim(TPDc_realms)
    rm(TPDc_realms)
    t1 <- Sys.time()
    cat(paste0("\r GROUP:",targetGroupName," :",bo,"/",boot," .... Progress: ",round(bo/boot,4)*100,"% in ",hms_span(t0, t1),"\r"))
  }
  
  if(!"all"%in%namesOrd){
    for (i in 1:length(alphaFunc)){
      alphaFunc[[i]][[1]]<-lapply(alphaFunc[[i]][[1]],function(x){x[namesOrd,namesOrd]})
    }
  }
  saveRDS(alphaFunc,saveFile) 
  return(alphaFunc)
}


make.functional.below<-function(species,traitsUSE,TPDsAux,PCAImpute,
                                level,Rank,taxo,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA,dim=1){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  colAbove <- "darkolivegreen"
  colBelow <- "darkorange4"
  
  if(dim == 1) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 2)]
  if(dim == 2) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 3)]
  if(dim == 3) loadingsUse <- PCAImpute$PCA$loadings[, c(1, 4)]
  if(dim == 4) loadingsUse <- PCAImpute$PCA$loadings[, c(2, 3)]
  if(dim == 5) loadingsUse <- PCAImpute$PCA$loadings[, c(2, 4)]
  if(dim == 6) loadingsUse <- PCAImpute$PCA$loadings[, c(3, 4)]
  
  
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,]
  var<-NULL
  var<-PCAImpute$variance
  if(is.null(var)){var<-PCAImpute$Variance}
  if(dim==1){var<-var[c(1,2)]}
  if(dim==6){var<-var[c(3,4)]}

  if(level == "All"){
    clade<-rownames(dataAnalysis)
  }  
  if(level == "Order"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$order==Rank,]$scientificName)),])
  }
  if(level == "Family"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$family==Rank,]$scientificName)),])
  }
  if(level == "Genus"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$genus==Rank,]$scientificName)),])
  }
  
  
  commNames <- c("All","Clade")
  commMatrix <- matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                       dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] <- 1
  commMatrix["Clade",clade ] <- 1
  
  commMatrix<-commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux <- TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab <- paste0("PC1 (", round(100 * var[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * var[2], 2), "%)")
  imageMat <-imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges <- unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges <- unique(TPDcAux$data$evaluation_grid[,2])
  
  traitEdge<-melt(imageMat[,,"All"])
  colnames(traitEdge) <- c('x','y','z')
  
  traitEdgeClade<-melt(imageMat[,,"Clade"])
  colnames(traitEdgeClade) <- c('x','y','z')
  
  xmin <- 0
  xmax <- 1
  ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
  ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
                              round( (max(imageMat[, , "All"], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
  
  
  dataPlot<-dataAnalysis[rownames(dataAnalysis)%in%TPDsAux$data$species,]
  dataPlot<-cbind.data.frame(dataPlot,sizeP=rep(1,nrow(dataPlot)),colorP=rep(as.character("grey38"),nrow(dataPlot)))
  dataPlot$colorP<-as.character(dataPlot$colorP)
  dataPlot$names<-rownames(dataPlot)

  colnames(dataPlot)[c(1,2)]<-c("x","y")
  dataPlot$z<-rep(1,nrow(dataPlot))
  multiplierArrow<-2.1
  
  if (level !="All"){
    dataPlot[clade,c("sizeP")]<-2
    dataPlot[clade,c("colorP")]<-"darkgreen"
    p <- ggplot(data = traitEdge, aes(x=x,y=y,z=z)) +
      geom_point(aes(text=rownames(dataPlot)),
                 data=dataPlot[,c("x","y","z")], colour = dataPlot$colorP, alpha = 1/2,size=dataPlot$sizeP)+
      geom_contour(breaks =c(.99),color="grey77") +
      #geom_contour(data=traitEdgeClade,breaks =c(.50,.99),color="darkgreen") +
      theme_bw()
  }else{
    
      dataPlot[gsub(" ","_",Rank),c("sizeP")]<-3
      dataPlot[gsub(" ","_",Rank),c("colorP")]<-"red"
    
    p <- ggplot(data = traitEdge, aes(x=x,y=y,z=z)) +
      geom_point(aes(text=rownames(dataPlot)),
                 data=dataPlot[,c("x","y","z")], colour = dataPlot$colorP, alpha = 1/2,size=dataPlot$sizeP)+
      geom_contour(breaks =c(.50,.99),color="grey47") +
      theme_bw()
  }
  ggplotly(p, tooltip="text")
}




make.functional.3d<-function(species,traitsUSE,TPDsAux,PCAImpute,
                                level,Rank,taxo,namesOrd,imageTPD,
                                ncolors,limX,limY,ColorRamp,title=NA,dim=1){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
 
  var<-NULL
  var<-PCAImpute$variance
  if(is.null(var)){var<-PCAImpute$Variance}
  if(dim==1){var<-var[c(1,2)]}
  if(dim==6){var<-var[c(3,4)]}
  
  colnames(traitsUSE)<-c("PC1","PC2","PC3","PC4")
  
  dataAnalysis<-traitsUSE[rownames(traitsUSE)%in%species,namesOrd]

   
  if(level == "All"){
    clade<-rownames(dataAnalysis)
  }  
  if(level == "Order"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$order==Rank,]$scientificName)),])
  }
  if(level == "Family"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$family==Rank,]$scientificName)),])
  }
  if(level == "Genus"){
    clade<-rownames(dataAnalysis[rownames(dataAnalysis)%in%gsub(" ","_",as.character(taxo[taxo$genus==Rank,]$scientificName)),])
  }
  
  xlab <- paste0("PC1 (", round(100 * var[1], 2), "%)") 
  ylab <- paste0("PC2 (", round(100 * var[2], 2), "%)")
  
  
  dataPlot<-dataAnalysis[rownames(dataAnalysis)%in%TPDsAux$data$species,]
  dataPlot<-cbind.data.frame(dataPlot,sizeP=rep(1,nrow(dataPlot)),colorP=rep(as.character("grey38"),nrow(dataPlot)))
  dataPlot$colorP<-as.character(dataPlot$colorP)
  dataPlot$names<-rownames(dataPlot)
  
  
  colnames(dataPlot)[c(1,2,3)]<-c("x","y","z")
  multiplierArrow<-2.1
  
  if (level !="All"){
    dataPlot[clade,c("sizeP")]<-2
    dataPlot[clade,c("colorP")]<-"darkgreen"
    
    
    plot_ly(dataPlot, x = ~x, y = ~y, z = ~z, color = ~colorP,colors = c("darkgreen","grey38"),
            hoverinfo = 'text', alpha = 0.75,
            text = rownames(dataPlot)) %>% layout(showlegend = FALSE,
                                                  scene = list(xaxis = list(title = namesOrd[1]),
                                                               yaxis = list(title = namesOrd[2]),
                                                               zaxis = list(title = namesOrd[3]))) 
  }else{
    
    dataPlot[gsub(" ","_",Rank),c("sizeP")]<-2
    dataPlot[gsub(" ","_",Rank),c("colorP")]<-"red"
    
    plot_ly(dataPlot, x = ~x, y = ~y, z = ~z, color = ~colorP,colors = c("grey38","red"),
            hoverinfo = 'text', alpha = 0.75,
            text = rownames(dataPlot)) %>% layout(showlegend = FALSE,
                                                  scene = list(xaxis = list(title = namesOrd[1]),
                                                               yaxis = list(title = namesOrd[2]),
                                                               zaxis = list(title = namesOrd[3]))) 
  }
  
  
}






 # end (00_PCA_functions.R) 


 # start (00_TD_functions.R)

 
descript_pieplot<-function(df, clade){
  
  
  if(clade == "order"){
    dfOrder<-data.frame(sort(table(df$order),decreasing = T))
    hc <-  billboarder() %>% bb_donutchart(dfOrder) %>% bb_legend(show = FALSE) %>%
      bb_donut(title=clade,label=list(expand=T))
  }

  
  if(clade == "family"){
    dffamily<-data.frame(sort(table(df$family),decreasing = T))
    dffamily<-dffamily[which(dffamily$Freq>20),]
    hc <- billboarder() %>% bb_donutchart(dffamily) %>% bb_legend(show = FALSE) %>%
      bb_donut(title=clade,label=list(expand=T)) 
  }

  
  if(clade == "genus"){
    dfgenus<-data.frame(sort(table(df$genus),decreasing = T))
    dfgenus<-dfgenus[which(dfgenus$Freq>30),]
    hc <- billboarder() %>% bb_donutchart(dfgenus) %>% bb_legend(show = FALSE) %>%
      bb_donut(title=clade,label=list(expand=T))
  }
hc
}


descript_barplot<-function(df,IUCN,clade){
  df2<-cbind(df[names(IUCN),],IUCN)
  df2$IUCN<-as.character(df2$IUCN)
  df2$IUCN<-factor(df2$IUCN)
  df2<-na.omit(df2)
  
  if(clade == "order"){
  hc <- df2 %>% 
    group_by(order)%>%
    count(IUCN)%>%
    mutate(freq = n/sum(n))%>%
    hchart(
      'column', hcaes(x = 'order',y=freq, group = 'IUCN'),
      stacking = "normal"
    ) %>%
    hc_yAxis(max = 1)
  }
  
  if(clade == "family"){
  hc <- df2 %>% 
    group_by(family)%>%
    count(IUCN)%>%
    mutate(freq = n/sum(n))%>%
    hchart(
      'column', hcaes(x = 'family',y=freq, group = 'IUCN'),
      stacking = "normal"
    ) %>%
    hc_yAxis(max = 1)
  }
  
  
  if(clade == "genus"){
  hc <- df2 %>% 
    group_by(genus)%>%
    count(IUCN)%>%
    mutate(freq = n/sum(n))%>%
    hchart(
      'column', hcaes(x = 'genus',y=freq, group = 'IUCN'),
      stacking = "normal"
    ) %>%
    hc_yAxis(max = 1)
  }
  
hc
}

 # end (00_TD_functions.R) 


 # start (00_Traits_functions.R)

 
hist_traits<-function(trait,taxo,nm,clade,Rank="Order"){
  trt <- na.omit(cbind(trait,taxo[rownames(trait),]))
  if(Rank=="Order"){
    sam <- rownames(trt[trt$order==clade,])
  }
  if(Rank=="Family"){
    sam <- rownames(trt[trt$family==clade,])
  }
  if(Rank=="Genus"){
    sam <- rownames(trt[trt$genus==clade,])
  }
  hc<-hchart(density(trt[,nm]),type="area", name = "all",color = "lightblue")%>% 
    hc_add_series(
      density(trt[sam,nm]), type = "area",
      color = "#B71C1C", 
      name = clade
    )
  hc

}

hist_traits_up<-function(trait,traitGlob,taxo,nm){
  trt <- na.omit(cbind(trait,taxo[rownames(trait),]))
  hchart(density(traitGlob[,nm]),type="area", name = "Global data base",color = "lightblue")%>%
  hc_add_series(
    density(trt[,nm]), type = "area",
    color = "#B71C1C",
    name = "Your data"
  )
}

box_traits<-function(trait,taxo,nm,clade,rank="Order"){
  if(rank=="Order"){
  trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
  sam <- rownames(trt[trt$order%in%clade,])
  trtsam<-cbind.data.frame(value=c(trt[sam,nm],trt[,nm]),clade=c(rep(clade,length(sam)),rep("All",nrow(trt))))
  hc<-hcboxplot(
    x = trtsam$value,
    var = trtsam$clade,
    name = nm, 
    color = "#2980b9"
  ) 
  }
  
  if(rank=="Family"){
    trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
    sam <- rownames(trt[trt$family%in%clade,])
    trtsam<-cbind.data.frame(value=c(trt[sam,nm],trt[,nm]),clade=c(rep(clade,length(sam)),rep("All",nrow(trt))))
    hc<-hcboxplot(
      x = trtsam$value,
      var = trtsam$clade,
      name = nm, 
      color = "#2980b9"
    ) 
  }
  
  if(rank=="Genus"){
    trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
    sam <- rownames(trt[trt$genus%in%clade,])
    trtsam<-cbind.data.frame(value=c(trt[sam,nm],trt[,nm]),clade=c(rep(clade,length(sam)),rep("All",nrow(trt))))
    hc<-hcboxplot(
      x = trtsam$value,
      var = trtsam$clade,
      name = nm, 
      color = "#2980b9"
    ) 
  }
  hc
}


make.stattabletrait.rea<-function(trait,taxo,nm,clade,rank="Order"){
  if(rank=="Order"){
    trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
    sam <- rownames(trt[trt$order%in%clade,])
  }
  
  if(rank=="Family"){
    trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
    sam <- rownames(trt[trt$family%in%clade,])
  }
  
  if(rank=="Genus"){
    trt <-  na.omit(cbind(trait,taxo[rownames(trait),]))
    sam <- rownames(trt[trt$genus%in%clade,])
  }
  
  test<-t.test(trt[,nm],trt[sam,nm])
  pval<-ifelse(test$p.value<0.05,"*",round(test$p.value,2))
  
  tab<-rbind(distrib(trt[,nm],3),distrib(trt[sam,nm],3))
  tab<-cbind(tab,stat=rbind(NA,pval))
  rownames(tab)<-c("All species",clade)
  colnames(tab)[8]<-"Test"
  tab
}

 # end (00_Traits_functions.R) 


 # start (Aux_Functions.R)

 
### densityProfile: function to estimate functional richness profiles (functional space occupied by probabilistic quantiles) from kde object
#### (x is a 'kde' object)
densityProfile <- function(x, probs=seq(0, 0.99, by=0.01)){
  TPD <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSize <- prod(cellEdge)
  TPD <- TPD/sum(TPD)
  alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
  FRicFunct <- function(TPD, alpha = 0.99) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM KDE OBJECT
    FRic <- numeric(length(alpha))
    for(i in 1:length(alpha)){
      TPDAux <- TPD
      greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) > alpha[i])])
      TPDAux[TPDAux < greater_prob] <- 0
      FRic[i] <- sum(TPDAux>0)* cellSize
    }
    names(FRic) <- alpha
    return(FRic)
  }
  result <- FRicFunct(TPD = TPD, alpha = probs)
  return(result)
}

### densityProfileTPD: function to estimate functional richness profiles (functional space occupied by probabilistic quantiles) from TPDc object
#### (x is a 'TPDc' object)
densityProfileTPD <- function(x, probs=seq(0, 1, by=0.01)){
  TPDList <- x$TPDc$TPDc
  results <- matrix(NA, nrow=length(TPDList), ncol=length(probs),
                    dimnames = list(names(TPDList), probs))
  for(comm in 1:length(TPDList)){
    TPD <- TPDList[[comm]]
    cellSize <- x$data$cell_volume
    alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
    FRicFunct <- function(TPD, alpha = 1) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM TPD OBJECT
      FRic <- numeric(length(alpha))
      for(i in 1:length(alpha)){
        TPDAux <- TPD
        greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) >= alpha[i])])
        TPDAux[TPDAux < greater_prob] <- 0
        FRic[i] <- sum(TPDAux>0)* cellSize
      }
      names(FRic) <- alpha
      return(FRic)
    }
    results[comm,] <- FRicFunct(TPD = TPD, alpha = probs)  
  }
  return(results)
}
### densityProfileTPD_large: function to estimate (functional space occupied by probabilistic quantiles) from TPDc_large object
#### x is a TPDc object created with the "large" version
densityProfileTPD_large <- function(x, probs=seq(0, 1, by=0.01)){
  TPDList <- x$TPDc$TPDc
  results <- matrix(NA, nrow=length(TPDList), ncol=length(probs),
                    dimnames = list(names(TPDList), probs))
  for(comm in 1:length(TPDList)){
    TPD <- TPDList[[comm]][,"notZeroProb"]
    cellSize <- x$data$cell_volume
    alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
    FRicFunct <- function(TPD, alpha = 1) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM TPD OBJECT
      FRic <- numeric(length(alpha))
      for(i in 1:length(alpha)){
        TPDAux <- TPD
        greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) >= alpha[i])])
        TPDAux[TPDAux < greater_prob] <- 0
        FRic[i] <- sum(TPDAux>0)* cellSize
      }
      names(FRic) <- alpha
      return(FRic)
    }
    results[comm,] <- FRicFunct(TPD = TPD, alpha = probs)  
  }
  return(results)
}

### overlapF: function to estimate probabilistic overlap between kde objects
overlapF <- function(x,y){
  if(!identical(x$eval.points, y$eval.points)) stop("Not identical eval points!")
  TPDx <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSizex <- prod(cellEdge)
  TPDx <- TPDx * cellSizex
  TPDx <- TPDx/sum(TPDx)
  
  TPDy<- y$estimate
  TPDy <- TPDy * cellSizex
  TPDy <- TPDy/sum(TPDy)
  OL <- sum(pmin(TPDx, TPDy))
  return(OL)
}

### percentileTPD: function to transform probabilities from TPDc object into quantiles
### x is a TPDc object 
percentileTPD <- function(x){
  TPDList <- x$TPDc$TPDc
  results <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, nrow(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"notZeroProb"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
    results[[comm]] <- TPDList[[comm]]
  }
  names(results) <- names(TPDList)
  return(results)
}

### imageTPD: fucntion to plot a TPDc object of two dimensions
# x is a TPDc object 
imageTPD <- function(x, thresholdPlot = 0.99){
  TPDList <- x$TPDc$TPDc
  imageTPD <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
                            prob = TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
    imageTPD[[comm]] <- TPDList[[comm]]
  }
  names(imageTPD) <- names(TPDList)
  spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
                             ncol = length(TPDList), 
                             dimnames = list(1:nrow(x$data$evaluation_grid),
                                             names(TPDList)))
  trait1Edges <- unique(x$data$evaluation_grid[,1])
  trait2Edges <- unique(x$data$evaluation_grid[,2])
  imageMat <- array(NA, c(length(trait1Edges), 
                          length(trait2Edges),
                          length(imageTPD)),
                    dimnames = list(trait1Edges, trait2Edges, names(TPDList)))
  for(comm in 1:length(TPDList)){
    percentileSpace <- x$data$evaluation_grid
    percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
    percentileSpace[, "percentile"] <-imageTPD[[comm]][,"percentile"]
    for(i in 1:length(trait2Edges)){
      colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
      imageMat[, i, comm] <- colAux$percentile 
    } 
    imageMat[imageMat > thresholdPlot] <- NA
  }
  return(imageMat)
}

### quantileTPD: function to transform probabilities from TPDc object into quantiles
### x is a TPDc object 
quantileTPD <- function(x, thresholdPlot = 0.99){
  TPDList <- x$TPDc$TPDc
  quantileTPD <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
                            prob = TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
    quantileTPD[[comm]] <- TPDList[[comm]]
  }
  names(quantileTPD) <- names(TPDList)
  spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
                             ncol = length(TPDList), 
                             dimnames = list(1:nrow(x$data$evaluation_grid),
                                             names(TPDList)))
  trait1Edges <- unique(x$data$evaluation_grid[,1])
  trait2Edges <- unique(x$data$evaluation_grid[,2])
  imageMat <- array(NA, c(length(trait1Edges), 
                          length(trait2Edges),
                          length(quantileTPD)),
                    dimnames = list(trait1Edges, trait2Edges, names(TPDList)))
  for(comm in 1:length(TPDList)){
    percentileSpace <- x$data$evaluation_grid
    percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
    percentileSpace[, "percentile"] <-quantileTPD[[comm]][,"percentile"]
    for(i in 1:length(trait2Edges)){
      colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
      imageMat[, i, comm] <- 1 - colAux$percentile 
    } 
    imageMat[imageMat < (1 - thresholdPlot)] <- 0
  }
  return(imageMat)
}

### MAE_TPD: function to estimate differences of quantiles between two quantileTPD objects
MAE_TPD <- function(x, y){# x and y are made with quantileTPD function
  diffXY <- abs(x - y)
  MAE <- mean(diffXY)
  return(MAE)
}  

### quantileTPD_large: function to transform probabilities from TPDc_large object into quantiles
### x is a TPDc_large object 
quantileTPD_large <- function(x){
  TPDList <- x$TPDc$TPDc
  results <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, nrow(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- 1- cumsum(TPDList[[comm]][,"notZeroProb"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
    results[[comm]] <- TPDList[[comm]]
  }
  names(results) <- names(TPDList)
  return(results)
}

### MAE_TPD_large: function to estimate differences of quantiles between two quantileTPD_large objects
MAE_TPD_large <- function(x, y){# x and y are made with quantileTPD_large function
  mergeXY <- merge(x, y, by="notZeroIndex", all=T)
  mergeXY[is.na(mergeXY)] <- 0
  mergeXY$change <- abs(mergeXY$percentile.x - mergeXY$percentile.y)
  MAE <- mean(mergeXY$change)
  return(MAE)
}

### TPDsMean_large: function to estimate TPDs functions using a single average value per species and a given bandwidth (standard deviation). 
# This function is equivalent to TPD::TPDsMean, but does not record the cells with zero probability, making it more suitable for higher dimensions.
TPDsMean_large<- function(species, means, sds, alpha = 0.95, samples = NULL,
                          trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {
  
  # INITIAL CHECKS:
  #   1. Compute the number of dimensions (traits):
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
         number of dimensions")
  }
  #   2. sds and means must have the same dimensions:
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  #   3. species and means must have the same "dimensions":
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means'
         and 'sds'")
  }
  #	4. NA's not allowed in means, sds & species:
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  #	5. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means'
         and 'sds'")
    }
    if (any(is.na(samples))) {
      stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  }
  
  #	6. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	6. Create the grid of cells in which the density function is evaluated:
  if (is.null(n_divisions)) {
    n_divisions_choose<- c(1000, 200, 50, 25)
    n_divisions<- n_divisions_choose[dimensions]
  }
  grid_evaluate<-list()
  edge_length <- list()
  cell_volume<-1
  for (dimens in 1:dimensions){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] -
      grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions==1){
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  # Creation of lists to store results:
  results <- list()
  # DATA: To store data and common features
  results$data <- list()
  results$data$evaluation_grid <- evaluation_grid
  results$data$cell_volume <- cell_volume
  results$data$edge_length <- edge_length
  results$data$species <- species
  results$data$means <- means
  results$data$sds <- sds
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }
  
  results$data$alpha <- alpha
  results$data$pop_means <- list()
  results$data$pop_sds <- list()
  results$data$pop_sigma <- list()
  results$data$dimensions <- dimensions
  results$data$type <- type
  results$data$method <- "mean"
  
  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()
  
  
  ########Multivariate normal density calculation
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi]<-
      unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      
      multNormAux <- mvtnorm::dmvnorm(x = evaluation_grid,
                                      mean = results$data$pop_means[[spi]],
                                      sigma = results$data$pop_sigma[[spi]])
      multNormAux <- multNormAux / sum(multNormAux)
      # Now, we extract the selected fraction of volume (alpha), if necessary
      extract_alpha <- function(x){
        # 1. Order the 'cells' according to their probabilities:
        alphaSpace_aux <- x[order(x, decreasing=T)]
        # 2. Find where does the accumulated sum of ordered probabilities becomes
        #   greater than the threshold (alpha):
        greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux ) > alpha) [1]]
        # 3. Substitute smaller probabilities with 0:
        x[x < greater_prob] <- 0
        # 5. Finally, reescale, so that the total density is finally 1:
        x <- x / sum(x)
        return(x)
      }
      if (alpha < 1){
        multNormAux <- extract_alpha(multNormAux)
      }
      notZeroIndex <- which(multNormAux != 0)
      notZeroProb <- multNormAux[notZeroIndex]
      results$TPDs[[spi]] <- cbind(notZeroIndex, notZeroProb)
      
      
    }
    if (dimensions == 1) stop("This function is intended for > 1 dimension")
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}

### TPDc_large: function to estimate Trait Probability Density of Communities based on TPDs_large objects
TPDc_large <- function(TPDs, sampUnit){
  sampUnit <- as.matrix(sampUnit)
  # 1. species names:
  if (is.null(colnames(sampUnit)) | any(is.na(colnames(sampUnit)))) {
    stop("colnames(sampUnit), must contain the names of the species.
      NA values are not allowed")
  }
  # 2. communities names:
  if (is.null(rownames(sampUnit)) | any(is.na(rownames(sampUnit)))){
    stop("rownames(sampUnit), must contain the names of the sampling units.
      NA values are not allowed")
  }
  # 3. Data values will be inherithed from TPDs, which must be of class TPD
  if (class(TPDs) != "TPDsp") {
    stop("TPDs must be an object of class 'TPDsp', created with the function
      'TPDs'")
  }
  species <- samples <- abundances <- numeric()
  for (i in 1:nrow(sampUnit)){
    samples <- c(samples, rep(rownames(sampUnit)[i], ncol(sampUnit)))
    species <- c(species, colnames(sampUnit))
    abundances <- c(abundances, sampUnit[i, ])
  }
  nonZero <- which(abundances > 0)
  samples <- samples[nonZero]
  species <- species[nonZero]
  abundances <- abundances[nonZero]
  
  # Creation of lists to store results:
  results <- list()
  results$data <- TPDs$data
  results$data$sampUnit <- sampUnit
  type <- results$data$type
  # All the species or populations in 'species' must be in the species or
  #   populations of TPDs:
  if (type == "Multiple populations_One species" |
      type == "Multiple populations_Multiple species") {
    species_base <- paste(species, samples, sep = ".")
    if (!all(unique(species_base) %in% unique(results$data$populations))) {
      non_found_pops <- which(unique(species_base) %in%
                                unique(results$data$populations) == 0)
      stop("All the populations TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_pops], collapse=" / "))
    }
  }
  if (type == "One population_One species" |
      type == "One population_Multiple species") {
    species_base <- species
    if (!all(unique(species_base) %in% unique(results$data$species))) {
      non_found_sps <- which(unique(species_base) %in%
                               unique(results$data$species) == 0)
      stop("All the species TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_sps], collapse=" / "))
    }
  }
  # END OF INITIAL CHECKS
  # TPDc computation
  results$TPDc <- list()
  results$TPDc$species <- list()
  results$TPDc$abundances <- list()
  results$TPDc$speciesPerCell <- list()
  # results$TPDc$RTPDs <- list()
  results$TPDc$TPDc <- list()
  
  for (samp in 1:length(unique(samples))) {
    selected_rows <- which(samples == unique(samples)[samp])
    species_aux <- species_base[selected_rows]
    abundances_aux <- abundances[selected_rows] / sum(abundances[selected_rows])
    RTPDsAux <- rep(0, nrow(results$data$evaluation_grid))
    TPDs_aux <- TPDs$TPDs[names(TPDs$TPDs) %in% species_aux]
    cellsOcc <- numeric()
    for (sp in 1:length(TPDs_aux)) {
      selected_name <- which(names(TPDs_aux) == species_aux[sp])
      cellsToFill <- TPDs_aux[[selected_name]][,"notZeroIndex"]
      cellsOcc <- c(cellsOcc, cellsToFill)
      probsToFill <- TPDs_aux[[selected_name]][,"notZeroProb"] * abundances_aux[sp]
      RTPDsAux[cellsToFill] <- RTPDsAux[cellsToFill] + probsToFill
    }
    TPDc_aux <- RTPDsAux
    notZeroIndex <- which(TPDc_aux != 0)
    notZeroProb <- TPDc_aux[notZeroIndex]
    results$TPDc$TPDc[[samp]] <- cbind(notZeroIndex, notZeroProb)
    results$TPDc$species[[samp]] <- species_aux
    results$TPDc$abundances[[samp]] <- abundances_aux
    results$TPDc$speciesPerCell[[samp]] <- table(cellsOcc)
    names(results$TPDc$TPDc)[samp] <-
      names(results$TPDc$species)[samp] <- names(results$TPDc$abundances)[samp] <-
      names(results$TPDc$speciesPerCell)[samp] <- unique(samples)[samp]
  }
  class(results) <- "TPDcomm"
  return(results)
}

### TPDRichness: Function to estimate functional richness from TPD object
TPDRichness <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      # message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    } else {
      results$populations <- list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### TPDRichness_large: Function to estimate functional richness from TPD_large object
TPDRichness_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      # message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    } else {
      results$populations <- list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### REND_large: Functional Evenness, Richness and Divergence from TPD_large object
REND_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # 2. Functional Evenness
  Calc_FEve <- function(x) {
    results_FE <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_eve <- rep((1 / length(TPD_aux)), times = length(TPD_aux))
      results_FE[i] <- sum(pmin(TPD_aux, TPD_eve))
    }
    names(results_FE) <- names_aux
    return(results_FE)
  }
  # 3. Functional Divergence
  Calc_FDiv <- function(x) {
    results_FD <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) { 
      notZeroCells <- TPD[[i]][,"notZeroIndex"]
      functional_volume <- evaluation_grid[notZeroCells, , drop=F]
      # Functional volume has to be standardised so that distances are
      # independent of the scale of the axes
      for (j in 1:ncol(functional_volume)){
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) /
          (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      # 1. Calculate the center of gravity
      COG <- colMeans(functional_volume, na.rm=T)
      # 2. Calculate the distance of each point in the functional volume to the
      #   COG:
      dist_COG <- function(x, COG) {
        result_aux<-stats::dist(rbind(x, COG))
        return(result_aux)
      }
      COGDist <- apply(functional_volume, 1, dist_COG, COG)
      # 3. Calculate the mean of the COGDist's
      meanCOGDist <- mean(COGDist)
      # 4. Calculate the sum of the abundance-weighted deviances for distaces
      #   from the COG (AWdistDeviances) and the absolute abundance-weighted
      #   deviances:
      distDeviances <- COGDist - meanCOGDist
      AWdistDeviances <- sum(TPD_aux * distDeviances)
      absdistDeviances <- abs(COGDist - meanCOGDist)
      AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
      #Finally, calculate FDiv:
      results_FD[i] <- (AWdistDeviances + meanCOGDist) /
        ( AWabsdistDeviances +  meanCOGDist)
    }
    names(results_FD) <- names_aux
    return(results_FD)
  }
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
    message("Calculating FEvenness of communities")
    results$communities$FEvenness <- Calc_FEve(TPDc)
    message("Calculating FDivergence of communities")
    results$communities$FDivergence <- Calc_FDiv(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of species")
      results$species$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of species")
      results$species$FDivergence <- Calc_FDiv(TPDs)
    } else {
      results$populations <- list()
      message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of populations")
      results$populations$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of populations")
      results$populations$FDivergence <- Calc_FDiv(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### redundancy_large:  Functional Redundancy of Communities from TPD_large object
redundancy_large <- function(TPDc = NULL) {
  if (class(TPDc) != "TPDcomm") {
    stop("TPDc must be an object of class TPDcomm generated with the TPDc
		    function")
  }
  x <- TPDc
  results <- list()
  results$redundancy <- results$richness <- numeric()
  for (i in 1:length(x$TPDc$TPDc)) {
    TPDc_aux <- x$TPDc$TPDc[[i]][, "notZeroProb"]
    M <- x$TPDc$speciesPerCell[[i]]
    results$redundancy[i] <- sum(M * TPDc_aux) - 1
    results$richness[i] <- sum(x$TPDc$abundances[[i]] >0)
  }
  results$redundancyRelative <- results$redundancy / (results$richness -1)
  names(results$redundancy) <- names(results$richness) <-
    names(results$redundancyRelative) <- names(x$TPDc$TPDc)
  return(results)
}

dissim_large <- function(x = NULL) {
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (class(x) == "TPDcomm"){
    TPDType<-"Communities"
    TPDc <- x
  } else{
    if (class(x) == "TPDsp"){
      TPDType<-"Populations"
      TPDs <- x
    } else{
      stop("x must be an object of class TPDcomm or TPDsp")
    }
  }
  results <- list()
  Calc_dissim <- function(x) {
    # 1. BetaO (functional dissimilarity)
    results_samp <- list()
    if (TPDType == "Communities") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
    }
    if (TPDType == "Populations") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
    }
    results_samp$dissimilarity <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                         dimnames = list(names_aux, names_aux))
    results_samp$P_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                    dimnames = list(names_aux, names_aux))
    results_samp$P_non_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                        dimnames = list(names_aux, names_aux))
    for (i in 1:length(TPD)) {
      TPD_i <- TPD[[i]]
      for (j in 1:length(TPD)) {
        if (i > j) {
          TPD_j <- TPD[[j]]
          commonTPD <- rbind(TPD_i, TPD_j)
          #Now, select only those cells that appear in both i and j
          duplicatedCells <- names(which(table(commonTPD[,"notZeroIndex"])==2))
          doubleTPD <- commonTPD[which(commonTPD[,"notZeroIndex"] %in% duplicatedCells), ]
          O_aux <- sum(tapply(doubleTPD[,"notZeroProb"], doubleTPD[,"notZeroIndex"], min))
          A_aux <- sum(tapply(doubleTPD[,"notZeroProb"], doubleTPD[,"notZeroIndex"], max)) - O_aux
          only_in_i_aux <- which(TPD_i[,"notZeroIndex"] %in% setdiff(TPD_i[,"notZeroIndex"], 
                                                                     TPD_j[,"notZeroIndex"]))	
          B_aux <- sum(TPD_i[only_in_i_aux,"notZeroProb"])
          only_in_j_aux <- which(TPD_j[,"notZeroIndex"] %in% setdiff(TPD_j[,"notZeroIndex"], 
                                                                     TPD_i[,"notZeroIndex"]))	
          C_aux <- sum(TPD_j[only_in_j_aux,"notZeroProb"])
          results_samp$dissimilarity[i, j] <-
            results_samp$dissimilarity[j, i] <-	1 - O_aux
          if (results_samp$dissimilarity[j, i] == 0) {
            results_samp$P_non_shared[i, j] <- NA
            results_samp$P_non_shared[j, i] <- NA
            results_samp$P_shared[i, j] <- NA
            results_samp$P_shared[j, i] <- NA
          }	else {
            results_samp$P_non_shared[i, j] <-
              results_samp$P_non_shared[j, i] <-
              (2 * min(B_aux, C_aux)) / (A_aux + 2 * min(B_aux, C_aux))
            results_samp$P_shared[i, j] <- results_samp$P_shared[j, i] <-
              1 - results_samp$P_non_shared[i, j]
          }
        }
        if (i == j) {
          results_samp$dissimilarity[i, j] <- 0
        }
      }
    }
    return(results_samp)
  }
  # IMPLEMENTATION
  if (TPDType == "Communities") {
    message("Calculating dissimilarities between ", length(TPDc$TPDc$TPDc)," communities. It might take some time")
    results$communities <- Calc_dissim(TPDc)
  }
  if (TPDType == "Populations") {
    message("Calculating dissimilarities between ", length(TPDs$TPDs)," populations. It might take some time")
    results$populations <- Calc_dissim(TPDs)
  }
  class(results) <- "OverlapDiss"
  return(results)
}










# 
# ####################
# #################### TRANSFORMING PROBABILITIUES INTO QUANTILES
# 
# 
# # x is a TPDc object created with the "large" version
# PlotPercentileTPD_large <- function(x, whichPlot = 1, ncolors=1000, colorsPick="YlOrBr", 
#                               xlab="Trait 1", ylab="Trait2",
#                               cont=c(50,90), xlim=NULL, ylim=NULL){
#   TPDList <- x$TPDc$TPDc
#   percentilesTPD <- list()
#   for(comm in 1:length(TPDList)){
#     percentile <- rep(NA, nrow(TPDList[[comm]]))
#     TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
#     orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
#     TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
#     TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"notZeroProb"]) *100
#     TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
#     percentilesTPD[[comm]] <- TPDList[[comm]]
#   }
#   names(percentilesTPD) <- names(TPDList)
#   spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
#                              ncol = length(TPDList), 
#                              dimnames = list(1:nrow(x$data$evaluation_grid),
#                                              names(TPDList)))
#   trait1Edges <- unique(x$data$evaluation_grid[,1])
#   trait2Edges <- unique(x$data$evaluation_grid[,2])
#   imageMat <- array(NA, c(length(trait1Edges), 
#                           length(trait2Edges),
#                           length(percentilesTPD)))
#   ncolors <- 1000
#   ColorRamp <- colorRampPalette(RColorBrewer::brewer.pal(9,colorsPick)[c(9:7, 5:3)])(ncolors)
#   
#   for(comm in whichPlot){
#     percentileSpace <- x$data$evaluation_grid
#     percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
#     percentileSpace[percentilesTPD[[comm]][,"notZeroIndex"], "percentile"] <- 
#       percentilesTPD[[comm]][,"percentile"]
#     for(i in 1:length(trait2Edges)){
#       colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
#       imageMat[, i, comm] <- colAux$percentile 
#     } 
#     xmin <- 0
#     xmax <- 100
#     ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
#     ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
#                                 round( (max(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
#     if(is.null(xlim)) xlim <- range(trait1Edges)
#     if(is.null(ylim)) ylim <- range(trait2Edges)
#     image(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], col=ColorRamp_ex,
#           xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=xlim, ylim=ylim,
#           main = names(TPDList)[comm])
#     box(which="plot")
#     axis(1,tcl=0.3,lwd=0.8)
#     axis(2, las=1, tcl=0.3,lwd=0.8)
#     contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], levels=cont,
#             drawlabels = T, labcex = 0.8, lwd=0.5, lty=1, col="black", add=T)
#   }
# }
# 
# # x is a TPDc object created with the "large" version
# PlotPercentileTPD <- function(x, whichPlot = 1, ncolors=1000, colorsPick="YlOrBr", 
#                                     xlab="Trait 1", ylab="Trait2",
#                                     cont=c(0.5,0.9), xlim=NULL, ylim=NULL,
#                               thresholdPlot = 0.99, gradientColorsF=NULL){
#   TPDList <- x$TPDc$TPDc
#   percentilesTPD <- list()
#   for(comm in whichPlot){
#     percentile <- rep(NA, length(TPDList[[comm]]))
#     TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
#                             prob = TPDList[[comm]], percentile)
#     orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
#     TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
#     TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
#     TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
#     percentilesTPD[[comm]] <- TPDList[[comm]]
#   }
#   names(percentilesTPD) <- names(TPDList)[whichPlot]
#   spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
#                              ncol = length(whichPlot), 
#                              dimnames = list(1:nrow(x$data$evaluation_grid),
#                                              names(TPDList)[whichPlot]))
#   trait1Edges <- unique(x$data$evaluation_grid[,1])
#   trait2Edges <- unique(x$data$evaluation_grid[,2])
#   imageMat <- array(NA, c(length(trait1Edges), 
#                           length(trait2Edges),
#                           length(percentilesTPD)))
#   ncolors <- 1000
#   ColorRamp <- rev(gradientColorsF(ncolors))
#   
#   for(comm in whichPlot){
#     percentileSpace <- x$data$evaluation_grid
#     percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
#     percentileSpace[, "percentile"] <-percentilesTPD[[comm]][,"percentile"]
#     for(i in 1:length(trait2Edges)){
#       colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
#       imageMat[, i, comm] <- colAux$percentile 
#     } 
#     imageMat[imageMat >thresholdPlot] <- NA
#     xmin <- 0
#     xmax <- 1
#     ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
#     ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
#                                 round( (max(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
#     if(is.null(xlim)) xlim <- range(trait1Edges)
#     if(is.null(ylim)) ylim <- range(trait2Edges)
#     image(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], col=ColorRamp_ex,
#           xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=xlim, ylim=ylim,
#           main = "", asp = 1)
#     box(which="plot")
#     axis(1,tcl=0.3,lwd=0.8)
#     axis(2, las=1, tcl=0.3,lwd=0.8)
#     contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], levels=cont,
#             drawlabels = T, labcex = 0.8, lwd=0.5, lty=1, col="black", add=T)
#   }
# }
# 
# 
# 
# ###################
# ###################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
distrib<-function(x,k=6) {
  # NA removing
  x<-x[is.na(x)==F]
  
  Nx<-length(x)
  sum<-round(summary(x),k)
  res<-c(Nx,sum) ; names(res)<-c("N",names(sum))
  return(res)
} #end of function distrib


hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

sesandpvalues_bis <- function (obs,rand,nreps,probs=c(0.025,0.975),rnd = 2){
  SES <- (obs - mean(rand)) / sd(rand)
  pValsSES <- rank(c(obs,rand))[1] / (length(rand) + 1)  
  results <- round(c(obs,SES, mean(rand),quantile(rand,prob=probs),pValsSES,nreps), rnd)
  names(results)<- c("Observed","SES","MeanRd","CI025Rd","CI975Rd","Pval","Nreps")
  return(results)
}

colGradient <- c("white",  "yellow", "red")
gradientColorsF <- colorRampPalette(colGradient, space = "Lab")
ncolors <- 1000
ColorRamp <- rev(gradientColorsF(ncolors))
contourLevels <- c(0.5, 0.99)


 # end (Aux_Functions.R) 


 # start (Aux_PlotVariables.R)

 
library(colorspace)
###### SET GRAPHIC OPTIONS

niceplot<-function(limX=c(0,1),limY=c(0,1), marg=c(4.5,4.5,2,2),
                   tick=-0.4,   labX=c(),labY=c(),    nmlabX=c(),nmlabY=c(),
                   lasX=1,lasY=1,    lineX=-0.1, lineY=0.1,   cexX=0.9, cexY=0.9,
                   nmX="X",nmY="Y",   lineXt=lineX+2, lineYt=lineY+2.5,   cexXt=1, cexYt=1,main=c(),cex=.5  ) {
  
  par(mar=marg)   # margins
  plot(limX,limY,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=limX,ylim=limY,main=main) # graphic window
  rect(limX[1],limY[1],limX[2],limY[2])   # border
  if (is.na(sum(labX))==F) {
    labx<-base::pretty(limX,n=7,min.n=5) ; labx<-labx[which(labx>=min(limX) & labx<=max(limX))] ; nmlabx<-labx  # default
    if (length(labX)>0) { labx<-labX ; nmlabx<-nmlabX }                                                   # customized
    axis(side=1, at=labx, labels=F, tcl=tick, pos=limY[1])  # thicks
    mtext(side=1, nmlabx, at=labx, line=lineX, cex=cexX, las=lasX) # labels
    mtext(side=1,nmX,cex=cexXt,line=lineXt,font=2) # title   
  } # end of if labels
  if (is.na(sum(labY))==F) {
    laby<-base::pretty(limY,n=7,min.n=5) ; laby<-laby[which(laby>=min(limY) & laby<=max(limY))] ; nmlaby<-laby  # default
    if (length(labY)>0) { laby<-labY ; nmlaby<-nmlabY }                                                   # customized
    axis(side=2, at=laby, labels=F, tcl=tick, pos=limX[1]) # thicks 
    mtext(side=2, nmlaby, at=laby, line=lineY, cex=cexY, las=lasY) # labels
    mtext(side=2,nmY,cex=cexYt,line=lineYt,font=2) # title  
  } # end of if labels
  
}
### COLORS
### Transparency
alphaCol<-60
# COLORS FOR THREAT, NO THREAT, WORLD CONTOUR AND DENSITY GRADIENT
# colNoThreatAux <- rgb(red=0, green=84, blue=150, alpha=alphaCol, maxColorValue = 255)
# colThreatAux <- rgb(red=202, green=108, blue=24, alpha=alphaCol, maxColorValue = 255)
# colThreat <- c(colNoThreatAux, colThreatAux)
colorsPick <- "RdYlGn"
colThreat <- brewer.pal(11,colorsPick)[c(9,3)]
colThreatLines <- brewer.pal(11,colorsPick)[c(11,1)]
colWorld <- "black" #rgb(red=123, green=10, blue=107, alpha=255, maxColorValue = 255)
colGradient <- c("white", "white", "yellow", "red")

for (i in 1:length(colThreat)){
  colThreat[i] <- rgb(t(col2rgb(colThreat)/255), alpha=alphaCol/255)[i]
}
### CONTOUR LEVELS:
conLines<-c(50, 99)
lineTypeWorld <- 1
lineTypeContour <- 1
lineTypeWorldExt <- 1
thickCountour <- 2+(length(conLines):1)
thickWorldExt <- max(thickCountour)

ptSize<-30
cexMain <- 1.25

#### x is a kde object
densityProfile <- function(x, probs=seq(0, 0.99, by=0.01)){
  TPD <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSize <- prod(cellEdge)
  TPD <- TPD/sum(TPD)
  alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
  FRicFunct <- function(TPD, alpha = 0.99) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM KDE OBJECT
    FRic <- numeric(length(alpha))
    for(i in 1:length(alpha)){
      TPDAux <- TPD
      greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) > alpha[i])])
      TPDAux[TPDAux < greater_prob] <- 0
      FRic[i] <- sum(TPDAux>0)* cellSize
    }
    names(FRic) <- alpha
    return(FRic)
  }
  result <- FRicFunct(TPD = TPD, alpha = probs)
  return(result)
}
overlapF <- function(x,y){
  if(!identical(x$eval.points, y$eval.points)) stop("Not identical eval points!")
  TPDx <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSizex <- prod(cellEdge)
  TPDx <- TPDx * cellSizex
  TPDx <- TPDx/sum(TPDx)
  
  TPDy<- y$estimate
  TPDy <- TPDy * cellSizex
  TPDy <- TPDy/sum(TPDy)
  OL <- sum(pmin(TPDx, TPDy))
  return(OL)
}


### PCA Arrows:
multArrow <- c(0.8, 1.5, 2.2, 1, 0.8, 1)
multiplierTextGr <- c(1.2, 1.2, 1.25, 1.15, 1.3, 1)
### legend position:
legPos <- c("bottomleft", "bottomright", "topright", "bottomleft", "bottomleft", "topleft")
linesLegPos <- c("bottomright", "bottomleft", "topleft", "bottomright", "bottomright", "topright")
# 
# 
# limXlist <- list(list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 8.5), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 9), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 8), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 8.5), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 8), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 5.5), #Plants
#                        c(-4.6, 8), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7)),
#                  list( c(-5.7, 6.5), #Plants
#                        c(-4.6, 8.5), #Mammals
#                        c(-4.7, 10.1), #Aves
#                        c(-5.8, 7), #reptiles
#                        c(-5, 5.9),#Amphibians
#                        c(-7, 7))) #FWFish
# #
# limYlist <- list(list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)),
#                  list( c(-8.5, 6.5), #Plants
#                        c(-3.5, 4), #Mammals
#                        c(-3, 6.1), #Aves
#                        c(-6.5, 5.5), #reptiles
#                        c(-5.5,5),
#                        c(-7, 7)))

 # end (Aux_PlotVariables.R) 

