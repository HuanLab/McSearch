#########################################
# This script is for MCS library construction.
# Shipei Xing, Mar 12, 2020
# Copyright @ The University of British Columbia

#########################################
# Parameter setting part
database.path <- 'E:/McSearch/library'    # Data path for the MS/MS spectral database. Use "/" instead of "\".
db.name <- 'test.msp'                     # The MS/MS spectral database name ending with ".msp"
mz.tol <- 0.01                            # The m/z tolerance used for MCS library construction (default 0.01 Da)
topNo <- 30                               # The number of highest peaks in each spectrum will be used, default 30.
HNL.threshold <- 36                       # The mass threshold for HNL values reserved in MCS library, default 36 Da.
exclusion.after.premass <- FALSE          # Logical: if TRUE, all the fragments larger than the precursor mass will be excluded. FALSE by default.

#########################################

# Main program
setwd(database.path)
library("CHNOSZ")
library('metaMS')
db <- read.msp(db.name, only.org = FALSE,org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

element.table <- data.frame(matrix(ncol=3,nrow = 10))
colnames(element.table) <- c('element','mass','number')
element.table[,1] <- c('C','H','N','O','P','S','F','Cl','Br','I')
element.table[,2] <- c(12,1.007825,14.003074,15.994915,30.973762,31.972071,18.998403,34.968853,78.918337,126.904473)
element.table[,3] <- rep(0,10)
formulafilter <- function(std.formula){
  table <- count.elements(std.formula)
  a <- element.table
  for(m in 1:length(table)){a[a[,1]==names(table)[m],3] <- table[m]}
  validHNL <- vector(mode = 'numeric')
  k <- 1
  for(c in 3:a[1,3]){for(h in 1:min(3*c,a[2,3])){for(n in 0:min(1.3*c,a[3,3])){for(o in 0:min(1.2*c,a[4,3])){for(p in 0:a[5,3]){for(s in 0:a[6,3])
    {for(f in 0:a[7,3]){for(Cl in 0:a[8,3]){for(Br in 0:a[9,3]){for(I in 0:a[10,3]){
    validHNL[k] <- c*a[1,2]+h*a[2,2]+n*a[3,2]+o*a[4,2]+p*a[5,2]+s*a[6,2]+f*a[7,2]+Cl*a[8,2]+Br*a[9,2]+I*a[10,2]
    k <- k +1
  }}}}}}}}}}
  return(validHNL)
}

new.db <- db
badNo.vector <- vector(mode = 'numeric')
z<-1
for(i in 1:length(db)){
  ms2 <- as.data.frame(db[[i]]$pspectrum)
  if(nrow(ms2)==1){
    if(decimalplaces(ms2[1,1])<=2){
      badNo.vector[z] <- i
      z<-z+1
    }
  }
  if(nrow(ms2)>1){
    if(decimalplaces(ms2[1,1])<=2 & decimalplaces(ms2[2,1])<=2){
      badNo.vector[z] <- i
      z<-z+1
    }
  }
  if(is.element(i,badNo.vector)==FALSE){
    if(is.null(db[[i]]$Formula)){
      badNo.vector[z] <- i
      z<-z+1
    }
  }
  if(is.null(db[[i]]$PrecursorMZ)==FALSE){premass <- round(db[[i]]$PrecursorMZ,4)}
  if(is.null(db[[i]]$Formula)==FALSE){
    formula <- db[[i]]$Formula
    if(grepl('\\[',formula)){
      formula <- substring(formula,regexpr("\\[",formula)+1,regexpr("\\]",formula)-1)
    }
    else if(grepl('\\+',formula) | grepl('\\-',formula)){
      formula <- substring(formula,1,nchar(formula)-1)
    }
  }

  ms2[,2] <- 100*ms2[,2]/max(ms2[,2])
 
  # exclude after premass
  if(is.null(db[[i]]$PrecursorMZ)==FALSE & exclusion.after.premass==TRUE){
    ms2 <- ms2[ms2[,1]<=(premass+mz.tol),]
  }
  
  # top 'topNo' peaks
  if(nrow(ms2)>topNo){
    ms2 <- ms2[ms2[,2] > sort(ms2[,2],decreasing=TRUE)[topNo+1],]
  }
  
  # square root transformation
  ms2[,2] <- sqrt(ms2[,2])

  # positive mode: (1,0)
  # negative mode: (-1,0)
  if(is.null(db[[i]]$Ion_mode)==FALSE){
    if(db[[i]]$Ion_mode %in% c('P','p')){ms2 <- rbind(c(1.007276,0),ms2)}
    if(db[[i]]$Ion_mode %in% c('n','N')){ms2 <- rbind(c(-1.007276,0),ms2)}
  }
  if(is.null(db[[i]]$Ion_mode)==TRUE){
    if(is.null(db[[i]]$Precursor_type)==FALSE){
      str <- substr(db[[i]]$Precursor_type,nchar(db[[i]]$Precursor_type),nchar(db[[i]]$Precursor_type))
      if(str=="+"){ms2 <- rbind(c(1.007276,0),ms2)}
      if(str=="-"){ms2 <- rbind(c(-1.007276,0),ms2)}
      if(str!='+' & str!='-'){
        if(is.element(i,badNo.vector)==FALSE){
          badNo.vector[z] <- i
          z<-z+1
        }
      }
    }
    if(is.null(db[[i]]$Precursor_type)){
      if(is.element(i,badNo.vector)==FALSE){
        badNo.vector[z] <- i
        z<-z+1
      }
    }
  }
  
  # if premass is not in the MS2 spectrum, add (premass, 0)
  if(is.null(db[[i]]$PrecursorMZ)==FALSE){
    if(decimalplaces(db[[i]]$PrecursorMZ)>=3 & min(abs(ms2[,1]-premass)) > mz.tol){ms2 <- rbind(ms2,c(premass,0))}
  }
  
  HNL <- data.frame(matrix(ncol=6))
  colnames(HNL) <- c('HNL','mz.a','mz.b','int.a','int.b','HNL.int')
  h <- 1
  for(m in 1:(nrow(ms2)-1)){
    for(n in (m+1):nrow(ms2)){
      HNL[h,1] <- ms2[n,1] - ms2[m,1]
      HNL[h,2] <- ms2[n,1]
      HNL[h,3] <- ms2[m,1]
      HNL[h,4] <- ms2[n,2]
      HNL[h,5] <- ms2[m,2]
      if(m==1){HNL[h,6] <- ms2[n,2]}
      if(m!=1){HNL[h,6] <- 0.5*(HNL[h,4]+HNL[h,5])} # average int is used to be the new int of HNL 
      h <- h+1
    }
  }
  HNL <- HNL[HNL[,1]>=HNL.threshold,] # HNL threshold
  
  # formula filter
  if(is.null(db[[i]]$Formula)==FALSE){
    if(nrow(HNL)>(nrow(ms2))){
      valid <- formulafilter(formula)
      for(m in nrow(ms2):nrow(HNL)){
        if(min(abs(HNL[m,1]-valid))>mz.tol){
          HNL[m,1] <- NA
        }
      }
      HNL <- HNL[complete.cases(HNL),]
    }
  }
  
  new.ms2 <- HNL[order(HNL[,1]),c(1,6)]
  new.ms2 <- new.ms2[new.ms2[,2]>0,]
  rownames(new.ms2) <- NULL
  new.db[[i]]$pspectrum <- as.data.frame(new.ms2)
}
filename <- paste0('low quality spectra indices_',substring(db.name,1,regexpr("\\.",db.name)-1),'.csv')
write.csv(badNo.vector,file = filename,row.names = FALSE)
write.msp(new.db, file='MCS library.msp',newFile = TRUE)

