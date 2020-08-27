#########################################
# This script is to do the batch McSearch against HNL library (mzXML input).
# Shipei Xing, Aug 20, 2020
# Copyright @ The University of British Columbia

#########################################
# Parameter setting
data.path <-'E:/McSearch/R_script'                    # Data path for "batch_search" folder. Use "/" instead of "\".
db.name <- 'Fiehn HILIC_HNL library.msp'              # The MS/MS spectral database name ending with ".msp"
mzXMLfile.name <- 'input_batch_search.mzXML'          # The name of input mzXML file.
pre.tol <- 0.01                                       # The mass tolerance for precursor search.
mass.range <- 200                                     # Precursor searching mass range, default 200 m/z. 
mz.tol <- 0.01                                        # The mass tolerance used for HNL matching (default 0.01 m/z).
topNo <- 100                                          # The maximum number of HNL peaks used for spectral matching, default 100.
HNL.threshold <- 36                                   # The minimum mass for HNL values, default 36 Da.
ion.mode <- 'n'                                       # Ionization mode: 'p' for positive, 'n' for negative

#########################################


# input: mzXML file
# output: a folder of csv files
## directed CCS score (library compounds as templates)

library('readMzXmlData')
library('clue')
library('metaMS')
library("CHNOSZ")
library("ChemmineR") 
library('ChemmineOB')
library("fmcsR")
library('stringr')
library('dplyr')
library('Rdisop')

setwd(data.path)
#input mzXML
data <- readMzXmlFile(mzXMLfile.name)

setwd(paste0(data.path,'/files'))
#input smiles.db, msp library
if(db.name=='MoNA_HNL library.msp'){smiles.db <- read.csv('smiles_db.csv',stringsAsFactors = FALSE)}
database <- read.msp(db.name, only.org = FALSE,
                     org.set = c('C','H','N','O','P','S','F','Cl','Br','I'), noNumbers = NULL)
bad.spectra <- read.csv(paste0('low quality spectra indices_',substring(db.name,1,regexpr("\\.",db.name)-1),'.csv'),stringsAsFactors = FALSE)
bad.No <- bad.spectra[,1]

setwd(paste0(data.path,'/biotrans'))
br.sdflist <- list.files(pattern= ".sdf")
bp <- read.csv('biotrans-plus.csv',stringsAsFactors = FALSE)
bm <- read.csv('biotrans-minus.csv',stringsAsFactors = FALSE)
brlist <- read.csv('biotrans mass change list.csv',stringsAsFactors = FALSE)
bruniquemass <- unique(brlist[,15])
#read and store the biotrans sdf info in variables called 'brsdf.x', x is the biotrans No.
for(m in 1:nrow(bp)){
  sdf.seq <- which(paste(m,'.sdf',sep='')==br.sdflist)
  if(length(sdf.seq)>0){
    assign(paste0('brsdf.',m),read.SDFset(br.sdflist[sdf.seq])) #read SDF
  }
}

# element table
element.table <- data.frame(matrix(0,ncol=2,nrow = 10))
colnames(element.table) <- c('element','number')
element.table[,1] <- c('C','H','N','O','P','S','F','Cl','Br','I')

CSS.score <- function(HNL.q, ms2.l){
  if(nrow(ms2.l)>topNo){
    ms2.l <- ms2.l[ms2.l[,2] > sort(ms2.l[,2],decreasing=TRUE)[topNo+1],]
  }
  # create HNL matrix for Hungarian algorithm
  # create HNL alignment
  HNL.alignment <- data.frame(matrix(ncol=4))
  colnames(HNL.alignment) <- c('HNL.q','int.q','HNL.l','int.l')
  for(m in 1:nrow(HNL.q)){
    mz.diff <- abs(HNL.q[m,1]-ms2.l[,1])
    if(min(mz.diff)<= mz.tol){
      HNL.alignment.individual <- cbind(HNL.q[m,1],HNL.q[m,6],ms2.l[mz.diff<=mz.tol,1],ms2.l[mz.diff<=mz.tol,2])
      colnames(HNL.alignment.individual) <- colnames(HNL.alignment)
      HNL.alignment <- rbind(HNL.alignment,HNL.alignment.individual)
    }
  }
  HNL.alignment <- HNL.alignment[complete.cases(HNL.alignment),]
  
  CSS <- 0
  mp <- 0
  if(nrow(HNL.alignment)>0){
    uniqueHNL.q <- unique(HNL.alignment[,1])
    uniqueHNL.l <- unique(HNL.alignment[,3])
    max.length <- max(length(uniqueHNL.q),length(uniqueHNL.l))
    matrix <- data.frame(matrix(0,ncol=(max.length+1),nrow=(max.length+1)))
    # fill the first row and column with HNL values
    matrix[2:nrow(matrix),1] <- c(uniqueHNL.q ,rep(0,(max.length-length(uniqueHNL.q))))
    matrix[1,2:nrow(matrix)] <- c(uniqueHNL.l ,rep(0,(max.length-length(uniqueHNL.l))))
    # fill in HNL.i * HNL.j
    for(m in 1:nrow(HNL.alignment)){
      matrix[matrix[,1]==HNL.alignment[m,1], matrix[1,]==HNL.alignment[m,3]] <- 
        HNL.alignment[m,2] * HNL.alignment[m,4]
    }
    if(length(uniqueHNL.q) > length(uniqueHNL.l)){matrix[2:nrow(matrix),(length(uniqueHNL.l)+2):ncol(matrix)] <- 0}
    if(length(uniqueHNL.q) < length(uniqueHNL.l)){matrix[(length(uniqueHNL.q)+2):nrow(matrix),2:ncol(matrix)] <- 0}
    matrix.B <-as.matrix(matrix[2:nrow(matrix),2:nrow(matrix)])
    #LSAP problem in 'clue' library (Hungarian algorithm)
    optimal <- solve_LSAP(matrix.B, maximum = TRUE)
    # calculate CSS score
    sum.q <- 0
    for(m in 1:max.length){
      CSS <- CSS + matrix.B[m,optimal[m]]
      if(matrix.B[m,optimal[m]]>0){
        sum.q <- sum.q + max(HNL.q[HNL.q[,1]==matrix[m+1,1],6])^2
        mp <- mp + 1
      }
    }
    CSS <- CSS/(sum.q*sum(ms2.l[,2]^2))^0.5
  }
  CSSreturn <- list(CSS,mp,nrow(ms2.l))
  return(CSSreturn)
}


for(i in 1:length(data)){
  # mslevel = 2
  if(data[[i]]$metaData$msLevel !=2) next
  # ms2
  if(length(data[[i]]$spectrum$mass)==0)next
  
  # information of the query spectrum
  premass.Q <- data[[i]]$metaData$precursorMz
  ms2.Q <- as.data.frame(cbind(data[[i]]$spectrum$mass, data[[i]]$spectrum$intensity))
  rt.Q <- data[[i]]$metaData$retentionTime

  # mass.Q
  if(ion.mode %in% c('P','p')){ mass.Q <- premass.Q - 1.007276}
  if(ion.mode %in% c('N','n')){ mass.Q <- premass.Q + 1.007276}
  
  # relative int
  ms2.Q[,2] <- 100*ms2.Q[,2]/max(ms2.Q[,2])
  # 1% threshold
  ms2.Q <- ms2.Q[ms2.Q[,2] >= 1,]
  # exclude after premass
  ms2.Q <- ms2.Q[ms2.Q[,1]<=(premass.Q+mz.tol),]
  # top 30 peaks
  if(nrow(ms2.Q)>30){ms2.Q <- ms2.Q[ms2.Q[,2] > sort(ms2.Q[,2],decreasing=TRUE)[31],]}
  # square root transformation
  ms2.Q[,2] <- sqrt(ms2.Q[,2])
  # introduce pseudo m/z '+1/-1' with int 0
  if(ion.mode %in% c('P','p')){
    ms2.Q <- rbind(c(1.007276,0),ms2.Q)
    adduct <- '[M+H]+'
  }
  if(ion.mode %in% c('N','n')){
    ms2.Q <- rbind(c(-1.007276,0),ms2.Q)
    adduct <- '[M-H]-'
  }
  # if premass is not in the MS2 spectrum, add (premass, 0)
  if(min(abs(ms2.Q[,1]-premass.Q)) > mz.tol){ms2.Q <- rbind(ms2.Q,c(premass.Q,0))}
  # HNL matrix
  HNL.Q <- data.frame(matrix(ncol=6))
  colnames(HNL.Q) <- c('HNL','mz.a','mz.b','int.a','int.b','HNL.int')
  h <- 1
  for(m in 1:(nrow(ms2.Q)-1)){
    for(n in (m+1):nrow(ms2.Q)){
      HNL.Q[h,1] <- ms2.Q[n,1] - ms2.Q[m,1]
      HNL.Q[h,2] <- ms2.Q[n,1]
      HNL.Q[h,3] <- ms2.Q[m,1]
      HNL.Q[h,4] <- ms2.Q[n,2]
      HNL.Q[h,5] <- ms2.Q[m,2]
      if(m==1){HNL.Q[h,6] <- ms2.Q[n,2]} # original fragment ions
      if(m!=1){HNL.Q[h,6] <- 0.5*(HNL.Q[h,4]+HNL.Q[h,5])} # average int is used to be the new int of HNL 
      h <- h+1
    }
  }
  HNL.Q <- HNL.Q[HNL.Q[,1]>=HNL.threshold,] # HNL threshold
  
  # topNo HNL peaks of query compound are used
  HNL.Q.1st <- HNL.Q[1:(nrow(ms2.Q)-1),]
  HNL.Q.2nd <- HNL.Q[nrow(ms2.Q):nrow(HNL.Q),]
  if(nrow(HNL.Q.2nd)>topNo){
    HNL.Q.2nd <- HNL.Q.2nd[HNL.Q.2nd[,6] > sort(HNL.Q.2nd[,6],decreasing=TRUE)[topNo+1],]
    HNL.Q <- rbind(HNL.Q.1st, HNL.Q.2nd)
  }
  
  # calculate the CSS scores with the stds(|premass diff| <= mass.range) in the database
  score.matrix <- as.data.frame(matrix(ncol=9))
  colnames(score.matrix) <- c('CSS intensity score','matched HNL No.','std HNL No.','name','formula','database No.','SMILES','InChIKey','ion mode')
  h <- 1
  for(l in 1:length(database)){
    if(is.element(l,bad.No)) next
    if(is.null(database[[l]]$PrecursorMZ)==FALSE){if(abs(database[[l]]$PrecursorMZ-premass.Q) > mass.range) next}
    
    formula.L <- database[[l]]$Formula
    if(grepl('\\[',formula.L)){
      formula.L <- substring(formula.L,regexpr("\\[",formula.L)+1,regexpr("\\]",formula.L)-1)
    }
    if(grepl('\\[',formula.L)==FALSE){
      if(grepl('\\+',formula.L) | grepl('\\-',formula.L)){
        formula.L <- substring(formula.L,1,nchar(formula.L)-1)
      }
    }
    
    a <- element.table
    a.table <- count.elements(formula.L)
    for(m in 1:length(a.table)){a[a[,1]==names(a.table)[m],2] <- a.table[m]}
    mass.L <- getMolecule(formula.L, z=0)$exactmass
    mass.list <- mass.L + bruniquemass - mass.Q
    if(min(abs(mass.list)) > pre.tol) next
    
    name.L <- database[[l]]$Name
    
    if(is.null(database[[l]]$Ion_mode)==FALSE){ionmode.L <- database[[l]]$Ion_mode}
    if(is.null(database[[l]]$Ion_mode)){
      if(is.null(database[[l]]$Precursor_type)==FALSE){
        str <- substr(database[[l]]$Precursor_type,nchar(database[[l]]$Precursor_type),nchar(database[[l]]$Precursor_type))
        if(str=="+"){ionmode.L <- 'P'}
        if(str=="-"){ionmode.L <- 'N'}
      }
      if(is.null(database[[l]]$Precursor_type)){ionmode.L <- 'Unknown'}
    }
    
    if(grepl("computed SMILES=", database[[l]]$Comments)){
      a <- substring(database[[l]]$Comments, regexpr("computed SMILES=", database[[l]]$Comments) + 16)
      smiles.L <- strsplit(a, '\"')[[1]][1]
    }
    if(grepl("computed SMILES=", database[[l]]$Comments)==FALSE){
      a <- substring(database[[l]]$Comments, regexpr("SMILES=", database[[l]]$Comments) + 7)
      smiles.L <- strsplit(a, '\"')[[1]][1]
    }
    
    if(is.null(database[[l]]$InChIKey)==FALSE){inchikey.L <- database[[l]]$InChIKey}
    if(is.null(database[[l]]$InChIKey)){inchikey.L <- paste0('No InChIKey info:',l)}
    
    ms2.L <- as.data.frame(database[[l]]$pspectrum)
    ms2.L[,2] <- 10*ms2.L[,2]/max(ms2.L[,2])
    
    CSS.list <- CSS.score(HNL.Q,ms2.L)
    score.matrix[h,1] <- as.numeric(CSS.list[1])
    score.matrix[h,2] <- as.numeric(CSS.list[2])
    score.matrix[h,3] <- as.numeric(CSS.list[3])
    score.matrix[h,4] <- name.L
    score.matrix[h,5] <- formula.L
    score.matrix[h,6] <- l
    if(db.name=='MoNA_HNL library.msp'){score.matrix[h,7] <- smiles.db[l,2]}
    if(db.name!='MoNA_HNL library.msp'){score.matrix[h,7] <- smiles.L}
    score.matrix[h,8] <- inchikey.L
    score.matrix[h,9] <- ionmode.L
    h <- h + 1 
  }
  score.matrix <- score.matrix[complete.cases(score.matrix),]
  score.matrix <- score.matrix[score.matrix[,1]>0,] # score > 0
  score.matrix <- score.matrix[score.matrix[,3]>1,] # std.HNL > 1
  score.matrix <- score.matrix[score.matrix[,2]>0,] # mp > 0
  score.matrix <- score.matrix[order(-(70*(score.matrix[,2]/topNo)/(0.5*log10(100*score.matrix[,3]/topNo))+5*score.matrix[,1])),]
  if(nrow(score.matrix) > 500){score.matrix <- score.matrix[1:500,]}
  score.matrix <- score.matrix[complete.cases(score.matrix),]
  if(nrow(score.matrix)== 0) next
  # same structure removal
  for(m in 1:(nrow(score.matrix)-1)){ 
    for(n in (m+1):nrow(score.matrix)){
      if(score.matrix[n,7]==score.matrix[m,7]|str_to_lower(score.matrix[n,4])==str_to_lower(score.matrix[m,4])|score.matrix[n,8]==score.matrix[m,8]){ # same calculated SMILES or name or inchikey
        score.matrix[n,1] <- NA
      }
    }
  }
  score.matrix <- score.matrix[complete.cases(score.matrix),]
  if(nrow(score.matrix) > 100){score.matrix <- score.matrix[1:100,]}
  
  
  # biotrans part
  output <- as.data.frame(matrix(ncol=20))
  colnames(output) <- c(colnames(score.matrix),'Adduct type','Heavy atom No.',"Reaction.1","Description.1","Reaction1_No.","Reaction.2","Description.2",
                        "Reaction2_No.","Formula.change","Mass error",'Final formula')
  for(p in 1:nrow(score.matrix)){
    #std sdf
    std.sdf <- smiles2sdf(score.matrix[p,7])
    
    a <- element.table
    a.table <- count.elements(score.matrix[p,5])
    for(m in 1:length(a.table)){a[a[,1]==names(a.table)[m],2] <- a.table[m]}
    std.mass <- getMolecule( score.matrix[p,5], z = 0)$exactmass
    options(digits = 6)
    #if same mass, dont go through following biotrans
    if(abs(std.mass-mass.Q)<=pre.tol){
      output.individual <- cbind(score.matrix[p,],adduct,matrix(0,ncol=8,nrow = 1),abs(std.mass-mass.Q),score.matrix[p,5])
      colnames(output.individual) <- colnames(output)
      output <- rbind(output, output.individual)
      next
    }
    
    result <- cbind(brlist,0,0)
    colnames(result)[17:18] <- c('Mass error','Final formula')
    result[,17] <- std.mass + brlist[,15] - mass.Q
    result <- result[abs(result[,17]) < pre.tol,]
    if(nrow(result)==0)next
    
    ## valid biotrans No. vector
    br.l <- vector(mode='numeric')
    h <- 1
    nonbr.l <- vector(mode='numeric')
    k <- 1
    for(m in 1:nrow(result)){
      #requirement element table
      b.req <- count.elements(result[m,1])
      b <- element.table
      for(n in 1:length(b.req)){b[b[,1]==names(b.req)[n],2] <- b.req[n]}
      # formula filter to see whether a biotrans can possibly occur
      if((all(a[,2]>=b[,2])==FALSE) | (all(a[,2]==b[,2])) | (all(a[,2]>=b[,2]) & a[1,2]==b[1,2]) ){result[m,1] <- NA}
      if((all(a[,2]>=b[,2])) & (all(a[,2]==b[,2])==FALSE)){
        # element change table
        c.change <- count.elements(result[m,2])
        c <- element.table
        for(n in 1:length(c.change)){c[c[,1]==names(c.change)[n],2] <- c.change[n]}
        # new element table
        new.et <- element.table
        new.et[,2] <- a[,2] + c[,2]
        new.formula <- character()
        for(n in 1:10){if(new.et[n,2]!=0){new.formula <- paste0(new.formula,new.et[n,1],new.et[n,2])}}
        result[m,18] <- new.formula
        #seven golden rules
        if(new.et[1,2]==0){result[m,1] <- NA}
        if(new.et[1,2]!=0 & (new.et[2,2]/new.et[1,2]>3.1|new.et[2,2]/new.et[1,2]<0.2|new.et[3,2]/new.et[1,2]>1.3
                             |new.et[4,2]/new.et[1,2]>1.2|new.et[5,2]/new.et[1,2]>0.3|new.et[6,2]/new.et[1,2]>0.8
                             |new.et[7,2]/new.et[1,2]>1.5|new.et[8,2]/new.et[1,2]>0.8|new.et[9,2]/new.et[1,2]>0.8)){result[m,1] <- NA}
        if(result[m,8]==0){
          ## structural similarity filter
          if(result[m,7] <= nrow(bp))next 
          if(is.element(result[m,7],br.l)|is.element(result[m,7],nonbr.l)) next
          if(exists(paste0('brsdf.',result[m,7]-nrow(bp)))){
            br.sdf <- get(paste0('brsdf.',result[m,7]-nrow(bp)))
            bond.tol <- bp[result[m,7]-nrow(bp),7]
            oc <- fmcs(std.sdf[[1]],br.sdf,bu=bond.tol,matching.mode='static')[['stats']][["Overlap_Coefficient"]]
            oc.limit <- (result[m,3]-1)/result[m,3]
            if(oc >= oc.limit){
              br.l[h] <- result[m,7]
              h <- h + 1
            }
            if(oc < oc.limit){
              nonbr.l[k] <- result[m,7]
              k <- k + 1
              result[m,1] <- NA
            }
          }
        }
        #2nd requirement element table
        if(result[m,8]!=0){
          d.req <- count.elements(result[m,8])
          d <- element.table
          for(n in 1:length(d.req)){d[d[,1]==names(d.req)[n],2] <- d.req[n]}
          # formula filter to see whether a biotrans can possibly occur
          if(all(new.et[,2]>=d[,2])==FALSE){result[m,1] <- NA}
          if(all(new.et[,2]>=d[,2])){
            # element change table
            e.change <- count.elements(result[m,16])
            e <- element.table
            for(n in 1:length(e.change)){e[e[,1]==names(e.change)[n],2] <- e.change[n]}
            # final formula
            final.et <- element.table
            final.et[,2] <- a[,2] + e[,2]
            new.formula <- character()
            for(n in 1:10){if(final.et[n,2]!=0){new.formula <- paste0(new.formula,final.et[n,1],final.et[n,2])}}
            result[m,18] <- new.formula
            ##seven golden rules
            if(final.et[1,2]==0){result[m,1] <- NA}
            if(final.et[1,2]!=0 & (final.et[2,2]/final.et[1,2]>3.1|final.et[2,2]/final.et[1,2]<0.2|final.et[3,2]/final.et[1,2]>1.3
                                   |final.et[4,2]/final.et[1,2]>1.2|final.et[5,2]/final.et[1,2]>0.3|final.et[6,2]/final.et[1,2]>0.8
                                   |final.et[7,2]/final.et[1,2]>1.5|final.et[8,2]/final.et[1,2]>0.8|final.et[9,2]/final.et[1,2]>0.8)){result[m,1] <- NA}
            ## structural similarity filter
            if((result[m,7] <= nrow(bp)) & (result[m,14]<= nrow(bp)))next 
            if(is.element(result[m,7],c(br.l,nonbr.l)) & is.element(result[m,14],c(br.l,nonbr.l))) next
            if((result[m,7] > nrow(bp)) & (is.element(result[m,7],c(br.l,nonbr.l))==FALSE)){
              if(exists(paste0('brsdf.',result[m,7]-nrow(bp)))){
                br.sdf <- get(paste0('brsdf.',result[m,7]-nrow(bp)))
                bond.tol <- bp[result[m,7]-nrow(bp),7]
                oc <- fmcs(std.sdf[[1]],br.sdf,bu=bond.tol,matching.mode='static')[['stats']][["Overlap_Coefficient"]]
                oc.limit <- (result[m,3]-1)/result[m,3]
                if(oc >= oc.limit){
                  br.l[h] <- result[m,7]
                  h <- h + 1
                }
                if(oc < oc.limit){
                  nonbr.l[k] <- result[m,7]
                  k <- k + 1
                  result[m,1] <- NA
                }
              }
            }
            if((result[m,14] > nrow(bp)) & (is.element(result[m,14],c(br.l,nonbr.l))==FALSE)){
              if(exists(paste0('brsdf.',result[m,14]-nrow(bp)))){
                br.sdf <- get(paste0('brsdf.',result[m,14]-nrow(bp)))
                bond.tol <- bp[result[m,14]-nrow(bp),7]
                oc <- fmcs(std.sdf[[1]],br.sdf,bu=bond.tol,matching.mode='static')[['stats']][["Overlap_Coefficient"]]
                oc.limit <- (result[m,10]-1)/result[m,10]
                if(oc >= oc.limit){
                  br.l[h] <- result[m,14]
                  h <- h + 1
                }
                if(oc < oc.limit){
                  nonbr.l[k] <- result[m,14]
                  k <- k + 1
                  result[m,1] <- NA
                }
              }
            }
          }
        }
      }
    }
    result <- result[complete.cases(result),]
    if(nrow(result)>0){
      for(m in 1:nrow(result)){
        if(is.element(result[m,7],nonbr.l) | is.element(result[m,14],nonbr.l)){result[m,1] <- NA}
      }
    }
    result <- result[complete.cases(result),]
    
    if(nrow(result)==0)next

    heavyatom <- result[,3] + result[,10]
    output.individual <- cbind(score.matrix[p,],adduct,heavyatom,result[,c(4,6,7,11,13,14,16,17,18)])
    colnames(output.individual) <- colnames(output)
    output <- rbind(output, output.individual)
  }
  output <- output[complete.cases(output),]
  
  if(nrow(output)>=1){
    # McSearch score
    output <- cbind(output,0,0)
    colnames(output)[21:22] <- c('modified matched ratio','McSearch score')
    output[,21] <- (output[,2]/topNo)/(0.5*log10(100*output[,3]/topNo))
    output[,22] <- 70*output[,21] + 5*output[,1] + 10*(1-abs(output[,19])/pre.tol) + 10*15/(output[,11]+15)
    for(j in 1:nrow(output)){
      if(output[j,15]==0){output[j,22] <- output[j,22] + 2}
      low <- c(48:51,107:110)
      if(output[j,11] %in% low==FALSE & output[j,14] %in% low==FALSE){output[j,22] <- output[j,22] + 3}
    }
    output <- output[order(-output[,22]),]
    output <- output[complete.cases(output),]
    
    if(nrow(output)==1) {output.rank <- cbind(output,1)}
    if(nrow(output)>1){
      # McSearch rank
      for(m in 1:nrow(output)){if(output[m,8]==0){output[m,8] <- output[m,6]}}
      row.names(output) <- 1:nrow(output)
      output.rank.individual <- data.frame(matrix(ncol=23))
      colnames(output.rank.individual) <- c(colnames(output),'rank')
      output.rank <- output.rank.individual
      rank <- 1
      h <- 1
      unique.hit <- output[row.names(unique(output[,c(8,20)])),]
      for(m in 1:nrow(unique.hit)){
        output.rank.individual <- cbind(unique.hit[m,],rank)
        row.No <- as.numeric(row.names(unique(output[,c(8,20)]))[m])+1
        for(n in row.No:nrow(output)){
          if(row.No >= nrow(output)) next #for the last row
          if(output[n,8]==unique.hit[m,8] & output[n,20]==unique.hit[m,20]){
            output.rank.individual <- rbind(output.rank.individual,cbind(output[n,],rank))
          }
        }
        rank <- rank + 1
        output.rank <- rbind(output.rank, output.rank.individual)
      }
      output.rank <- output.rank[complete.cases(output.rank),]
    }
  }
  output.rank <- output.rank[,c(23,1:22)]
  
  setwd(paste0(data.path,'/output'))
  filename <- paste0('output_premass',premass.Q,'_rt',rt.Q,'.csv')
  write.csv(output.rank,file=filename,row.names = FALSE)
}


