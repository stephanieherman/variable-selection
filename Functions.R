######################
#
#   Grep name in x
#
#####################
specgrep <- function(x,name) {
  x=x[,grep(name, names(x))]
  return(x)
}


######################
#
#   Remove name in x
#
#####################
removegrep <- function(x,name) {
  x=x[,-grep(name, names(x))]
  return(x)
}


######################
#
#   BlankFilter - find features which are "highly" present in blank (based on a ratio cutoff)
#   to.remove=AdvancedBlankFilter(Blanks,samples,0.01)
#
#####################
BlankFilter <- function(blanks, samples, cutoff) {
  blanks[is.na(blanks)] <- 0
  samples[is.na(samples)] <- 0
  
  blanks <- apply(blanks,1,median,na.rm=TRUE)
  samples <- apply(samples,1,max,na.rm=TRUE)
  
  to.remove <- which(blanks/samples >= cutoff)
  return(to.remove)
}


######################
#
#   coverage - compute features with coverage of cutoff or above
#   features=coverage(samples,1)
#
#####################
coverage <- function(intensities, coverage) {
  amount = dim(intensities)
  cutoff = amount[2]*coverage
  features=rep(0,amount[1])
  
  for (i in 1:amount[1]) {
    temp=intensities[i,]
    a=!is.na(temp)
    a=which(a)
    if (length(a)>=cutoff) {
      features[i]=i
    }
  }
  features=features[which(features>0)]
  return(features)
}


######################
#
#   extract.names - extract unique names
#   names=extract.names(samples)
#
#####################
extract.names <- function(x){
  names=gsub("intensity_","",names(x))
  names=gsub("_Rep1","",names)
  names=gsub("_Rep2","",names)
  names=gsub("_Rep3","",names)
  
  names=gsub("_2","",names)
  names=unique(names)
  return(names)
}


######################
#
#   remove.narows - removes rows/features which are completely missing (NA)
#   samples=remove.narows(samples)
#
#####################
remove.narows<-function(X)  {
  ind <- apply(X, 1, function(x) all(is.na(x)))
  return(ind)
}


######################
#
#   Outersect
#
#####################
outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

######################
#
#   coverage - compute features with coverage of cutoff or above
#   features=coverage(samples,0.75)
#
#####################
coverage <- function(intensities, coverage) {
  amount = dim(intensities)
  cutoff = amount[2]*coverage
  features=rep(0,amount[1])
  
  for (i in 1:amount[1]) {
    temp=intensities[i,]
    a=!is.na(temp)
    a=which(a)
    if (length(a)>=cutoff) {
      features[i]=i
    }
  }
  features=features[which(features>0)]
  return(features)
}


#####################
#
#   Normalize against references (spiked-ins)
#
#####################
norm2ref <- function(data,ref) {
  m <- rowMeans(ref)
  for (i in 1:ncol(data)) {
    factor <- c()
    for (j in 1:length(m)) {
      factor[j]<-ref[j,i]-m[j]
    }
    factor <- mean(factor)
    data[,i] <- as.matrix(data[,i])-rep(factor,nrow(data))
  }
  return(data)
}

