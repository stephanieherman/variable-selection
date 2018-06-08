rm(list=ls())
mainDir <- "~/MS/data/pos"
newDir <- "170626"
unlink(file.path(mainDir,newDir))
dir.create(file.path(mainDir,newDir))
setwd(file.path(mainDir,newDir))

load(file.path(mainDir,"df.RData")) # load raw data
load(file.path(mainDir,"mz_cf.RData")) # load consensus mz
load(file.path(mainDir,"rt_cf.RData")) # load consensus rt
load(file.path(mainDir,"spikes.RData")) # load spiked indices
source('Functions.R') # source function scripts

################################################
#
#               Start of pipeline
#
################################################

### Catch output log
sink('log.txt')

### Remove the first five samples in analysis sequence
names(df[,c(1:5)])
df <- df[,-c(1:5)]
df[df==0] <- NA

### Extract and remove
mz_spikes <- mz_cf[spikes]
mz_cf <- mz_cf[-spikes]

rt_spikes <- rt_cf[spikes]
rt_cf <- rt_cf[-spikes]

ref_spikes <- df[spikes,]
df <- df[-spikes,]
cat("Initial number of features: ", nrow(df), "\n")

### Remove NA rows
ind = remove.narows(df)
if (any(ind)) {  
  df <- df[!ind, ]
  mz_cf <- mz_cf[!ind]
  rt_cf <- rt_cf[!ind]
  }

### Blankfiltering - removing contaminants
blanks <- specgrep(df,"Blank")
samples.bf <- removegrep(df, "Blank")

to.remove <- BlankFilter(blanks,samples.bf,0.01)
samples.bf <- samples.bf[-to.remove,]
mz_cf <- mz_cf[-to.remove]
rt_cf <- rt_cf[-to.remove]
cat("Blankfiltering", "\n", "Removed by BlankFilter: ",length(to.remove), "\n")
rm(to.remove)

### Extract QC samples
qcs <- specgrep(samples.bf,"QC")
samples.bf.cov <- removegrep(samples.bf,"QC")

### Dilution based filtering
d <- specgrep(samples.bf.cov,"D")
samples.bf.cov.dil <- removegrep(samples.bf.cov,"D")

conc <- c(0.5,1,2,4,8,16,32)

# Calculate correlation between QC concentration and feature
c = list()
for (i in 1:dim(d)[1]) {
  x <- as.numeric(d[i,])
  y <- conc
  ind <- which(is.na(x))
  if (length(ind) > 0) {
    x <- x[-ind]
    y <- y[-ind]
  }
  if (length(x) > 2) {
    c[[i]] = cor.test(x,y,use="p")
  } else {
    c[[i]] = 0
  }
  
}

# Find features with an absolute correlation higher than 0.5
significant = list()
to.keep = c()
n = 1
for (j in 1:length(c)) {
  if (length(c[[j]])>1){
    if (c[[j]]$p.value < 0.1 && !is.na(c[[j]]$p.value)) {
      if (c[[j]]$estimate < -0.5 | c[[j]]$estimate > 0.5) {
        significant[[n]] <- c[[j]]
        to.keep[n] <- j
        n = n + 1
      }
    }
  }
}
cat("Dilution based filtering", "\n", "Features which correlate with the dilution serie: ", length(to.keep), "\n")
samples.bf.cov.dil = samples.bf.cov.dil[to.keep,]
mz_cf <- mz_cf[to.keep]
rt_cf <- rt_cf[to.keep]
rm(to.keep)

### log2 transformation
cat("log2 transformation", "\n")
samples.bf.cov.dil.log = log2(samples.bf.cov.dil)

### Outlier exclusion
cat("Outlier exclusion", "\n")

TIC <- colSums(samples.bf.cov.dil.log, na.rm=T)
par(mfrow = c(1,2))
boxplot(TIC, ylab = "Extracted peaks area")
hist(TIC, breaks = 100, main = "", xlab = "Extracted peaks area")

samples.bf.cov.dil.log.oe <- samples.bf.cov.dil.log[,!(colnames(samples.bf.cov.dil.log) 
                                                       %in% c(colnames(samples.bf.cov.dil.log)[TIC<mean(TIC)*0.7]))]
cat("Number of samples: ", ncol(samples.bf.cov.dil.log), "\n")
cat("Non-outlier samples: ", ncol(samples.bf.cov.dil.log.oe),"\n")

TIC <- colSums(samples.bf.cov.dil.log.oe, na.rm=T)
par(mfrow = c(1,2))
boxplot(TIC, ylab = "Extracted peaks area")
hist(TIC, breaks = 100, main = "", xlab = "Extracted peaks area")

### Remove NA rows
cat("Remove NA rows", "\n")
ind = remove.narows(samples.bf.cov.dil.log.oe)
if (any(ind)) {  
  samples.bf.cov.dil.log.oe <- samples.bf.cov.dil.log.oe[!ind, ]
  mz_cf <- mz_cf[!ind]
  rt_cf <- rt_cf[!ind]
  cat("Number of NA rows: ", length(which(ind)), "\n")
}

### Normalize to reference a.k.a. spiked-ins

# Extract remaining samples and log2 tranform the spiked-ins intensities
cat("Normalize to reference spikes", "\n")
ref <- ref_spikes[,names(samples.bf.cov.dil.log.oe)]
ref <- log2(ref)

samples.bf.cov.dil.log.oe.norm <- norm2ref(samples.bf.cov.dil.log.oe, ref)

par(mfrow = c(2,1))
boxplot(samples.bf.cov.dil.log.oe, main = "Before normalization")
boxplot(samples.bf.cov.dil.log.oe.norm, main = "After normalization")


### Merge replicates
cat("Merge replicates", "\n")

substrRight <- function(x){
  substr(x, nchar(x)-1, nchar(x))
}

samplesIDs <- gsub("_Rep1", "", names(samples.bf.cov.dil.log.oe.norm))
samplesIDs <- gsub("_Rep2", "", samplesIDs)
names <- unique(samplesIDs)
samplesIDs <- unlist(lapply(names,substrRight))

y = specgrep(samples.bf.cov.dil.log.oe.norm,paste(samplesIDs[1], "_", sep=""))
c = cor(y, use="complete.obs")[1,2]
corr = c
y = apply(y, 1, mean, na.rm = TRUE)
samples.bf.cov.dil.log.oe.norm.merged <- y

for (i in 2:length(samplesIDs)) {
  y = specgrep(samples.bf.cov.dil.log.oe.norm,paste(samplesIDs[i], "_", sep=""))
  if (dim(y)[2] > 2) {
    warning("More than two samples were grepped")
  }
  c = cor(y, use="complete.obs")[1,2]
  if (c > 0.80) {
    corr = c(corr,c)
    y = apply(y, 1, mean, na.rm=TRUE)
    samples.bf.cov.dil.log.oe.norm.merged <- cbind(samples.bf.cov.dil.log.oe.norm.merged,y)
  } else {
    corr = c(corr,c)
    cat("Correlation between duplicates in sample", samplesIDs[i], "had a correlation of: ", c, "\n")
    names <- names[-i]
  }
}
samples.bf.cov.dil.log.oe.norm.merged <- data.frame(samples.bf.cov.dil.log.oe.norm.merged)
names(samples.bf.cov.dil.log.oe.norm.merged) <- names
cat("Minimum duplicate correlation: ", min(corr), "\n")
cat("Average duplicate correlation: ", mean(corr), "\n")

### Albumin ratio filtering

#load albumindata
albumin <- read.table("~/Projects/MS/Albuminkvot.csv", sep=";", header=T)
row.names(albumin) < -albumin$Provtagningsidentitet
albumin <- albumin[names(samples.bf.cov.dil.log.oe.norm.merged),]

correlating <- data.frame(matrix(NA,nrow=nrow(samples.bf.cov.dil.log.oe.norm.merged),ncol=3))
names(correlating) <- c("ind", "corr","p-value")
row.names(correlating) <- row.names(samples.bf.cov.dil.log.oe.norm.merged)

for (i in 1:nrow(samples.bf.cov.dil.log.oe.norm.merged)) {
  alb <- albumin$Albuminkvot
  f <- as.numeric(samples.bf.cov.dil.log.oe.norm.merged[i,])
  ind <- unique(c(which(is.na(alb)), which(is.na(f))))
  if (length(ind)>0) {
    alb <- alb[-ind]
    f <- f[-ind]
  }
  if (length(f) > 15) {
    c <- cor.test(f,alb,method = "spearman")
    if (c$p.value < 0.05 & abs(c$estimate) > 0.5) {
      correlating[i,1] <- i
      correlating[i,2] <- c$estimate
      correlating[i,3] <- c$p.value
    }
  }
}

# remove na rows
ind = remove.narows(correlating)
if (any(ind)) {  
  correlating <- correlating[!ind, ]
}

cat("Features correlating with the albumin ratio: ", length(correlating$ind), "\n")
samples.bf.cov.dil.log.oe.norm.merged.alb <- samples.bf.cov.dil.log.oe.norm.merged[-correlating$ind,]
df <- samples.bf.cov.dil.log.oe.norm.merged.alb

# coverage cutoff
df <- df[coverage(df,0.75),]

save.image("workspace_processed.RData")
save(df, file="processed_Msdata.RData")
save(mz_cf, file="mz_cf.RData")
save(rt_cf, file="rt_cf.RData")
cat("Final number of features: ", nrow(df), "\n")
sink()
