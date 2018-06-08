rm(list=ls())
library(foreach)
library(doParallel)
library(ropls)
mainDir <- "~/Projects/MS/data/"
setwd(file.path(mainDir,"doParallelCV"))
neg <- "neg/170206"
pos <- "pos/170206"

# load data
load(file.path(file.path(mainDir,neg),"negDF.RData"))
load(file.path(file.path(mainDir,pos),"posDF.RData"))
load("crpData.rdata")




# ----------------------------------- Fix metabolites -----------------------------------

# merge positive and negative ion mode data
row.names(pos.df) <- paste(c(1:nrow(pos.df)),"pos", sep="")
row.names(neg.df) <- paste(c(1:nrow(neg.df)),"neg", sep="")
metabolites <- data.frame(rbind(pos.df,neg.df))

# positive
pos <- as.numeric(gsub("pos","",row.names(metabolites)[grep("pos",row.names(metabolites))]))
load ("C:\\Users\\stehe524\\Documents\\Projects\\MS\\data\\pos\\170626\\mz_cf.RData")
load ("C:\\Users\\stehe524\\Documents\\Projects\\MS\\data\\pos\\170626\\rt_cf.RData")
pos_cf <- data.frame(cbind(mz_cf[pos],rt_cf[pos]))
rm(mz_cf)
rm(rt_cf)

# negative
neg <- as.numeric(gsub("neg","",row.names(metabolites)[grep("neg",row.names(metabolites))]))
load ("C:\\Users\\stehe524\\Documents\\Projects\\MS\\data\\neg\\170626\\mz_cf.RData")
load ("C:\\Users\\stehe524\\Documents\\Projects\\MS\\data\\neg\\170626\\rt_cf.RData")
neg_cf <- data.frame(cbind(mz_cf[neg],rt_cf[neg]))

metabolites <- cbind(rbind(pos_cf,neg_cf),metabolites)
names(metabolites)[c(1,2)] <- c("mz","rt")

# match in silico fragmented IDs
ppmCal<-function(run,ppm) {
  return((run*ppm)/1000000)
}

IDs_pos <- read.csv("~/Projects/MS/data/MS2/Galaxy79-[metfrag_identification_results]_pos.csv")
IDs_neg <- read.csv("~/Projects/MS/data/MS2/Galaxy192-[metfrag_identification_results]_neg.csv")

identified <- c()
allinfo <- data.frame(matrix(NA,nrow=1,ncol=ncol(IDs_neg)))
names(allinfo) <- names(IDs_neg)
for (i in 1:nrow(metabolites)) {
  mz <- metabolites[i,1]
  rt <- metabolites[i,2]
  if (grepl("pos",row.names(metabolites)[i])) {
    ind_mz <- which(IDs_pos$parentMZ>mz-ppmCal(mz,10) & IDs_pos$parentMZ<mz+ppmCal(mz,10))
    ind_rt <- which(IDs_pos$parentRT>rt-60 & IDs_pos$parentRT<rt+60)
    ind <- intersect(ind_mz,ind_rt)
    
    if (length(ind)>0) {
      identified <- c(identified,i)
      temp <-IDs_pos[ind,]
      allinfo <- rbind(allinfo, temp[order(temp$FragmenterScore,decreasing=TRUE)[1:3],])
    }
  } else {
    ind_mz <- which(IDs_neg$parentMZ>mz-ppmCal(mz,10) & IDs_neg$parentMZ<mz+ppmCal(mz,10))
    ind_rt <- which(IDs_neg$parentRT>rt-60 & IDs_neg$parentRT<rt+60)
    ind <- intersect(ind_mz,ind_rt)
    
    if (length(ind)>0) {
      identified <- c(identified,i)
      temp <-IDs_neg[ind,]
      allinfo <- rbind(allinfo, temp[order(temp$FragmenterScore,decreasing=TRUE)[1:3],])
    }
  }
}

metabolites <- metabolites[identified,]

# correct for age
ind <- which(crp$Type=="SPMS")

metabolites_corr <- data.frame(matrix(NA,nrow=(nrow(metabolites)),ncol=ncol(metabolites)))
names(metabolites_corr)<-names(metabolites)
row.names(metabolites_corr)<-row.names(metabolites)
metabolites_corr$Type <- metabolites$Type
n <- 0
for (i in 2:ncol(metabolites))  {
  c <- cor.test(metabolites[-ind,i], crp$Age[-ind], na.action = "na.omit")
  if (c$p.value < 0.05) { # check if the corr is significant
    n <- n + 1
    model <- lm(metabolites[-ind,i] ~ crp$Age[-ind], na.action = "na.omit") # x = C*age
    metabolites_corr[,i] <- metabolites[,i]-model$coefficients[2]*crp$Age # x_corrected = x - C*age
  } else {
    metabolites_corr[,i] <- metabolites[,i]
  }
}

cat("Number of age corrected metabolites: ", n, "\n")
rm(metabolites)
metabolites <- metabolites_corr

# impute by the column mean
colmissing <- which(apply(is.na(metabolites), 2, any))
for (i in 1:length(colmissing)) {
  temp <- metabolites[,colmissing[i]] 
  temp[is.na(temp)] <- mean(temp,na.rm=T)
  metabolites[,colmissing[i]] <- temp
}

# extract transitioning patients
trans <- metabolites[c("Transition.01","Transition.02","Transition.03","Transition.04"),]
metabolites <- metabolites[!(row.names(metabolites) %in% c("Transition.01","Transition.02","Transition.03","Transition.04")),]

order < -names(metabolites)
metabolites$group <- ''
metabolites <- metabolites[,c("group",order)]


# ----------------------------------- Fix CRP -----------------------------------

# correct for age
crp_corr <- data.frame(matrix(NA,nrow=(nrow(crp)),ncol=ncol(crp)))
names(crp_corr)<-names(crp)
row.names(crp_corr)<-row.names(crp)
crp_corr$Type <- crp$Type
MRI <- c("Size_ventricle","Rmsize","TotalT1", "TotalT2", "TotalGd")
n=c()
for (i in 3:48)  {
  c <- cor.test(crp[-ind,i], crp$Age[-ind], na.action = "na.omit")
  if (names(crp)[i] %in% MRI) {
    crp_corr[,i] <- crp[,i]
  } else if (c$p.value < 0.05) { # check if the corr is significant
    n = c(n,i)
    model <- lm(crp[-ind,i] ~ crp$Age[-ind], na.action = "na.omit") # x = C*age
    crp_corr[,i] <- crp[,i]-model$coefficients[2]*crp$Age  # x_corrected = x - C*age
  } else {
    crp_corr[,i] <- crp[,i]
  }
}

cat("Age corrected CRP variables: ", names(crp)[n], "\n")
rm(crp)
crp <- crp_corr[,!(names(crp_corr) %in% c("Age"))]

# impute by the column mean
colmissing <- which(apply(is.na(crp), 2, any))
for (i in 1:length(colmissing)) {
  temp <- crp[,colmissing[i]] 
  temp[is.na(temp)] <- mean(temp,na.rm=T)
  crp[,colmissing[i]] <- temp
}


# extract transitioning patients
trans <- crp[c("Transition.01","Transition.02","Transition.03","Transition.04"),]
crp <- crp[!(row.names(crp) %in% c("Transition.01","Transition.02","Transition.03","Transition.04")),]
table(crp$Type)

names<-names(crp)
crp$group <- ''
crp <- crp[,c("group",names)]

# ----------------------------------- Cross validation -----------------------------------

# define CV blocks
control <- c(1,1,2,2,3,3,4,4,5,5)
RRMS <- c(rep(1,6),rep(2,5),rep(3,5),rep(4,5),rep(5,5))
SPMS <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,4))

CV <- 10

# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer
registerDoParallel(cl)

ERmatrix <- foreach(j=1:CV, .combine=cbind, .packages = "ropls", .verbose = T) %dopar% {
  set.seed(j)
  
  # calcER - calculate class specific error rates
  calcER <- function(label, true, prediction) {
    loc <- grep(label,true)
    er <- length(which(prediction[loc]!=true[loc]))/length(loc)
    return(er)
  }
  
  newgroups <- c(sample(control),sample(RRMS),sample(SPMS))
  metabolites$group <- newgroups
  crp$group <- newgroups
  
  library('ropls')
  ER <- data.frame(matrix(NA, ncol= 5, nrow=20))
  for(i in 1:5) {
    
    # partition the data into training and testing set
    Metabo_test <- metabolites[which(metabolites$group==i),]
    CRP_test <- crp[which(crp$group==i),]

    Metabo_train <- metabolites[-which(metabolites$group==i),]
    CRP_train <- crp[-which(crp$group==i),]
    
    train <- cbind(Metabo_train[,-1],CRP_train[,-c(1,2)])
    test <- cbind(Metabo_test[,-1],CRP_test[,-c(1,2)])
    
    names(test) <- gsub("X","",names(test))
    
    # ----------------- Metabolomics models -----------------
    
    ## All vs All
    # build metabolomics model - all vs all
    data.plsda <- opls(Metabo_train[,-c(1,2)], Metabo_train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    metabo_AvsA <- names(sort(data.plsda@vipVn,decreasing=T)[1:10]) # extract top 10 variables
    
    # predict test set
    prediction <- predict(data.plsda, Metabo_test[,-c(1,2)])
    
    # calculate error rates
    SP.er <- calcER("SPMS",Metabo_test$Type,prediction)
    RR.er <- calcER("RRMS",Metabo_test$Type,prediction)
    con.er <- calcER("control",Metabo_test$Type,prediction)
    BER <- (SP.er+RR.er+con.er)/3
    
    ER[1,i] <- SP.er                                      ##OUTPUT
    ER[2,i] <- BER                                        ##OUTPUT
    
    ## RRMS vs SPMS
    # temporarily remove the control group
    tmp_Metabo_train <- Metabo_train
    tmp_Metabo_test <- Metabo_test
    
    tmp_Metabo_train <- tmp_Metabo_train[-which(tmp_Metabo_train$Type=="control"),]
    tmp_Metabo_test <- tmp_Metabo_test[-which(tmp_Metabo_test$Type=="control"),]
    
    # build metabolomics model - RRMS vs SPMS
    data.plsda <- opls(tmp_Metabo_train[,-c(1,2)], tmp_Metabo_train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    metabo_RvsS <- names(sort(data.plsda@vipVn,decreasing=T)[1:10]) # extract top 10 variables
    
    
    # predict test set
    prediction <- predict(data.plsda, tmp_Metabo_test[,-c(1,2)])
    
    # calculate error rates
    SP.er<- calcER("SPMS",tmp_Metabo_test$Type,prediction)
    RR.er<- calcER("RRMS",tmp_Metabo_test$Type,prediction)
    BER <- (SP.er+RR.er)/2
    
    ER[3,i] <- SP.er                                      ##OUTPUT
    ER[4,i] <- BER                                        ##OUTPUT
    
    
    # ----------------- Reduced metabolomics models -----------------
    
    top_metabo <- intersect(metabo_AvsA,metabo_RvsS) 
    
    if (length(top_metabo)>1)  {
      
      Metabo_train <- Metabo_train[,c("Type",top_metabo)]
      Metabo_test <- Metabo_test[,c("Type",top_metabo)]
      
      # build reduced metabolomics model - all vs all
      data.plsda <- opls(Metabo_train[,-1], Metabo_train$Type, printL=FALSE,
                         predI = 2, scaleC='standard', plotL=FALSE)
      
      # predict test set
      prediction <- predict(data.plsda, Metabo_test[,-1])
      
      # calculate error rates
      SP.er<- calcER("SPMS",Metabo_test$Type,prediction)
      RR.er<- calcER("RRMS",Metabo_test$Type,prediction)
      con.er <- calcER("control",Metabo_test$Type, prediction)
      BER <- (SP.er+RR.er+con.er)/3
      
      ER[5,i] <- SP.er                                      ##OUTPUT
      ER[6,i] <- BER                                        ##OUTPUT

      # build reduced metabolomics model - RRMS vs SPMS
      Metabo_train <- Metabo_train[-which(Metabo_train$Type=="control"),]
      Metabo_test <- Metabo_test[-which(Metabo_test$Type=="control"),]
      
      data.plsda <- opls(Metabo_train[,-1], Metabo_train$Type, printL=FALSE,
                         predI = 2, scaleC='standard', plotL=FALSE)
      
      # predict test set
      prediction <- predict(data.plsda, Metabo_test[,-1])
      
      # calculate error rates
      SP.er<- calcER("SPMS",Metabo_test$Type,prediction)
      RR.er<- calcER("RRMS",Metabo_test$Type,prediction)
      BER <- (SP.er+RR.er)/2
      
      ER[7,i] <- SP.er                                      ##OUTPUT
      ER[8,i] <- BER                                        ##OUTPUT
    }
    
    rm(Metabo_train)
    rm(Metabo_test)
    
    # ----------------- CRP models -----------------
    
    ## All vs All
    # build CRP model - all vs all
    data.plsda <- opls(CRP_train[,-c(1,2)], CRP_train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    crp_AvsA <- names(sort(data.plsda@vipVn,decreasing=T)[1:10]) # extract top 10 variables
    
    # predict test set
    prediction <- predict(data.plsda, CRP_test[,-c(1,2)])
    
    # calculate error rates
    SP.er<- calcER("SPMS",CRP_test$Type,prediction)
    RR.er<- calcER("RRMS",CRP_test$Type,prediction)
    con.er <- calcER("control",CRP_test$Type, prediction)
    BER <- (SP.er+RR.er+con.er)/3
    
    ER[9,i] <- SP.er                                      ##OUTPUT
    ER[10,i] <- BER                                       ##OUTPUT
    
    ## RRMS vs SPMS
    # temporarily remove the control group
    tmp_CRP_train <- CRP_train
    tmp_CRP_test <- CRP_test
    
    tmp_CRP_train <- tmp_CRP_train[-which(tmp_CRP_train$Type=="control"),]
    tmp_CRP_test <- tmp_CRP_test[-which(tmp_CRP_test$Type=="control"),]
    
    # build CRP model - RRMS vs SPMS
    data.plsda <- opls(tmp_CRP_train[,-c(1,2)], tmp_CRP_train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    crp_RvsS <- names(sort(data.plsda@vipVn,decreasing=T)[1:10]) # extract top 10 variables
    
    # predict test set
    prediction <- predict(data.plsda, tmp_CRP_test[,-c(1,2)])
    
    # calculate error rates
    SP.er<- calcER("SPMS",tmp_CRP_test$Type,prediction)
    RR.er<- calcER("RRMS",tmp_CRP_test$Type,prediction)
    BER <- (SP.er+RR.er)/2
    
    ER[11,i] <- SP.er                                     ##OUTPUT
    ER[12,i] <- BER                                       ##OUTPUT
    
    # -----------------  Reduced CRP models ----------------- 
    
    top_crp <- intersect(crp_AvsA,crp_RvsS)
    
    if (length(top_crp)>1) {
      
      CRP_train <- CRP_train[,c("Type",top_crp)]
      CRP_test <- CRP_test[,c("Type",top_crp)]
      
      # build reduced CRP model - all vs all
      data.plsda <- opls(CRP_train[,-1], CRP_train$Type, printL=FALSE,
                         predI = 2, scaleC='standard', plotL=FALSE)
      
      # predict test set
      prediction <- predict(data.plsda, CRP_test[,-1])
      
      # calculate error rates
      SP.er<- calcER("SPMS",CRP_test$Type,prediction)
      RR.er<- calcER("RRMS",CRP_test$Type,prediction)
      con.er <- calcER("control",CRP_test$Type, prediction)
      BER <- (SP.er+RR.er+con.er)/3
      
      ER[13,i] <- SP.er                                     ##OUTPUT
      ER[14,i] <- BER                                       ##OUTPUT
      
      # build reduced crp model - RRMS vs SPMS
      CRP_train <- CRP_train[-which(CRP_train$Type=="control"),]
      CRP_test <- CRP_test[-which(CRP_test$Type=="control"),]
      
      data.plsda <- opls(CRP_train[,-1], CRP_train$Type, printL=FALSE,
                         predI = 2, scaleC='standard', plotL=FALSE)
      
      # predict test set
      prediction <- predict(data.plsda, CRP_test[,-1])
      
      # calculate error rates
      SP.er<- calcER("SPMS",CRP_test$Type,prediction)
      RR.er<- calcER("RRMS",CRP_test$Type,prediction)
      BER <- (SP.er+RR.er)/2
      
      ER[15,i] <- SP.er                                     ##OUTPUT
      ER[16,i] <- BER                                       ##OUTPUT
      }
    
    
    # -----------------  CRPM models ----------------- 
    
    toptop <- c(top_metabo,top_crp)
    toptop <- c("Type",toptop)
    test <- test[,toptop]
    train <- train[,toptop]
    
    ## All vs All
    # build CRPM model - all vs all
    data.plsda <- opls(train[,-1], train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    # predict test set
    prediction <- predict(data.plsda, test[,-1])
    
    # calculate error rates
    SP.er<- calcER("SPMS",test$Type,prediction)
    RR.er<- calcER("RRMS",test$Type,prediction)
    con.er <- calcER("control",test$Type, prediction)
    BER <- (SP.er+RR.er+con.er)/3
    
    ER[17,i] <- SP.er                                     ##OUTPUT
    ER[18,i] <- BER                                       ##OUTPUT
    
    ## RRMS vs SPMS
    # remove control group
    train <- train[-which(train$Type=="control"),]
    test <- test[-which(test$Type=="control"),]
    
    # build CRPM model - RRMS vs SPMS
    data.plsda <- opls(train[,-1], train$Type, printL=FALSE,
                       predI = 2, scaleC='standard', plotL=FALSE)
    
    # predict test set
    prediction <- predict(data.plsda, test[,-1])
    
    # calculate error rates
    SP.er<- calcER("SPMS",test$Type,prediction)
    RR.er<- calcER("RRMS",test$Type,prediction)
    BER <- (SP.er+RR.er)/2
    
    ER[19,i] <- SP.er                                     ##OUTPUT
    ER[20,i] <- BER                                       ##OUTPUT
  }
  return(ER)
}
stopCluster(cl)

row.names(ERmatrix) <- c("spms_AvA_metabo","global_AvA_metabo",
                         "spms_RvS_metabo","global_RvS_metabo",
                         
                         "spms_AvA_metabo_reduced","global_AvA_metabo_reduced",
                         "spms_RvS_metabo_reduced","global_RvS_metabo_reduced",
                         
                         "spms_AvA_crp","global_AvA_crp",
                         "spms_RvS_crp","global_RvS_crp",
                         
                         "spms_AvA_crp_reduced","global_AvA_crp_reduced",
                         "spms_RvS_crp_reduced","global_RvS_crp_reduced",
                         
                         "spms_AvA_comb","global_AvA_comb",
                         "spms_RvS_comb","global_RvS_comb")
names(ERmatrix) <- 1:50

ERmatrix<-data.frame(t(ERmatrix)) 
ERs <- data.frame(rbind((colMeans(ERmatrix)),apply(ERmatrix,2,sd)))
