#####################################################################################
#SETTING ENVIRONMENT
#####################################################################################

#Load packages needes for the analysis
lapply(c("ggplot2", "psych", "RCurl", "irr", "nortest", "moments","nFactors",
         "psych","GPArotation", "gmodels","sqldf","gdata","sem"), library, character.only=T)


#####################################################################################
#IMPORTING DATA
#####################################################################################

#uploading data ------------------------------------------------------------------------
#Functions to pull the dara from the internet file 
#see http://goo.gl/mQwxO on how to get this link

data<-read.csv("/Users/joaovissoci/Desktop/mergeRC1.csv",header=T)


###########################################################################################
#EXPLORATORY DATA ANALYSIS
###########################################################################################
#Exploratory Data Anlysis
#dim(data)
str (data) #Wll give you the classification of each variable
#head(data) #Will show you the first 6 observations of the dataset
#tail(data) #Will show you the ast 6 observations of the dataset
#names(data) #Will show you the names for all the variables
summary(data)#This comand will provide a whole set of descriptive results for each variables
#ad.test() # Anderson-Darling test for normality
#skewness() #Will provide skweness analysis
#kurtosis() - 3 #Will provide kurtosis analysis
#qplot() # histogram plot
#pwr.anova.test(k = , n = , f = , sig.level = , power = )#Power analysis for variance analysis
#boxplot() #will provide a boxplot for the variables to analysis potential outliers

###########################################################################################
#DATA MANAGEMENT
###########################################################################################
#Excluding ST1, MD6, SH5 due ot missing data
data<-with(data,data.frame(hindex,noOfDocs,IDPT,IM,MR,
  AFR,BR,DR,GDP,HEP,HET))

summary(data)
data$hindex<-as.numeric(as.character(data$hindex))
data$noOfDocs<-as.numeric(as.character(data$noOfDocs))
data$IM<-as.numeric(as.character(data$IM))
data$MR<-as.numeric(as.character(data$MR))
data$AFR<-as.numeric(as.character(data$AFR))
data$BR<-as.numeric(as.character(data$BR))
data$DR<-as.numeric(as.character(data$DR))
data$GDP<-as.numeric(as.character(data$GDP))
data$HEP<-as.numeric(as.character(data$HEP))
data$HET<-as.numeric(as.character(data$HET))

str(data)
datanew<-na.omit(data)
str(datanew)

write.csv(data, "/Users/joaovissoci/Desktop/LODHealthAquamethData.csv")

#AI<-data.frame(data$AI1,data$AI2,data$AI3,data$AI4)
#dim(AI)
#SH<-data.frame(data$SH2,data$SH3,data$SH4)
#MD<-data.frame(data$MD1,data$MD2,data$MD3,data$MD4,data$MD5)
#EI<-data.frame(data$EI1,data$EI2)
#Health<-data.frame(data$H1,data$H2,data$H3)
#HI<-data.frame(SH,MD)

###########################################################################################
#TABLE 1: EFA
###########################################################################################
#Correlations

qgraph(cor(datanew,method="spearman"),cut=1,minimum=.20,layout="spring",
  directed=FALSE,label.scale=FALSE,edge.labels=TRUE,label.cex = 0.80)


#efadata<-remove.vars(data,names=c("Country","Year","Year_1","X"))
#efadata<-na.omit(efadata)
#str(efadata)

##Eigen-Values and Scree plot
par(mfrow=c(2,2)) #Command to configure the plot area for the scree plot graph
ev <- eigen(cor(datanew)) # get eigenvalues - insert the data you want to calculate the scree plot for
ev # Show eigend values
ap <- parallel(subject=nrow(datanew),var=ncol(datanew),rep=100,cent=.05) #Calculate the acceleration factor
nS <- nScree(ev$values) #Set up the Scree Plot 
plotnScree(nS) # Plot the ScreePlot Graph

#KMO
kmo = function(data){
  
  library(MASS)
  X <- cor(as.matrix(data))
  iX <- ginv(X)
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a)
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy
  
  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the
  # correlation matrix. That is the
  # negative of the partial correlations,
  # partialling out all other variables.
  
  kmo <- BB/(AA+BB)                     # overall KMO statistic
  
  # Reporting the conclusion
  if (kmo >= 0.00 && kmo < 0.50){
    test <- 'The KMO test yields a degree of common variance
    unacceptable for FA.'
  } else if (kmo >= 0.50 && kmo < 0.60){
    test <- 'The KMO test yields a degree of common variance miserable.'
  } else if (kmo >= 0.60 && kmo < 0.70){
    test <- 'The KMO test yields a degree of common variance mediocre.'
  } else if (kmo >= 0.70 && kmo < 0.80){
    test <- 'The KMO test yields a degree of common variance middling.'
  } else if (kmo >= 0.80 && kmo < 0.90){
    test <- 'The KMO test yields a degree of common variance meritorious.'
  } else {
    test <- 'The KMO test yields a degree of common variance marvelous.'
  }
  
  ans <- list(  overall = kmo,
                report = test,
                individual = MSA,
                AIS = AIS,
                AIR = AIR )
  return(ans)
  
}    # end of kmo()
kmo(datanew)

## EFA para o banco todo

health_efa<-with(datanew,data.frame(IDPT,IM,MR,
  AFR,BR,HEP))

efa_HO<-fa(health_efa,1,rotate="promax")


###########################################################################################
#SEM MODELS
###########################################################################################
#SEM Model C - Working but poor fit indicators
datanew<-remove.vars(datanew,names=c("DR"))
#attach(semdata)

#FINAL MODEL

semmodel<- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).
#Latent Variables
HealthOut->IDPT,efa14,NA
HealthOut->IM,efa11,NA
HealthOut->MR,efa12,NA
HealthOut->AFR,efa13,NA
HealthOut->BR,NA,1
HealthOut->HEP,efa16,NA
hindex->GDP,ob1,NA
HET->GDP,ob3,NA
hindex->HealthOut,lat12,NA
HET->HealthOut,lat14,NA
HealthOut<->HealthOut,laterro1,NA
GDP<->GDP,erro1,NA
hindex<->hindex,erro3,NA
noOfDocs<->noOfDocs,erro4,NA
HET<->HET,erro5,NA
MR<->AFR,erro12,NA
HET<->hindex,erro13,NA
#End of the model

# Insert de covariance matrix
var<-var(semdata)
cov<-cov(datanew)
cor<-cor(datanew)


#Running SEM model
sem <- sem::sem(semmodel,cor, N=781)
summary(sem,fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", 
                          "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc"))
modIndices(sem)

qgraph(sem,cut = 0.8,gray=TRUE)

#SEM Model B - Working but poor fit indicators

semmodel<- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).
#Latent Variables
HealthOut->IDPT,efa14,NA
HealthOut->IM,efa11,NA
HealthOut->MR,efa12,NA
HealthOut->AFR,efa13,NA
HealthOut->BR,NA,1
HealthOut->HEP,efa16,NA
hindex->GDP,ob1,NA
HET->GDP,ob3,NA
hindex->HealthOut,lat12,NA
HET->HealthOut,lat14,NA
HealthOut<->HealthOut,laterro1,NA
GDP<->GDP,erro1,NA
hindex<->hindex,erro3,NA
noOfDocs<->noOfDocs,erro4,NA
HET<->HET,erro5,NA
#End of the model

# Insert de covariance matrix
var<-var(semdata)
cov<-cov(datanew)
cor<-cor(datanew)


#Running SEM model
sem <- sem::sem(semmodel,cor, N=781)
summary(sem,fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", 
                          "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc"))
modIndices(sem)

qgraph(sem,cut = 0.8,gray=TRUE)

#SEM Model A - Working but poor fit indicators

semmodel<- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).
#Latent Variables
HealthOut->IDPT,efa14,NA
HealthOut->IM,efa11,NA
HealthOut->MR,efa12,NA
HealthOut->AFR,efa13,NA
HealthOut->BR,NA,1
HealthOut->HEP,efa16,NA
RD->hindex,NA,1
RD->HET,efa22,NA
RD->noOfDocs,efa23,NA
RD->GDP,ob1,NA
RD->HealthOut,lat14,NA
HealthOut<->HealthOut,laterro1,NA
GDP<->GDP,erro1,NA
RD<->RD,erro2,NA
hindex<->hindex,erro3,NA
noOfDocs<->noOfDocs,erro4,NA
HET<->HET,erro5,NA
#End of the model

# Insert de covariance matrix
var<-var(semdata)
cov<-cov(datanew)
cor<-cor(datanew)


#Running SEM model
sem <- sem::sem(semmodel,cor, N=781)
summary(sem,fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", 
                          "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc"))
modIndices(sem)

qgraph(sem,cut = 0.8,gray=TRUE)

#SEM Model D - Working but poor fit indicators

semmodel<- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).
#Latent Variables
HealthOut->IDPT,efa14,NA
HealthOut->IM,efa11,NA
HealthOut->MR,efa12,NA
HealthOut->AFR,efa13,NA
HealthOut->BR,NA,1
HealthOut->HEP,efa16,NA
GDP->hindex,ob1,NA
GDP->HET,ob3,NA
hindex->HealthOut,lat12,NA
HET->HealthOut,lat14,NA
HealthOut<->HealthOut,laterro1,NA
GDP<->GDP,erro1,NA
hindex<->hindex,erro3,NA
noOfDocs<->noOfDocs,erro4,NA
HET<->HET,erro5,NA
#End of the model

# Insert de covariance matrix
var<-var(semdata)
cov<-cov(datanew)
cor<-cor(datanew)


#Running SEM model
sem <- sem::sem(semmodel,cor, N=781)
summary(sem,fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", 
                          "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc"))
modIndices(sem)

