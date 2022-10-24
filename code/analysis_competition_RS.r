#################################################################################################
#
#   analysis_competition_RS.r
#   October 2022
#
#################################################################################################

# clear workspace
rm(list=ls())

# define locations
root <- "~/Managing-Self-Confidence-Replication"
figdir <- paste(root, "figures/", sep="/")
tabdir <- paste(root, "tables/", sep="/")
datadir <- paste(root, "data/", sep="/")
libdir <- paste(root, "code/", sep="/")

# load packages
library(erer)
library(gmm)
library(Hmisc)
library(gplots)
library(car)

# load user-created libraries
source(paste(libdir, "lib_multiregtable.r", sep="/"))
source(paste(libdir, "lib_clmclx.r", sep="/"))
source(paste(libdir, "lib_hacktex.r", sep="/"))

# create a utility function
rfind <- function(x) seq(along=x)[x != 0]

### load data and create variables and datasets

# load base data
basedata <- read.csv(paste(datadir, "main_session.csv", sep="/"))

# mean score of others participating in the same quiz
basedata$meanscore <- ave(basedata$quizscore, basedata$quiztype, FUN=mean)
basedata$quiz.n <- ave(basedata$quizscore, basedata$quiztype, FUN=NROW)
basedata$othermeanscore <- (basedata$meanscore * basedata$quiz.n - basedata$quizscore) / (basedata$quiz.n - 1)

# calculate conservatism measures
basedata$logit0 <- log(basedata$top50afterquiz/(100-basedata$top50afterquiz))
basedata$logit1 <- log(basedata$top50afters1/(100-basedata$top50afters1))
basedata$logit2 <- log(basedata$top50afters2/(100-basedata$top50afters2))
basedata$logit3 <- log(basedata$top50afters3/(100-basedata$top50afters3))
basedata$logit4 <- log(basedata$top50afters4/(100-basedata$top50afters4))
basedata$dev1 <- (basedata$logit1-basedata$logit0)/log(3)
basedata$dev2 <- (basedata$logit2-basedata$logit1)/log(3)
basedata$dev3 <- (basedata$logit3-basedata$logit2)/log(3)
basedata$dev4 <- (basedata$logit4-basedata$logit3)/log(3)
basedata$responsivenessdefined <- abs(basedata$logit0)!=Inf & abs(basedata$logit1)!=Inf & abs(basedata$logit2)!=Inf & abs(basedata$logit3)!=Inf & abs(basedata$logit4)!=Inf
basedata$temptophalf <- as.numeric(basedata$quiztemprank <= 50)
basedata$signalsumnormalized <- basedata$signalsum-1-2*basedata$temptophalf
basedata$goodnews<-basedata$signalsum-basedata$signalsum
ind<-rfind(basedata$temptophalf==0 & basedata$signalsum>1)
basedata$goodnews[ind]=1
ind<-rfind(basedata$temptophalf==1 & basedata$signalsum>2)
basedata$goodnews[ind]=1	

# load competition data
compdata <- read.csv(paste(datadir, "competition_session.csv", sep="/"))
names(compdata) <- tolower(names(compdata))
compdata$task2estimatedwin=c()

# calculate probability of winning tournament
for(i in 1:2){
  ind<-rfind(compdata$chargame==i)
  scores=compdata$task1score[ind]
  for(j in 1:length(ind)){
    ind_lower<-rfind(scores<scores[j])
    ind_same<-rfind(scores==scores[j])
    A=length(ind_lower)
    B=length(ind_same)-1
    T=length(ind)-1
    prob=(A*(A-1)*(A-2)/6+A*(A-1)/2*B*0.5+A*B*(B-1)/2*1/3+B*(B-1)*(B-2)*0.25)/(T*(T-1)*(T-2)/6)
    compdata$task2estimatedwin[ind[j]]=prob*100
  }
}

compdata$task2estimatedwinlogit=log(compdata$task2estimatedwin/(100-compdata$task2estimatedwin))
ind<-rfind(compdata$task2estimatedwinlogit==Inf)
compdata$task2estimatedwinlogit[ind]=log(99);
ind<-rfind(compdata$task2estimatedwinlogit==-Inf)
compdata$task2estimatedwinlogit[ind]=-log(99);
compdata$task2preconfidencelogit=log(compdata$task2preconfidence/(100-compdata$task2preconfidence))
ind<-rfind(compdata$task2preconfidencelogit==Inf)
compdata$task2preconfidencelogit[ind]=log(99);
ind<-rfind(compdata$task2preconfidencelogit==-Inf)
compdata$task2preconfidencelogit[ind]=-log(99);

# merge data from base experiment and competition experiment
data <- merge(basedata, compdata, by="subjectid")

# define top half
data$temptophalf <- as.numeric(data$quiztemprank <= 50)

# restrict to observations with reasonable updating behavior
data <- data[data$restrict2==1,]
basedata <- basedata[basedata$restrict2==1,]

# subset for gender-specific analysis
mdata <- data[data$male==1,]
fdata <- data[data$male==0,]

# restrict to observations where responsiveness is defined	
respbasedata <- basedata[basedata$responsivenessdefined,]

# calculate average up/down responsiveness
dev=list()
devrank=list()
dev[[1]]<-respbasedata$dev1
dev[[2]]<-respbasedata$dev2
dev[[3]]<-respbasedata$dev3
dev[[4]]<-respbasedata$dev4
for(i in 1:4){
  devrank[[i]]<-dev[[i]]
  #up responsiveness
  ind<-rfind(dev[[i]]>0)
  devrank[[i]][ind]=rank(dev[[i]][ind])/length(ind)
  #down responsiveness
  ind<-rfind(dev[[i]]<=0)
  devrank[[i]][ind]=rank(abs(dev[[i]][ind]))/length(ind)
}
respbasedata$responsiveness <-(devrank[[1]]+devrank[[2]]+devrank[[3]]+devrank[[4]])/4	
respbasedata$conservative <-as.numeric(respbasedata$responsiveness<=0.5)
respdata <- merge(respbasedata, compdata, by="subjectid")
respdata$temptophalf <- as.numeric(respdata$quiztemprank <= 50)

### Universal variable labels ##########################################################################

# variable labels
varlabels <- pairlist("task2preconfidence" = "Confidence (Experiment 2)",
                      "task2preconfidencelogit" = "Confidence (Experiment 2)",
                      "signalsum"          = "Feedback (Experiment 1)",
                      "goodnews"          = "Feedback (Experiment 1)",
                      "temptophalf"        = "Ability (Experiment 1)",
                      "othermeanscore"     = "Mean score of others (Experiment 1)",
                      "task2estimatedwin"  = "Predicted winning prob.",
                      "task2estimatedwinlogit"  = "Predicted winning prob.",
                      "conservative" = "Conservative",
                      "task2preconfidence_lc" = "Confidence (Experiment 2) * LC",
                      "task2preconfidence_male" = "Confidence (Experiment 2) * Male")


###---GENDER ANALYSIS---###
### Table 3: How confidence affects competition for male vs female      

# male regressions
attach(mdata)  
fm.firststage.male <- lm(task2preconfidence ~ signalsum + temptophalf)
fm.firststage.male$var <- hccm(fm.firststage.male)
fm.redform.male <- lm(task2tournament ~ signalsum + temptophalf)
fm.redform.male$var <- hccm(fm.redform.male)
fm.gmm.male <- gmm(task2tournament ~ task2preconfidence + temptophalf, cbind(signalsum, temptophalf))
detach(mdata)

# female regressions
attach(fdata)
fm.firststage.fem <- lm(task2preconfidence ~ signalsum + temptophalf)
fm.firststage.fem$var <- hccm(fm.firststage.fem)
fm.redform.fem <- lm(task2tournament ~ signalsum + temptophalf)
fm.redform.fem$var <- hccm(fm.redform.fem)
fm.gmm.fem <- gmm(task2tournament ~ task2preconfidence + temptophalf, cbind(signalsum, temptophalf))
detach(fdata)

# output
vars.gmm.gender <- c("task2preconfidence","signalsum","temptophalf")
table.gmm.gender <- multiregtable(vars.gmm.gender, varlabels, list(fm.firststage.fem, 
                                                                               fm.firststage.male, 
                                                                               fm.redform.fem,
                                                                               fm.redform.male,
                                                                               fm.gmm.fem,
                                                                               fm.gmm.male), 3)
result <- hacktex(table.gmm.gender, 
                  file=paste(tabdir, "application_gender.tex", sep="/"),
                  table.env=FALSE,
                  caption.loc="top",
                  center="none",
                  rowname =NULL,
                  cgroup=c("","First Stage","Reduced Form","IV"),
                  n.cgroup=c(1,2,2,2),
                  colheads=c("","Female","Male","Female","Male","Female","Male"),
                  collabel.just=c("l","c","c","c","c","c","c"))

### Figure 1: Men or women have less accurate beliefs?
postscript(file=paste(figdir, "task2estimatedwin_preconfidence_gender_exp2.eps", sep=""), onefile=FALSE, horizontal=FALSE, width=7, height=7)
par(mar=c(5.1, 4.1, 1.1, 2.1))
fm.male <- lm(task2preconfidence ~ task2estimatedwin, data=respdata[respdata$male==1,])
fm.fem <- lm(task2preconfidence ~ task2estimatedwin, data=respdata[respdata$male==0,])
fv.male <- predict(fm.male, newdata=data.frame(task2estimatedwin=seq(min(respdata$task2estimatedwin), max(respdata$task2estimatedwin),1)))/100
fv.fem <- predict(fm.fem, newdata=data.frame(task2estimatedwin=seq(min(respdata$task2estimatedwin), max(respdata$task2estimatedwin),1)))/100
plot(seq(min(respdata$task2estimatedwin), max(respdata$task2estimatedwin),1), fv.male, type="l", lwd=2, ylim=c(0,1), col=grey(0.8), xlab="Predicted winning probability", ylab="Confidence") 
lines(seq(min(respdata$task2estimatedwin), max(respdata$task2estimatedwin),1), fv.fem, lwd=2, col=grey(0.3))
points(respdata$task2estimatedwin[respdata$male==1], respdata$task2preconfidence[respdata$male==1]/100, pch=16, col=grey(0.8))
points(respdata$task2estimatedwin[respdata$male==0], respdata$task2preconfidence[respdata$male==0]/100, pch=16, col=grey(0.3))
dev.off()



