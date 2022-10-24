#################################################################################################
#
#   analysis_main_RS.r
#   October 2022
#
#################################################################################################

# clear workspace
rm(list=ls())

# define locations
root <- "~/Managing-Self-Confidence-Replication"
libdir <- paste(root, "code/", sep="/")
datadir <- paste(root, "data/", sep="/")
tabdir <- paste(root, "tables/", sep="/")



# load packages
library(stats)
library(Hmisc)
library(sem)
library(rms)
library(AER)
library(car)
library(gmm)
library(quantreg)
library(plotrix)
library(stargazer)


# load user-created libraries
source(paste(libdir, "lib_multiregtable.r", sep="/"))
source(paste(libdir, "lib_clmclx.r", sep="/"))
source(paste(libdir, "lib_hacktex.r", sep="/"))

### load data and create variables and datasets

# load data from main experiment
data <- read.csv(paste(datadir, "main_session.csv", sep="/"))

# generate implied valuations
data$coarse <- data$talkingbob - data$silentbob
data$precise <- data$precisebob - data$silentbob
data$public_coarse <- data$talkingbob_hof - data$talkingbob
data$public_precise <- data$precisebob_hof - data$precisebob

# misc
data$temptophalf <- as.numeric(data$quiztemprank <= 50)
data$qtmedian <- ave(data$quizscore, as.factor(data$quiztype), FUN=median)
data$tophalf <- as.numeric(data$quizscore >= data$qtmedian)
data$sigsum <- data$signal1 + data$signal2 + data$signal3 + data$signal4



# mean score of others participating in the same quiz
data$meanscore <- ave(data$quizscore, data$quiztype, FUN=mean)
data$quiz.n <- ave(data$quizscore, data$quiztype, FUN=NROW)
data$othermeanscore <- (data$meanscore * data$quiz.n - data$quizscore) / (data$quiz.n - 1)


# reshaped dataset for graphs
graphdata <- data[data$restrict2==1,c("preconfidence1","preconfidence2","preconfidence3","preconfidence4",
                                      "confchange1","confchange2","confchange3","confchange4",
                                      "signal1","signal2","signal3","signal4","male","subjectid")]

graphdata <- reshape(graphdata, varying=c(1:12), 
                     idvar="subjectid", 
                     timevar="round", 
                     direction="long",
                     sep="")

graphdata$confcat <- cut(graphdata$preconfidence, breaks=(0:10)*10, right=FALSE)
graphdata$postconfidence.b <- NA
graphdata$postconfidence.b[graphdata$signal==1] <- 100*(graphdata$preconfidence[graphdata$signal==1] * 0.75) / (graphdata$preconfidence[graphdata$signal==1] * 0.75 + (100 - graphdata$preconfidence[graphdata$signal==1]) * 0.25)
graphdata$postconfidence.b[graphdata$signal==0] <- 100*(graphdata$preconfidence[graphdata$signal==0] * 0.25) / (graphdata$preconfidence[graphdata$signal==0] * 0.25 + (100 - graphdata$preconfidence[graphdata$signal==0]) * 0.75)
graphdata$confchange.b <- graphdata$postconfidence.b - graphdata$preconfidence

# reshaped dataset for regressions
regdata <- data[,c("preconfidence1","preconfidence2","preconfidence3","preconfidence4",
                   "postconfidence1","postconfidence2","postconfidence3","postconfidence4",
                   "signal1","signal2","signal3","signal4","male","subjectid","temptophalf","tophalf",
                   "restrict2","quiztype","quizscore","quiztemprank","year","othermeanscore")]


# restrict to balanced panel
regdata <- regdata[regdata$preconfidence1!=0 & regdata$preconfidence1!=100
                   & regdata$preconfidence2!=0 & regdata$preconfidence2!=100
                   & regdata$preconfidence3!=0 & regdata$preconfidence3!=100
                   & regdata$preconfidence4!=0 & regdata$preconfidence4!=100
                   & regdata$postconfidence4!=0 & regdata$postconfidence4!=100,]

regdata <- reshape(regdata, varying=c(1:12), 
                   idvar="subjectid", 
                   timevar="round", 
                   direction="long",
                   sep="")

regdata$id <- as.factor(regdata$subjectid)

# generate log linear vars
regdata$postconfidence <- regdata$postconfidence / 100
regdata$preconfidence <- regdata$preconfidence / 100
regdata$postodds <- regdata$postconfidence / (1 - regdata$postconfidence)
regdata$preodds <- regdata$preconfidence / (1 - regdata$preconfidence)
regdata$lpostodds <- log(regdata$postodds)
regdata$lpreodds <- log(regdata$preodds)
regdata$llr <- log(0.75 / 0.25)
regdata$llr[regdata$signal==0] <- log(0.25 / 0.75)

# generate interactions
regdata$lpreodds_male <- regdata$lpreodds * regdata$male
regdata$llr_male <- regdata$llr * regdata$male
regdata$lpreodds_tth <- regdata$lpreodds * regdata$temptophalf
regdata$llr_tth <- regdata$llr * regdata$temptophalf

# split LLR variables into +/-
regdata$llr1 <- regdata$llr * as.numeric(regdata$signal==1)
regdata$llr0 <- regdata$llr * as.numeric(regdata$signal==0)
regdata$llr1_male <- regdata$llr1 * regdata$male
regdata$llr0_male <- regdata$llr0 * regdata$male
regdata$llr1_tth <- regdata$llr1 * regdata$temptophalf
regdata$llr0_tth <- regdata$llr0 * regdata$temptophalf

# drop extreme values
regdata <- regdata[is.finite(regdata$lpostodds) & is.finite(regdata$lpreodds),]

### Table 1: belief updating

# create the base regression table
varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                      "lpreodds"            = "$\\delta$",
                      "llr1"                = "$\\beta_H$",
                      "llr0"                = "$\\beta_L$",
                      "male"                = "Male",
                      "lpreodds_male"       = "Log Prior Odds * Male",
                      "llr_male"            = "Log Likelihood Ratio * Male",
                      "lpreodds_tth"        = "Log Prior Odds * Top Half",
                      "llr_tth"             = "Log Likelihood Ratio * Top Half",
                      "temptophalf"         = "Top Half",
                      "Intercept"           = "Intercept",
                      "(Intercept)"         = "Intercept",
                      "d$lpreodds"          = "Log Prior Odds",
                      "d$llr"               = "Log Likelihood Ratio")

# Round-by-round plus pooled regressions
fm.r1 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==1 & regdata$restrict2==1,])
fm.r1$var <- hccm(fm.r1)
fm.r2 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==2 & regdata$restrict2==1,])
fm.r2$var <- hccm(fm.r2)
fm.r3 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==3 & regdata$restrict2==1,])
fm.r3$var <- hccm(fm.r3)
fm.r4 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==4 & regdata$restrict2==1,])
fm.r4$var <- hccm(fm.r4)
fm.ar <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$restrict2==1,])
fm.ar$var <- clx(fm.ar, 1, regdata$subjectid[regdata$restrict2==1])
fm.ar.ur <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata)
fm.ar.ur$var <- clx(fm.ar.ur, 1, regdata$subjectid)


# add tests for responsiveness and symmetry
baseaddrows <- matrix(NA, nrow=3, ncol=7)
baseaddrows[,1] <- c("$\\mathbb{P}(\\beta_H = 1)$","$\\mathbb{P}(\\beta_L = 1)$","$\\mathbb{P}(\\beta_H = \\beta_L)$")
baseaddrows[1,2] <- roundsig(linearHypothesis(fm.r1, c(1,0,0), vcov=fm.r1$var)["Pr(>F)"][[1]][2],3)
baseaddrows[1,3] <- roundsig(linearHypothesis(fm.r2, c(1,0,0), vcov=fm.r2$var)["Pr(>F)"][[1]][2],3)
baseaddrows[1,4] <- roundsig(linearHypothesis(fm.r3, c(1,0,0), vcov=fm.r3$var)["Pr(>F)"][[1]][2],3)
baseaddrows[1,5] <- roundsig(linearHypothesis(fm.r4, c(1,0,0), vcov=fm.r4$var)["Pr(>F)"][[1]][2],3)
baseaddrows[1,6] <- roundsig(linearHypothesis(fm.ar, c(1,0,0), vcov=fm.ar$var)["Pr(>F)"][[1]][2],3)
baseaddrows[1,7] <- roundsig(linearHypothesis(fm.ar.ur, c(1,0,0), vcov=fm.ar.ur$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,2] <- roundsig(linearHypothesis(fm.r1, c(0,1,0), vcov=fm.r1$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,3] <- roundsig(linearHypothesis(fm.r2, c(0,1,0), vcov=fm.r2$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,4] <- roundsig(linearHypothesis(fm.r3, c(0,1,0), vcov=fm.r3$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,5] <- roundsig(linearHypothesis(fm.r4, c(0,1,0), vcov=fm.r4$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,6] <- roundsig(linearHypothesis(fm.ar, c(0,1,0), vcov=fm.ar$var)["Pr(>F)"][[1]][2],3)
baseaddrows[2,7] <- roundsig(linearHypothesis(fm.ar.ur, c(0,1,0), vcov=fm.ar.ur$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,2] <- roundsig(linearHypothesis(fm.r1, c(1,-1,0), vcov=fm.r1$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,3] <- roundsig(linearHypothesis(fm.r2, c(1,-1,0), vcov=fm.r2$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,4] <- roundsig(linearHypothesis(fm.r3, c(1,-1,0), vcov=fm.r3$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,5] <- roundsig(linearHypothesis(fm.r4, c(1,-1,0), vcov=fm.r4$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,6] <- roundsig(linearHypothesis(fm.ar, c(1,-1,0), vcov=fm.ar$var)["Pr(>F)"][[1]][2],3)
baseaddrows[3,7] <- roundsig(linearHypothesis(fm.ar.ur, c(1,-1,0), vcov=fm.ar.ur$var)["Pr(>F)"][[1]][2],3)

vars.base <- c("lpreodds","llr1","llr0")
table.base <- multiregtable(vars.base, varlabels, list(fm.r1,fm.r2,fm.r3,fm.r4,fm.ar,fm.ar.ur), 3, addrows=baseaddrows)

# check first-stage F statistic
fm <- lm(lpreodds ~ as.factor(quiztype), data=regdata[regdata$restrict2==1,])
linearHypothesis(fm, cbind(c(0,0,0,0,0,0,0,0), diag(c(1,1,1,1,1,1,1,1))), rhs=c(0,0,0,0,0,0,0,0), vcov=clx(fm, 1, regdata$subjectid[regdata$restrict2==1]))
linearHypothesis(fm, cbind(c(0,0,0,0,0,0,0,0), diag(c(1,1,1,1,1,1,1,1))), rhs=c(0,0,0,0,0,0,0,0))
fm <- lm(lpreodds ~ othermeanscore, data=regdata[regdata$restrict2==1,])
linearHypothesis(fm, cbind(0,1), vcov=clx(fm, 1, regdata$subjectid[regdata$restrict2==1]))

# define helper functions to estimate GMM models and the associated F-statistics
estgmm <- function(restriction){
  
  lpostodds <- regdata$lpostodds[restriction]
  llr1 <- regdata$llr1[restriction]
  llr0 <- regdata$llr0[restriction]
  lpreodds <- regdata$lpreodds[restriction]
  quiztype <- regdata$quiztype[restriction]
  othermeanscore <- regdata$othermeanscore[restriction]
  return(gmm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, cbind(llr1, llr0, othermeanscore)))
  
}
gmmfstat <- function(restriction){
  firststage <- lm(lpreodds ~ 0 + othermeanscore + llr1 + llr0, data=regdata[restriction,])
  firststage$var <- clx(firststage, 1, regdata$subjectid[restriction])
  return((firststage$coefficient["othermeanscore"] / sqrt(firststage$var["othermeanscore","othermeanscore"]))^2)
}

# estimate GMM models
gmmmodels <- list()
gmmfstats <- list()
for (r in 1:4){
  gmmmodels[[r]] <- estgmm(regdata$restrict2==1 & regdata$round==r)
  gmmfstats[[r]] <- gmmfstat(regdata$restrict2==1 & regdata$round==r)
}
gmmmodels[[5]] <- estgmm(regdata$restrict2==1)
gmmfstats[[5]] <- gmmfstat(regdata$restrict2==1)
gmmmodels[[6]] <- estgmm(TRUE)
gmmfstats[[6]] <- gmmfstat(TRUE)

# add tests for symmetry
gmmaddrows <- matrix(NA, nrow=4, ncol=7)
gmmaddrows[,1] <- c("$\\mathbb{P}(\\beta_H = 1)$","$\\mathbb{P}(\\beta_L = 1)$","$\\mathbb{P}(\\beta_H = \\beta_L)$","First Stage $F$-statistic")
gmmaddrows[1,2] <- roundsig(linearHypothesis(gmmmodels[[1]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[1,3] <- roundsig(linearHypothesis(gmmmodels[[2]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[1,4] <- roundsig(linearHypothesis(gmmmodels[[3]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[1,5] <- roundsig(linearHypothesis(gmmmodels[[4]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[1,6] <- roundsig(linearHypothesis(gmmmodels[[5]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[1,7] <- roundsig(linearHypothesis(gmmmodels[[6]], c(1,0,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,2] <- roundsig(linearHypothesis(gmmmodels[[1]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,3] <- roundsig(linearHypothesis(gmmmodels[[2]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,4] <- roundsig(linearHypothesis(gmmmodels[[3]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,5] <- roundsig(linearHypothesis(gmmmodels[[4]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,6] <- roundsig(linearHypothesis(gmmmodels[[5]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[2,7] <- roundsig(linearHypothesis(gmmmodels[[6]], c(0,1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,2] <- roundsig(linearHypothesis(gmmmodels[[1]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,3] <- roundsig(linearHypothesis(gmmmodels[[2]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,4] <- roundsig(linearHypothesis(gmmmodels[[3]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,5] <- roundsig(linearHypothesis(gmmmodels[[4]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,6] <- roundsig(linearHypothesis(gmmmodels[[5]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[3,7] <- roundsig(linearHypothesis(gmmmodels[[6]], c(1,-1,0))["Pr(>Chisq)"][[1]][2],3)
gmmaddrows[4,2:7] <- roundsig(as.numeric(gmmfstats), 2)

vars.gmm <- c("lpreodds","llr1","llr0")
table.gmm <- multiregtable(vars.gmm, varlabels, gmmmodels, 3, addrows=gmmaddrows)

# output combined table
result <- hacktex(rbind(table.base, table.gmm), 
                  file=paste(tabdir, "base_combo.tex", sep="/"),
                  label="tab:base_combo",
                  table.env=FALSE,
                  caption.loc="top",
                  center="none",
                  rowlabel="",
                  rowname=rep("",23),
                  rgroup=c("Panel A: OLS","Panel B: IV"),
                  n.rgroup=c(11,12),
                  colheads=c("Regressor","Round 1","Round 2","Round 3","Round 4","All Rounds","Unrestricted"),
                  collabel.just=c("l","c","c","c","c","c","c"))

####----GENDER ANALYSIS----####
# Heterogeneity in updating by gender
# OLS
fm.g.r2 <- lm(lpostodds ~ 0 + llr1 + llr1_male + llr0 + llr0_male + lpreodds + lpreodds_male, data=regdata[regdata$restrict2==1,])
fm.g.r2$var <- clx(fm.g.r2, 1, regdata$subjectid[regdata$restrict2==1])

stargazer(fm.g.r2, type = "text")
# OLS with delta = 1 imposed
fm.g.r2 <- lm(lpostodds - lpreodds ~ 0 + llr1 + llr1_male + llr0 + llr0_male, data=regdata[regdata$restrict2==1,])
fm.g.r2$var <- clx(fm.g.r2, 1, regdata$subjectid[regdata$restrict2==1])

# create variables and interacted instruments
lpostodds <- regdata$lpostodds[regdata$restrict2==1]
llr1 <- regdata$llr1[regdata$restrict2==1]
llr0 <- regdata$llr0[regdata$restrict2==1]
lpreodds <- regdata$lpreodds[regdata$restrict2==1]
oms <- regdata$othermeanscore[regdata$restrict2==1]
lpreodds_male <- regdata$lpreodds_male[regdata$restrict2==1]
llr1_male <- regdata$llr1_male[regdata$restrict2==1]
llr0_male <- regdata$llr0_male[regdata$restrict2==1]
regdata$oms_male <- regdata$othermeanscore * regdata$male
oms_male <- regdata$oms_male[regdata$restrict2==1]


# create variables and interacted instruments
lpostodds <- regdata$lpostodds[regdata$restrict2==1]
llr1 <- regdata$llr1[regdata$restrict2==1]
llr0 <- regdata$llr0[regdata$restrict2==1]
lpreodds <- regdata$lpreodds[regdata$restrict2==1]
oms <- regdata$othermeanscore[regdata$restrict2==1]
lpreodds_male <- regdata$lpreodds_male[regdata$restrict2==1]
llr1_male <- regdata$llr1_male[regdata$restrict2==1]
llr0_male <- regdata$llr0_male[regdata$restrict2==1]
regdata$oms_male <- regdata$othermeanscore * regdata$male
oms_male <- regdata$oms_male[regdata$restrict2==1]

# estimation
gmm.g.r2 <- gmm(lpostodds ~ 0 + llr1 + llr0 + llr1_male + llr0_male + lpreodds + lpreodds_male, cbind(llr1,llr0,llr1_male,llr0_male,oms,oms_male))

# create the differential regression table
varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                      "llr1"                = "$\\beta_H$",
                      "llr0"                = "$\\beta_L$",
                      "llr1_male"           = "$\\beta_H^{Male}$",
                      "llr0_male"           = "$\\beta_L^{Male}$",
                      "lpreodds"            = "$\\delta$",
                      "lpreodds_male"       = "$\\delta^{Male}$")

vars.gender <- c("llr1","llr0","llr1_male","llr0_male")

table.gender <- multiregtable(vars.gender, varlabels, list(fm.g.r2), 3)

# Table 2: Gender heterogeneity in beleif updating
result <- hacktex(table.gender, 
                  file=paste(tabdir, "differential_gender.tex", sep="/"),
                  label="tab:differential_gender",
                  tabwidth="3in",
                  table.env=FALSE,
                  caption.loc="top",
                  rowname =NULL,
                  center="none",
                  colheads=c("Regressor","OLS"),
                  collabel.just=c("l","c"))

