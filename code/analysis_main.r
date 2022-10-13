#################################################################################################
#
#   analysis_main.r
#   July 2021
#
#################################################################################################

    # clear workspace
    rm(list=ls())
    
    # define locations
    root <- "D:/Dropbox/research/MSC/replication_package"
    figdir <- paste(root, "figures/", sep="/")
    tabdir <- paste(root, "tables/", sep="/")
    datadir <- paste(root, "data/", sep="/")
    libdir <- paste(root, "code/", sep="/")    
    
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

### make dataset for comparison with the makeup (MU) experiment

    muregdata <- data[,c("preconfidence1","preconfidence2","preconfidence3","preconfidence4",
                         "postconfidence1","postconfidence2","postconfidence3","postconfidence4",
                         "signal1","signal2","signal3","signal4",
                         "mupreconfidence1","mupreconfidence2","mupreconfidence3","mupreconfidence4",
                         "mupostconfidence1","mupostconfidence2","mupostconfidence3","mupostconfidence4",
                         "musignal1","musignal2","musignal3","musignal4","male","subjectid","temptophalf","tophalf",
                         "restrict2","murestrict2","restrict1","murestrict1")]
    
    # restrict to balanced panel
    muregdata <- muregdata[muregdata$preconfidence1!=0 & muregdata$preconfidence1!=100
                            & muregdata$preconfidence2!=0 & muregdata$preconfidence2!=100
                            & muregdata$preconfidence3!=0 & muregdata$preconfidence3!=100
                            & muregdata$preconfidence4!=0 & muregdata$preconfidence4!=100
                            & muregdata$postconfidence4!=0 & muregdata$postconfidence4!=100
                            & muregdata$mupreconfidence1!=0 & muregdata$mupreconfidence1!=100
                            & muregdata$mupreconfidence2!=0 & muregdata$mupreconfidence2!=100
                            & muregdata$mupreconfidence3!=0 & muregdata$mupreconfidence3!=100
                            & muregdata$mupreconfidence4!=0 & muregdata$mupreconfidence4!=100
                            & muregdata$mupostconfidence4!=0 & muregdata$mupostconfidence4!=100
                            & !is.na(muregdata$mupreconfidence1),]
    muregdata <- reshape(muregdata, varying=c(1:24), 
                                    idvar="subjectid", 
                                    timevar="round", 
                                    direction="long",
                                    sep="")
                                    
    # generate log linear vars
    muregdata$postconfidence <- muregdata$postconfidence / 100
    muregdata$preconfidence <- muregdata$preconfidence / 100
    muregdata$postodds <- muregdata$postconfidence / (1 - muregdata$postconfidence)
    muregdata$preodds <- muregdata$preconfidence / (1 - muregdata$preconfidence)
    muregdata$lpostodds <- log(muregdata$postodds)
    muregdata$lpreodds <- log(muregdata$preodds)
    muregdata$llr <- log(0.75 / 0.25)
    muregdata$llr[muregdata$signal==0] <- log(0.25 / 0.75)
    muregdata$llr1 <- muregdata$llr * as.numeric(muregdata$signal==1)
    muregdata$llr0 <- muregdata$llr * as.numeric(muregdata$signal==0)
    
    muregdata$mupostconfidence <- muregdata$mupostconfidence / 100
    muregdata$mupreconfidence <- muregdata$mupreconfidence / 100
    muregdata$mupostodds <- muregdata$mupostconfidence / (1 - muregdata$mupostconfidence)
    muregdata$mupreodds <- muregdata$mupreconfidence / (1 - muregdata$mupreconfidence)
    muregdata$mulpostodds <- log(muregdata$mupostodds)
    muregdata$mulpreodds <- log(muregdata$mupreodds)
    muregdata$mullr <- log(0.75 / 0.25)
    muregdata$mullr[muregdata$musignal==0] <- log(0.25 / 0.75)
    muregdata$mullr1 <- muregdata$mullr * as.numeric(muregdata$musignal==1)
    muregdata$mullr0 <- muregdata$mullr * as.numeric(muregdata$musignal==0)

    muregdata$loddschange <- muregdata$lpostodds - muregdata$lpreodds
    muregdata$muloddschange <- muregdata$mulpostodds - muregdata$mulpreodds


### Footnote 16: expected earnings
    
    # earnings from quiz: 25 cents per net correct answer, with a min of zero
    data$quiz_earnings <- pmax(data$quizscore * 0.25, 0)
    
    # (expected) earnings from belief updating
    exp_earnings_belief <- function(belief, winner){
      return(belief/100 * winner + (1 - belief/100)^2 / 2)
    }
    # for guys for whom we did not elicit top25...
    data$belief_earnings <- 3*(exp_earnings_belief(data$top50,data$temptophalf)
                            + exp_earnings_belief(data$top50afterquiz,data$temptophalf)
                            + exp_earnings_belief(data$top50afters1,data$temptophalf)
                            + exp_earnings_belief(data$top50afters2,data$temptophalf)
                            + exp_earnings_belief(data$top50afters3,data$temptophalf)
                            + exp_earnings_belief(data$top50afters4,data$temptophalf)) / 6
    # for guys for whom we *did* elicit top25...
    p75 <- function(x){return(quantile(x, probs=c(0.75)))}
    data$qt25 <- ave(data$quizscore, as.factor(data$quiztype), FUN=p75)
    data$intop25 <- as.numeric(data$quizscore >= data$qt25)
    data$belief_earnings[!is.na(data$top25)] <- 3*(exp_earnings_belief(data$top50[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top25[!is.na(data$top25)],data$intop25[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top50afterquiz[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top50afters1[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top50afters2[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top50afters3[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])
                                                 + exp_earnings_belief(data$top50afters4[!is.na(data$top25)],data$tophalf[!is.na(data$top25)])) / 7
    
    # (expected) earnings from information demand questions
    exp_earnings_demand <- function(bid){
      return(2*bid/100 - (bid/100)^2/2)
    }
    data$demand_earnings <- (exp_earnings_demand(data$silentbob)
                             + exp_earnings_demand(data$talkingbob)
                             + exp_earnings_demand(data$precisebob)
                             + exp_earnings_demand(data$talkingbob_hof)
                             + exp_earnings_demand(data$precisebob_hof)) / 5
    # we're missing demand variables for 14 people
    data[is.na(data$demand_earnings), c("silentbob","talkingbob","precisebob","talkingbob_hof","precisebob_hof")]
    
    # total earnings
    data$total_expected_earnings <- data$quiz_earnings + data$belief_earnings + data$demand_earnings
    summary(data$total_expected_earnings)
    sd(data$total_expected_earnings, na.rm=TRUE)
      
    
### Figure 2: belief distributions

    # post-quiz and post-signal4
    postscript(file=paste(figdir, "belief_cdfs04.eps", sep=""), width=6, height=4.5, horizontal=FALSE)
    par(lwd=2, mar=c(4.5,2,1,2))
    plot(ecdf(data$top50afterquiz[data$restrict2==1]), verticals=TRUE, pch=NA_integer_,
            xlim=c(0,100), xlab="Belief", main="", col=grey(0.8), lwd=2)
    lines(ecdf(data$top50afters4[data$restrict2==1]), verticals=TRUE, pch=NA_integer_, xlim=c(0,100), col=grey(0.3), lwd=2)
    legend(x=0, y=1, legend=c("Post Quiz","Post Signal 4"), lty=c("solid","solid"), col=grey(c(0.8,0.3)), lwd=c(2,2))
    dev.off()
    
### Section 3.1: subject pool

    # sample summary stats
    table(data$male) / nrow(data)
    table(data$year) / nrow(data)
    
### Section 3.2 / Table S-4: summary stats on quiz scores

    # function to calculate summary stats
    rowfun <- function(sd){
    
        qr <- sd$quizright
        qw <- sd$quizwrong
        return(c(nrow(sd),roundsig(c(mean(qr), sd(qr), mean(qw), sd(qw), mean(qr - qw), sd(qr - qw)), 1)))
    
    }
    
    # create the table
    sstab <- matrix("", nrow=16, ncol=7)
    sstab[2,] <- rowfun(data[data$restrict2==1,])
    sstab[3,] <- rowfun(data)
    qts <- unique(data$quiztype)
    for (i in 1:9){
        sstab[4+i,] <- rowfun(data[data$quiztype==qts[i] & data$restrict2==1,])
    }
    sstab[15,] <- rowfun(data[data$male==1 & data$restrict2==1,])
    sstab[16,] <- rowfun(data[data$male==0 & data$restrict2==1,])
    sstab <- data.frame(c("\\textbf{Overall}","\\hspace{1em}Restricted Sample","\\hspace{1em}Full Sample","\\textbf{By Quiz Type}","\\hspace{1em}1","\\hspace{1em}2","\\hspace{1em}3","\\hspace{1em}4","\\hspace{1em}5","\\hspace{1em}6","\\hspace{1em}7","\\hspace{1em}8","\\hspace{1em}9","\\textbf{By Gender}","\\hspace{1em}Male","\\hspace{1em}Female"), sstab)
    
    # output
    result <- hacktex(sstab, 
                    file=paste(tabdir, "quiz_scores.tex", sep=""),
                    label="tab:quiz_scores",
                    table.env=FALSE,
                    caption.loc="top",
                    rowname=NULL,
                    center="none",
                    align="lrrrrrrr",
                    colheads=c(" ","$N$","Mean","SD","Mean","SD","Mean","SD"),
                    collabel.just=c("l","c","c","c","c","c","c","c"),
                    cgroup=c("","Correct","Incorrect","Score"),
                    n.cgroup=c(2,2,2,2))

### Figures 3 & 4

    # find cell means
    poschanges <- tapply(graphdata$confchange[graphdata$signal==1], graphdata$confcat[graphdata$signal==1], FUN=mean)
    negchanges <- tapply(graphdata$confchange[graphdata$signal==0], graphdata$confcat[graphdata$signal==0], FUN=mean)
    poschanges.b <- tapply(graphdata$confchange.b[graphdata$signal==1], graphdata$confcat[graphdata$signal==1], FUN=mean)
    negchanges.b <- tapply(graphdata$confchange.b[graphdata$signal==0], graphdata$confcat[graphdata$signal==0], FUN=mean)
    poschanges.m <- tapply(graphdata$confchange[graphdata$signal==1 & graphdata$male==1], graphdata$confcat[graphdata$signal==1 & graphdata$male==1], FUN=mean)
    negchanges.m <- tapply(graphdata$confchange[graphdata$signal==0 & graphdata$male==1], graphdata$confcat[graphdata$signal==0 & graphdata$male==1], FUN=mean)
    poschanges.w <- tapply(graphdata$confchange[graphdata$signal==1 & graphdata$male==0], graphdata$confcat[graphdata$signal==1 & graphdata$male==0], FUN=mean)
    negchanges.w <- tapply(graphdata$confchange[graphdata$signal==0 & graphdata$male==0], graphdata$confcat[graphdata$signal==0 & graphdata$male==0], FUN=mean)

    # find confidence intervals
    meansd <- function(subset){
        sd <- tapply(graphdata$confchange[subset], graphdata$confcat[subset], FUN=sd)
        n <- tapply(graphdata$confchange[subset], graphdata$confcat[subset], FUN=NROW)
        return(sd / sqrt(n))
    }
    meansd.pos <- meansd(graphdata$signal==1)
    meansd.neg <- meansd(graphdata$signal==0)
    meansd.pos.m <- meansd(graphdata$signal==1 & graphdata$male==1)
    meansd.neg.m <- meansd(graphdata$signal==0 & graphdata$male==1)
    meansd.pos.w <- meansd(graphdata$signal==1 & graphdata$male==0)
    meansd.neg.w <- meansd(graphdata$signal==0 & graphdata$male==0)

    # conservatism graph
    postscript(file=paste(figdir,"conservatism_bars.eps",sep=""), width=6, height=4.5, horizontal=F)
    par(mar=c(4.5,2,1,2))
    barplot(rbind(poschanges.b, poschanges), beside=TRUE, ylim=c(-30,30), 
            legend.text=c("Bayes","Actual"),
            col=grey(c(0.8,0.3)),
            las=3,
            names.arg=c("0-9%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-100%"))
    barplot(rbind(negchanges.b, negchanges), beside=TRUE, add=TRUE, col=grey(c(0.8,0.3)), axisnames=FALSE)
    plotCI(x = (1:10)*3 - 0.5, y=poschanges, uiw = 2 * meansd.pos, add=TRUE, pch=NA)
    plotCI(x = (1:10)*3 - 0.5, y=negchanges, uiw = 2 * meansd.neg, add=TRUE, pch=NA)
    dev.off()
    
    # conservatism graph
    pdf(file=paste(figdir,"conservatism_bars.pdf",sep=""), width=6, height=4.5)
    par(mar=c(4.5,2,1,2))
    barplot(rbind(poschanges.b, poschanges), beside=TRUE, ylim=c(-30,30), 
            legend.text=c("Bayes","Actual"),
            col=grey(c(0.8,0.3)),
            las=3,
            names.arg=c("0-9%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-100%"))
    barplot(rbind(negchanges.b, negchanges), beside=TRUE, add=TRUE, col=grey(c(0.8,0.3)), axisnames=FALSE)
    plotCI(x = (1:10)*3 - 0.5, y=poschanges, uiw = 2 * meansd.pos, add=TRUE, pch=NA)
    plotCI(x = (1:10)*3 - 0.5, y=negchanges, uiw = 2 * meansd.neg, add=TRUE, pch=NA)
    dev.off()
    
    # asymmetry graph
    postscript(file=paste(figdir,"asymmetry_bars.eps",sep=""), width=6, height=4.5, horizontal=F)
    par(mar=c(4.5,2,1,2))
    barplot(rbind(poschanges, -rev(negchanges)), beside=TRUE, 
            legend.text=c("Positive","Negative"),
            col=grey(c(0.8,0.3)),
            las=3,
            ylim=c(0,max(poschanges + 2 * meansd.pos)),
            names.arg=c("0-9%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-100%"))
    plotCI(x = (1:10)*3 - 1.5, y=poschanges, uiw = 2 * meansd.pos, add=TRUE, pch=NA)
    plotCI(x = (1:10)*3 - 0.5, y=-rev(negchanges), uiw = 2 * meansd.neg, add=TRUE, pch=NA)
    dev.off()
    
    # asymmetry graph
    pdf(file=paste(figdir,"asymmetry_bars.pdf",sep=""), width=6, height=4.5)
    par(mar=c(4.5,2,1,2))
    barplot(rbind(poschanges, -rev(negchanges)), beside=TRUE, 
            legend.text=c("Positive","Negative"),
            col=grey(c(0.8,0.3)),
            las=3,
            ylim=c(0,max(poschanges + 2 * meansd.pos)),
            names.arg=c("0-9%","10-19%","20-29%","30-39%","40-49%","50-59%","60-69%","70-79%","80-89%","90-100%"))
    plotCI(x = (1:10)*3 - 1.5, y=poschanges, uiw = 2 * meansd.pos, add=TRUE, pch=NA)
    plotCI(x = (1:10)*3 - 0.5, y=-rev(negchanges), uiw = 2 * meansd.neg, add=TRUE, pch=NA)
    dev.off()
    
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

    # test for final-round weight on prior of 1
    2*(1 - pnorm(abs((1 - fm.r4$coefficients[3]) / sqrt(fm.r4$var[3,3]))))
    
    # comparison to naive regression of LPO on Bayesian LPO
    regdata$lbayesodds <- regdata$llr + regdata$lpreodds
    summary(lm(lpostodds ~ lbayesodds, data=regdata))
    
    # simple 2+/2- comparison
    data$netzero <- as.numeric(data$signal1+data$signal2+data$signal3+data$signal4==2)
    NROW(data$postconfidence4[data$netzero==1])  
    t.test(data$postconfidence4[data$netzero==1], data$preconfidence1[data$netzero==1], paired=TRUE)
    
    # calculation of relative asymmetry
    2*(0.037) / (0.336 + 0.037)

### Table S-6: updating regressions restricted to observations with non-zero changes in beliefs

    # toggle sample: exclude non-updates
    regdata$restrict3 <- regdata$restrict2 * as.numeric(regdata$lpostodds!=regdata$lpreodds)
    
    # Round-by-round plus pooled regressions
    fm.r1 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==1 & regdata$restrict3==1,])
    fm.r1$var <- hccm(fm.r1)
    fm.r2 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==2 & regdata$restrict3==1,])
    fm.r2$var <- hccm(fm.r2)
    fm.r3 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==3 & regdata$restrict3==1,])
    fm.r3$var <- hccm(fm.r3)
    fm.r4 <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$round==4 & regdata$restrict3==1,])
    fm.r4$var <- hccm(fm.r4)
    fm.ar <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds, data=regdata[regdata$restrict3==1,])
    fm.ar$var <- clx(fm.ar, 1, regdata$subjectid[regdata$restrict3==1])
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
    fm <- lm(lpreodds ~ as.factor(quiztype), data=regdata[regdata$restrict3==1,])
    linearHypothesis(fm, cbind(c(0,0,0,0,0,0,0,0), diag(c(1,1,1,1,1,1,1,1))), rhs=c(0,0,0,0,0,0,0,0), vcov=clx(fm, 1, regdata$subjectid[regdata$restrict3==1]))
    linearHypothesis(fm, cbind(c(0,0,0,0,0,0,0,0), diag(c(1,1,1,1,1,1,1,1))), rhs=c(0,0,0,0,0,0,0,0))
    fm <- lm(lpreodds ~ othermeanscore, data=regdata[regdata$restrict3==1,])
    linearHypothesis(fm, cbind(0,1), vcov=clx(fm, 1, regdata$subjectid[regdata$restrict3==1]))
    
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
      gmmmodels[[r]] <- estgmm(regdata$restrict3==1 & regdata$round==r)
      gmmfstats[[r]] <- gmmfstat(regdata$restrict3==1 & regdata$round==r)
    }
    gmmmodels[[5]] <- estgmm(regdata$restrict3==1)
    gmmfstats[[5]] <- gmmfstat(regdata$restrict3==1)
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
                      file=paste(tabdir, "base_combo_movers.tex", sep="/"),
                      label="tab:base_combo_movers",
                      table.env=FALSE,
                      caption.loc="top",
                      center="none",
                      rowlabel="",
                      rowname=rep("",23),
                      rgroup=c("Panel A: OLS","Panel B: IV"),
                      n.rgroup=c(11,12),
                      colheads=c("Regressor","Round 1","Round 2","Round 3","Round 4","All Rounds","Unrestricted"),
                      collabel.just=c("l","c","c","c","c","c","c"))
        
### Section 4.1: tests for stability of coefficients across rounds
    
    # test equality of coefficients across rounds
    regdata$llr0_r1 <- regdata$llr0 * as.numeric(regdata$round==1)
    regdata$llr0_r2 <- regdata$llr0 * as.numeric(regdata$round==2)
    regdata$llr0_r3 <- regdata$llr0 * as.numeric(regdata$round==3)
    regdata$llr0_r4 <- regdata$llr0 * as.numeric(regdata$round==4)
    regdata$llr1_r1 <- regdata$llr1 * as.numeric(regdata$round==1)
    regdata$llr1_r2 <- regdata$llr1 * as.numeric(regdata$round==2)
    regdata$llr1_r3 <- regdata$llr1 * as.numeric(regdata$round==3)
    regdata$llr1_r4 <- regdata$llr1 * as.numeric(regdata$round==4)
    regdata$lpreodds_r1 <- regdata$lpreodds * as.numeric(regdata$round==1)
    regdata$lpreodds_r2 <- regdata$lpreodds * as.numeric(regdata$round==2)
    regdata$lpreodds_r3 <- regdata$lpreodds * as.numeric(regdata$round==3)
    regdata$lpreodds_r4 <- regdata$lpreodds * as.numeric(regdata$round==4)
    regdata$othermeanscore_r1 <- regdata$othermeanscore * as.numeric(regdata$round==1)
    regdata$othermeanscore_r2 <- regdata$othermeanscore * as.numeric(regdata$round==2)
    regdata$othermeanscore_r3 <- regdata$othermeanscore * as.numeric(regdata$round==3)
    regdata$othermeanscore_r4 <- regdata$othermeanscore * as.numeric(regdata$round==4)

    ## OLS
    fm.ar.ct <- lm(lpostodds ~ 0 + llr1_r1 + llr1_r2 + llr1_r3 + llr1_r4 + llr0_r1 + llr0_r2 + llr0_r3 + llr0_r4 + lpreodds_r1 + lpreodds_r2 + lpreodds_r3 + lpreodds_r4, data=regdata[regdata$restrict2==1,])
    fm.ar.ct$var <- clx(fm.ar.ct, 1, regdata$subjectid[regdata$restrict2==1])
    restrictmat.llr1 <- rbind(c(1,-1,0,0,0,0,0,0,0,0,0,0), c(0,1,-1,0,0,0,0,0,0,0,0,0), c(0,0,1,-1,0,0,0,0,0,0,0,0))
    linearHypothesis(fm.ar.ct, restrictmat.llr1, vcov=fm.ar.ct$var)
    restrictmat.llr0 <- rbind(c(0,0,0,0,1,-1,0,0,0,0,0,0), c(0,0,0,0,0,1,-1,0,0,0,0,0), c(0,0,0,0,0,0,1,-1,0,0,0,0))
    linearHypothesis(fm.ar.ct, restrictmat.llr0, vcov=fm.ar.ct$var)
    restrictmat.lpreodds <- rbind(c(0,0,0,0,0,0,0,0,1,-1,0,0), c(0,0,0,0,0,0,0,0,0,1,-1,0), c(0,0,0,0,0,0,0,0,0,0,1,-1))
    linearHypothesis(fm.ar.ct, restrictmat.lpreodds, vcov=fm.ar.ct$var)

    ## IV
    attach(regdata[regdata$restrict2==1,])
    fm.ar.ct.iv <- gmm(lpostodds ~ 0 + llr1_r1 + llr1_r2 + llr1_r3 + llr1_r4 + llr0_r1 + llr0_r2 + llr0_r3 + llr0_r4 + lpreodds_r1 + lpreodds_r2 + lpreodds_r3 + lpreodds_r4, 
                        cbind(llr1_r1, llr1_r2, llr1_r3, llr1_r4, llr0_r1, llr0_r2, llr0_r3, llr0_r4, othermeanscore_r1, othermeanscore_r2, othermeanscore_r3, othermeanscore_r4))
    detach(regdata[regdata$restrict2==1,])
    # reactions to positive signals    
    restrictmat.llr1 <- rbind(c(1,-1,0,0,0,0,0,0,0,0,0,0), c(0,1,-1,0,0,0,0,0,0,0,0,0), c(0,0,1,-1,0,0,0,0,0,0,0,0))
    linearHypothesis(fm.ar.ct.iv, restrictmat.llr1, vcov=fm.ar.ct.iv$vcov)
    # reactions to negative signals       
    restrictmat.llr0 <- rbind(c(0,0,0,0,1,-1,0,0,0,0,0,0), c(0,0,0,0,0,1,-1,0,0,0,0,0), c(0,0,0,0,0,0,1,-1,0,0,0,0))
    linearHypothesis(fm.ar.ct.iv, restrictmat.llr0, vcov=fm.ar.ct.iv$vcov)
    # persistence of prior   
    restrictmat.lpreodds <- rbind(c(0,0,0,0,0,0,0,0,1,-1,0,0), c(0,0,0,0,0,0,0,0,0,1,-1,0), c(0,0,0,0,0,0,0,0,0,0,1,-1))
    linearHypothesis(fm.ar.ct.iv, restrictmat.lpreodds, vcov=fm.ar.ct.iv$vcov)
    # all pooled
    restrictmat <- rbind(restrictmat.llr1, restrictmat.llr0, restrictmat.lpreodds)
    linearHypothesis(fm.ar.ct.iv, restrictmat, vcov=fm.ar.ct.iv$vcov)

### Table 2: sufficient statistics 
    
    # round-by-round tests for memorylessness (OLS)
    wideregdata <- reshape(regdata, direction="wide", timevar="round")
    fm.r2.lag <- lm(lpostodds.2 ~ 0 + llr1.2 + llr0.2 + lpreodds.2 + llr.1, data=wideregdata[wideregdata$restrict2.2==1,])
    coefficients(fm.r2.lag) / sqrt(diag(hccm(fm.r2.lag)))
    fm.r3.lag <- lm(lpostodds.3 ~ 0 + llr1.3 + llr0.3 + lpreodds.3 + llr.1 + llr.2, data=wideregdata[wideregdata$restrict2.3==1,])
    coefficients(fm.r3.lag) / sqrt(diag(hccm(fm.r3.lag)))
    fm.r4.lag <- lm(lpostodds.4 ~ 0 + llr1.4 + llr0.4 + lpreodds.4 + llr.1 + llr.2 + llr.3, data=wideregdata[wideregdata$restrict2.4==1,])
    coefficients(fm.r4.lag) / sqrt(diag(hccm(fm.r4.lag)))

    # round-by-round tests for memorylessness (IV)
    wideregdata <- reshape(regdata, direction="wide", timevar="round")
    # round 2    
    attach(wideregdata[wideregdata$restrict2.2==1,])
    fm.r2.lag <- gmm(lpostodds.2 ~ 0 + llr1.2 + llr0.2 + lpreodds.2 + llr.1, cbind(llr1.2, llr0.2, othermeanscore_r2.2, llr.1))
    detach(wideregdata[wideregdata$restrict2.2==1,])
    # round 3    
    attach(wideregdata[wideregdata$restrict2.3==1,])
    fm.r3.lag <- gmm(lpostodds.3 ~ 0 + llr1.3 + llr0.3 + lpreodds.3 + llr.1 + llr.2, cbind(llr1.3, llr0.3, othermeanscore_r3.3, llr.1, llr.2))
    detach(wideregdata[wideregdata$restrict2.3==1,])
    # round 4
    attach(wideregdata[wideregdata$restrict2.4==1,])
    fm.r4.lag <- gmm(lpostodds.4 ~ 0 + llr1.4 + llr0.4 + lpreodds.4 + llr.1 + llr.2 + llr.3, cbind(llr1.4, llr0.4, othermeanscore_r4.4, llr.1, llr.2, llr.3))
    detach(wideregdata[wideregdata$restrict2.4==1,])

    # kludge to harmonize the names
    names(fm.r2.lag$coefficients) <- rownames(fm.r2.lag$vcov) <- colnames(fm.r2.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l1")
    names(fm.r3.lag$coefficients) <- rownames(fm.r3.lag$vcov) <- colnames(fm.r3.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l2", "llr.l1")
    names(fm.r4.lag$coefficients) <- rownames(fm.r4.lag$vcov) <- colnames(fm.r4.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l3", "llr.l2", "llr.l1")
    
    # create an output table
    vars.gmm.lags <- c("lpreodds","llr1","llr0","llr.l1","llr.l2","llr.l3")
    varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                          "lpreodds"            = "$\\delta$",
                          "llr1"                = "$\\beta_H$",
                          "llr0"                = "$\\beta_L$",
                          "llr.l1"             = "$\\beta_{-1}$",
                          "llr.l2"             = "$\\beta_{-2}$",
                          "llr.l3"             = "$\\beta_{-3}$")
    table.gmm.lags <- multiregtable(vars.gmm.lags, varlabels, list(fm.r2.lag, fm.r3.lag, fm.r4.lag), 3)
    
    # output combined table
    result <- hacktex(table.gmm.lags,
                    file=paste(tabdir, "gmm_lags.tex", sep="/"),
                    label="tab:gmm_lags",
                    table.env=FALSE,
                    caption.loc="top",
                    center="none",
                    rowlabel="",
                    rowname=rep("",14),
                    colheads=c("Regressor","Round 2","Round 3","Round 4"),
                    collabel.just=c("l","c","c","c"))

### Footnote 29: tests for stability: full sample
    
    # test equality of coefficients across rounds
    regdata$llr0_r1 <- regdata$llr0 * as.numeric(regdata$round==1)
    regdata$llr0_r2 <- regdata$llr0 * as.numeric(regdata$round==2)
    regdata$llr0_r3 <- regdata$llr0 * as.numeric(regdata$round==3)
    regdata$llr0_r4 <- regdata$llr0 * as.numeric(regdata$round==4)
    regdata$llr1_r1 <- regdata$llr1 * as.numeric(regdata$round==1)
    regdata$llr1_r2 <- regdata$llr1 * as.numeric(regdata$round==2)
    regdata$llr1_r3 <- regdata$llr1 * as.numeric(regdata$round==3)
    regdata$llr1_r4 <- regdata$llr1 * as.numeric(regdata$round==4)
    regdata$lpreodds_r1 <- regdata$lpreodds * as.numeric(regdata$round==1)
    regdata$lpreodds_r2 <- regdata$lpreodds * as.numeric(regdata$round==2)
    regdata$lpreodds_r3 <- regdata$lpreodds * as.numeric(regdata$round==3)
    regdata$lpreodds_r4 <- regdata$lpreodds * as.numeric(regdata$round==4)
    regdata$othermeanscore_r1 <- regdata$othermeanscore * as.numeric(regdata$round==1)
    regdata$othermeanscore_r2 <- regdata$othermeanscore * as.numeric(regdata$round==2)
    regdata$othermeanscore_r3 <- regdata$othermeanscore * as.numeric(regdata$round==3)
    regdata$othermeanscore_r4 <- regdata$othermeanscore * as.numeric(regdata$round==4)
    
    ## OLS
    fm.ar.ct <- lm(lpostodds ~ 0 + llr1_r1 + llr1_r2 + llr1_r3 + llr1_r4 + llr0_r1 + llr0_r2 + llr0_r3 + llr0_r4 + lpreodds_r1 + lpreodds_r2 + lpreodds_r3 + lpreodds_r4, data=regdata)
    fm.ar.ct$var <- clx(fm.ar.ct, 1, regdata$subjectid)
    restrictmat.llr1 <- rbind(c(1,-1,0,0,0,0,0,0,0,0,0,0), c(0,1,-1,0,0,0,0,0,0,0,0,0), c(0,0,1,-1,0,0,0,0,0,0,0,0))
    linearHypothesis(fm.ar.ct, restrictmat.llr1, vcov=fm.ar.ct$var)
    restrictmat.llr0 <- rbind(c(0,0,0,0,1,-1,0,0,0,0,0,0), c(0,0,0,0,0,1,-1,0,0,0,0,0), c(0,0,0,0,0,0,1,-1,0,0,0,0))
    linearHypothesis(fm.ar.ct, restrictmat.llr0, vcov=fm.ar.ct$var)
    restrictmat.lpreodds <- rbind(c(0,0,0,0,0,0,0,0,1,-1,0,0), c(0,0,0,0,0,0,0,0,0,1,-1,0), c(0,0,0,0,0,0,0,0,0,0,1,-1))
    linearHypothesis(fm.ar.ct, restrictmat.lpreodds, vcov=fm.ar.ct$var)
    
    ## IV
    remove(lpostodds)
    attach(regdata)
    fm.ar.ct.iv <- gmm(lpostodds ~ 0 + llr1_r1 + llr1_r2 + llr1_r3 + llr1_r4 + llr0_r1 + llr0_r2 + llr0_r3 + llr0_r4 + lpreodds_r1 + lpreodds_r2 + lpreodds_r3 + lpreodds_r4, 
                       cbind(llr1_r1, llr1_r2, llr1_r3, llr1_r4, llr0_r1, llr0_r2, llr0_r3, llr0_r4, othermeanscore_r1, othermeanscore_r2, othermeanscore_r3, othermeanscore_r4))
    detach(regdata)
    # reactions to positive signals    
    restrictmat.llr1 <- rbind(c(1,-1,0,0,0,0,0,0,0,0,0,0), c(0,1,-1,0,0,0,0,0,0,0,0,0), c(0,0,1,-1,0,0,0,0,0,0,0,0))
    linearHypothesis(fm.ar.ct.iv, restrictmat.llr1, vcov=fm.ar.ct.iv$vcov)
    # reactions to negative signals       
    restrictmat.llr0 <- rbind(c(0,0,0,0,1,-1,0,0,0,0,0,0), c(0,0,0,0,0,1,-1,0,0,0,0,0), c(0,0,0,0,0,0,1,-1,0,0,0,0))
    linearHypothesis(fm.ar.ct.iv, restrictmat.llr0, vcov=fm.ar.ct.iv$vcov)
    # persistence of prior   
    restrictmat.lpreodds <- rbind(c(0,0,0,0,0,0,0,0,1,-1,0,0), c(0,0,0,0,0,0,0,0,0,1,-1,0), c(0,0,0,0,0,0,0,0,0,0,1,-1))
    linearHypothesis(fm.ar.ct.iv, restrictmat.lpreodds, vcov=fm.ar.ct.iv$vcov)
    # all pooled
    restrictmat <- rbind(restrictmat.llr1, restrictmat.llr0, restrictmat.lpreodds)
    linearHypothesis(fm.ar.ct.iv, restrictmat, vcov=fm.ar.ct.iv$vcov)
    
### Table S-5: sufficient statistics: full sample
    
    # round-by-round tests for memorylessness (OLS)
    wideregdata <- reshape(regdata, direction="wide", timevar="round")
    fm.r2.lag <- lm(lpostodds.2 ~ 0 + llr1.2 + llr0.2 + lpreodds.2 + llr.1, data=wideregdata)
    coefficients(fm.r2.lag) / sqrt(diag(hccm(fm.r2.lag)))
    fm.r3.lag <- lm(lpostodds.3 ~ 0 + llr1.3 + llr0.3 + lpreodds.3 + llr.1 + llr.2, data=wideregdata)
    coefficients(fm.r3.lag) / sqrt(diag(hccm(fm.r3.lag)))
    fm.r4.lag <- lm(lpostodds.4 ~ 0 + llr1.4 + llr0.4 + lpreodds.4 + llr.1 + llr.2 + llr.3, data=wideregdata)
    coefficients(fm.r4.lag) / sqrt(diag(hccm(fm.r4.lag)))
    
    # round-by-round tests for memorylessness (IV)
    wideregdata <- reshape(regdata, direction="wide", timevar="round")
    # round 2    
    attach(wideregdata)
    fm.r2.lag <- gmm(lpostodds.2 ~ 0 + llr1.2 + llr0.2 + lpreodds.2 + llr.1, cbind(llr1.2, llr0.2, othermeanscore_r2.2, llr.1))
    detach(wideregdata)
    # round 3    
    attach(wideregdata)
    fm.r3.lag <- gmm(lpostodds.3 ~ 0 + llr1.3 + llr0.3 + lpreodds.3 + llr.1 + llr.2, cbind(llr1.3, llr0.3, othermeanscore_r3.3, llr.1, llr.2))
    detach(wideregdata)
    # round 4
    attach(wideregdata)
    fm.r4.lag <- gmm(lpostodds.4 ~ 0 + llr1.4 + llr0.4 + lpreodds.4 + llr.1 + llr.2 + llr.3, cbind(llr1.4, llr0.4, othermeanscore_r4.4, llr.1, llr.2, llr.3))
    detach(wideregdata)
    
    # kludge to harmonize the names
    names(fm.r2.lag$coefficients) <- rownames(fm.r2.lag$vcov) <- colnames(fm.r2.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l1")
    names(fm.r3.lag$coefficients) <- rownames(fm.r3.lag$vcov) <- colnames(fm.r3.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l2", "llr.l1")
    names(fm.r4.lag$coefficients) <- rownames(fm.r4.lag$vcov) <- colnames(fm.r4.lag$vcov) <- c("llr1", "llr0", "lpreodds", "llr.l3", "llr.l2", "llr.l1")
    
    # create an output table
    vars.gmm.lags <- c("lpreodds","llr1","llr0","llr.l1","llr.l2","llr.l3")
    varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                          "lpreodds"            = "$\\delta$",
                          "llr1"                = "$\\beta_H$",
                          "llr0"                = "$\\beta_L$",
                          "llr.l1"             = "$\\beta_{-1}$",
                          "llr.l2"             = "$\\beta_{-2}$",
                          "llr.l3"             = "$\\beta_{-3}$")
    table.gmm.lags <- multiregtable(vars.gmm.lags, varlabels, list(fm.r2.lag, fm.r3.lag, fm.r4.lag), 3)
    
    # output combined table
    result <- hacktex(table.gmm.lags,
                      file=paste(tabdir, "gmm_lags_fullsample.tex", sep="/"),
                      label="tab:gmm_lags_fullsample",
                      table.env=FALSE,
                      caption.loc="top",
                      center="none",
                      rowlabel="",
                      rowname=rep("",14),
                      colheads=c("Regressor","Round 2","Round 3","Round 4"),
                      collabel.just=c("l","c","c","c"))
    
### Table S-7: allowing for interactions between prior and signal

    # Round-by-round plus pooled OLS regressions
    fm.r1.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata[regdata$round==1 & regdata$restrict2==1,])
    fm.r1.per$var <- hccm(fm.r1.per)
    fm.r2.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata[regdata$round==2 & regdata$restrict2==1,])
    fm.r2.per$var <- hccm(fm.r2.per)
    fm.r3.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata[regdata$round==3 & regdata$restrict2==1,])
    fm.r3.per$var <- hccm(fm.r3.per)
    fm.r4.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata[regdata$round==4 & regdata$restrict2==1,])
    fm.r4.per$var <- hccm(fm.r4.per)
    fm.ar.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata[regdata$restrict2==1,])
    fm.ar.per$var <- clx(fm.ar.per, 1, regdata$subjectid[regdata$restrict2==1])
    fm.ar.ur.per <- lm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, data=regdata)
    fm.ar.ur.per$var <- clx(fm.ar.ur.per, 1, regdata$subjectid)

    # create OLS output table pane
    varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                          "lpreodds"            = "$\\delta$",
                          "llr1"                = "$\\beta_H$",
                          "llr0"                = "$\\beta_L$",
                          "lpreodds:signal"       = "$\\delta_H$")
    vars.base.per <- c("lpreodds","lpreodds:signal","llr1","llr0")
    table.base.per <- multiregtable(vars.base.per, varlabels, list(fm.r1.per,fm.r2.per,fm.r3.per,fm.r4.per,fm.ar.per,fm.ar.ur.per), 3)

    # IV regressions
    estgmm.per <- function(restriction){
    
        lpostodds <- regdata$lpostodds[restriction]
        llr1 <- regdata$llr1[restriction]
        llr0 <- regdata$llr0[restriction]
        lpreodds <- regdata$lpreodds[restriction]
        quiztype <- regdata$quiztype[restriction]
        othermeanscore <- regdata$othermeanscore[restriction]
        signal <- regdata$signal[restriction]
        return(gmm(lpostodds ~ 0 + llr1 + llr0 + lpreodds + lpreodds:signal, cbind(llr1, llr0, othermeanscore, othermeanscore*signal)))
        
    }
    
    # estimate GMM models
    gmmmodels.per <- list()
    gmmfstats.per <- list()
    for (r in 1:4){
        gmmmodels.per[[r]] <- estgmm.per(regdata$restrict2==1 & regdata$round==r)
    }
    gmmmodels.per[[5]] <- estgmm.per(regdata$restrict2==1)
    gmmmodels.per[[6]] <- estgmm.per(TRUE)

    # create output table
    vars.gmm.per <- c("lpreodds","lpreodds:signal","llr1","llr0")
    table.gmm.per <- multiregtable(vars.gmm.per, varlabels, gmmmodels.per, 3)
    
    # output combined table
    result <- hacktex(rbind(table.base.per, table.gmm.per), 
                    file=paste(tabdir, "base_combo_per.tex", sep="/"),
                    label="tab:base_combo_per",
                    table.env=FALSE,
                    caption.loc="top",
                    center="none",
                    rowlabel="",
                    rowname=rep("",16),
                    rgroup=c("Panel A: OLS","Panel B: IV"),
                    n.rgroup=c(10,10),
                    colheads=c("Regressor","Round 1","Round 2","Round 3","Round 4","All Rounds","Unrestricted"),
                    collabel.just=c("l","c","c","c","c","c","c"))
               
### Table 3 (heterogeneity in updating by ability) and Table S-1 (heterogeneity in updating by gender)

    # OLS
    fm.g.r2 <- lm(lpostodds ~ 0 + llr1 + llr1_male + llr0 + llr0_male + lpreodds + lpreodds_male, data=regdata[regdata$restrict2==1,])
    fm.g.r2$var <- clx(fm.g.r2, 1, regdata$subjectid[regdata$restrict2==1])
    fm.a.r2 <- lm(lpostodds ~ 0 + llr1 + llr1_tth + llr0 + llr0_tth + lpreodds + lpreodds_tth, data=regdata[regdata$restrict2==1,])
    fm.a.r2$var <- clx(fm.a.r2, 1, regdata$subjectid[regdata$restrict2==1])

    # OLS with delta = 1 imposed
    fm.g.r2 <- lm(lpostodds - lpreodds ~ 0 + llr1 + llr1_male + llr0 + llr0_male, data=regdata[regdata$restrict2==1,])
    fm.g.r2$var <- clx(fm.g.r2, 1, regdata$subjectid[regdata$restrict2==1])
    fm.a.r2 <- lm(lpostodds - lpreodds ~ 0 + llr1 + llr1_tth + llr0 + llr0_tth, data=regdata[regdata$restrict2==1,])
    fm.a.r2$var <- clx(fm.a.r2, 1, regdata$subjectid[regdata$restrict2==1])
    
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
    lpreodds_tth <- regdata$lpreodds_tth[regdata$restrict2==1]
    llr1_tth <- regdata$llr1_tth[regdata$restrict2==1]
    llr0_tth <- regdata$llr0_tth[regdata$restrict2==1]
    regdata$oms_tth <- regdata$othermeanscore * regdata$temptophalf
    oms_tth <- regdata$oms_tth[regdata$restrict2==1]
    
    # estimation
    gmm.g.r2 <- gmm(lpostodds ~ 0 + llr1 + llr0 + llr1_male + llr0_male + lpreodds + lpreodds_male, cbind(llr1,llr0,llr1_male,llr0_male,oms,oms_male))
    gmm.a.r2 <- gmm(lpostodds ~ 0 + llr1 + llr0 + llr1_tth + llr0_tth + lpreodds + lpreodds_tth, cbind(llr1,llr0,llr1_tth,llr0_tth,oms,oms_tth))
    
    # create the differential regression table
    varlabels <- pairlist("lpostodds"           = "Log Posterior Odds",
                          "llr1"                = "$\\beta_H$",
                          "llr0"                = "$\\beta_L$",
                          "llr1_male"           = "$\\beta_H^{Male}$",
                          "llr0_male"           = "$\\beta_L^{Male}$",
                          "llr1_tth"            = "$\\beta_H^{Able}$",
                          "llr0_tth"            = "$\\beta_L^{Able}$",
                          "lpreodds"            = "$\\delta$",
                          "lpreodds_male"       = "$\\delta^{Male}$",
                          "lpreodds_tth"        = "$\\delta^{Able}$")
    vars.gender <- c("llr1","llr0","llr1_male","llr0_male")
    table.gender <- multiregtable(vars.gender, varlabels, list(fm.g.r2), 3)
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
    vars.ability <- c("llr1","llr0","llr1_tth","llr0_tth")
    table.ability <- multiregtable(vars.ability, varlabels, list(fm.a.r2), 3)
    result <- hacktex(table.ability, 
                    file=paste(tabdir, "differential_ability.tex", sep="/"),
                    label="tab:differential_ability",
                    tabwidth="3in",
                    table.env=FALSE,
                    caption.loc="top",
                    rowname =NULL,
                    center="none",
                    colheads=c("Regressor","OLS"),
                    collabel.just=c("l","c"))
    
### Section S-1: miscellaneous tests and calculations derived from Table S-1
    
    # percent changes in conservatism and relative asymmetry
    perchangecons <- function(model){
        return(((model$coefficients["llr1"] + model$coefficients["llr1_male"])/model$coefficients["llr1"] + (model$coefficients["llr0"] + model$coefficients["llr0_male"])/model$coefficients["llr0"])/2 - 1)
    }
    perchangecons(fm.g.r2)
    perchangecons(gmm.g.r2)
    perchangeasym <- function(model){
        return((model$coefficients["llr1"] + model$coefficients["llr1_male"])/(model$coefficients["llr0"] + model$coefficients["llr0_male"]) - model$coefficients["llr1"]/model$coefficients["llr0"])
    }
    perchangeasym(fm.g.r2)
    perchangeasym(gmm.g.r2)
    
    # Wald tests for change in relative asymmetry
    relasymtest.gender <- function(model){
        h <- model$coefficients["llr1_male"]*model$coefficients["llr0"] - model$coefficients["llr0_male"]*model$coefficients["llr1"]
        gradh <- c(- model$coefficients["llr0_male"],model$coefficients["llr1_male"],model$coefficients["llr0"],-model$coefficients["llr1"])
        if (!is.null(model$var)){
            v <- t(gradh) %*% model$var[c("llr1","llr0","llr1_male","llr0_male"),c("llr1","llr0","llr1_male","llr0_male")] %*% gradh
        } else if (!is.null(model$vcov)){
            v <- t(gradh) %*% model$vcov[c("llr1","llr0","llr1_male","llr0_male"),c("llr1","llr0","llr1_male","llr0_male")] %*% gradh
        } else {
            return("ERROR: couldn't find covariance matrix.")
        }
        teststat <- h * (1 / v) * h
        1 - pchisq(teststat, df=1, ncp=0)
    }
    relasymtest.gender(fm.g.r2)
    relasymtest.gender(gmm.g.r2)
      
### Table 4: belief updating about own vs robot performance

    # drop guys whose signals != musignals
    bans <- unique(muregdata$subjectid[muregdata$signal != muregdata$musignal])
    pairdata <- muregdata[!is.element(muregdata$subjectid, bans),]
    pairdata$delta <- pairdata$muloddschange - pairdata$loddschange
    
    # stack data from the main and make-up experiments for estimation
    bdata <- pairdata[,c("subjectid","llr1","llr0","loddschange","male","temptophalf","restrict2","murestrict2","restrict1","murestrict1","round")]
    bdata$experiment <- "baseline"
    mdata <- pairdata[,c("subjectid","llr1","llr0","muloddschange","male","temptophalf","restrict2","murestrict2","restrict1","murestrict1","round")]
    mdata$experiment <- "mu"
    names(mdata)[4] <- "loddschange"
    stackeddata <- rbind(bdata, mdata)
    stackeddata$llr1_mu <- stackeddata$llr1 * as.numeric(stackeddata$experiment=="mu")
    stackeddata$llr0_mu <- stackeddata$llr0 * as.numeric(stackeddata$experiment=="mu")
    stackeddata$bothrestrict2 <- stackeddata$restrict2 & stackeddata$murestrict2
    stackeddata$bothrestrict1 <- stackeddata$restrict1 & stackeddata$murestrict1

    # regressions
    fm.rob.1 <- lm(loddschange ~ 0+ llr1 + llr0 + llr1_mu + llr0_mu, data=stackeddata[stackeddata$bothrestrict2==1,])
    fm.rob.1$var <- clx(fm.rob.1, 1, stackeddata$subjectid[stackeddata$bothrestrict2==1])
    fm.rob.2 <- lm(loddschange ~ 0+ llr1 + llr0 + llr1_mu + llr0_mu, data=stackeddata[stackeddata$bothrestrict1==1,])
    fm.rob.2$var <- clx(fm.rob.2, 1, stackeddata$subjectid[stackeddata$bothrestrict1==1])
    fm.rob.3 <- lm(loddschange ~ 0+ llr1 + llr0 + llr1_mu + llr0_mu, data=stackeddata)
    fm.rob.3$var <- clx(fm.rob.3, 1, stackeddata$subjectid)
    
    # hypothesis tests
    htests <- matrix(NA, nrow=4, ncol=4)
    htests[,1] <- c("$\\mathbb{P}(\\beta_H + \\beta_H^{Robot} = 1)$","$\\mathbb{P}(\\beta_L + \\beta_L^{Robot} = 1)$","$\\mathbb{P}(\\beta_H  = \\beta_L)$","$\\mathbb{P}(\\beta_H + \\beta_H^{Robot} = \\beta_L + \\beta_L^{Robot})$")
    htests[1,2] <- roundsig(linearHypothesis(fm.rob.1, c(1,0,1,0), rhs=c(1), vcov=fm.rob.1$var)["Pr(>F)"][[1]][2],3)
    htests[1,3] <- roundsig(linearHypothesis(fm.rob.2, c(1,0,1,0), rhs=c(1), vcov=fm.rob.2$var)["Pr(>F)"][[1]][2],3)
    htests[1,4] <- roundsig(linearHypothesis(fm.rob.3, c(1,0,1,0), rhs=c(1), vcov=fm.rob.3$var)["Pr(>F)"][[1]][2],3)
    htests[2,2] <- roundsig(linearHypothesis(fm.rob.1, c(0,1,0,1), rhs=c(1), vcov=fm.rob.1$var)["Pr(>F)"][[1]][2],3)
    htests[2,3] <- roundsig(linearHypothesis(fm.rob.2, c(0,1,0,1), rhs=c(1), vcov=fm.rob.2$var)["Pr(>F)"][[1]][2],3)
    htests[2,4] <- roundsig(linearHypothesis(fm.rob.3, c(0,1,0,1), rhs=c(1), vcov=fm.rob.3$var)["Pr(>F)"][[1]][2],3)
    htests[3,2] <- roundsig(linearHypothesis(fm.rob.1, c(1,-1,0,0), vcov=fm.rob.1$var)["Pr(>F)"][[1]][2],3)
    htests[3,3] <- roundsig(linearHypothesis(fm.rob.2, c(1,-1,0,0), vcov=fm.rob.2$var)["Pr(>F)"][[1]][2],3)
    htests[3,4] <- roundsig(linearHypothesis(fm.rob.3, c(1,-1,0,0), vcov=fm.rob.3$var)["Pr(>F)"][[1]][2],3)
    htests[4,2] <- roundsig(linearHypothesis(fm.rob.1, c(1,-1,1,-1), vcov=fm.rob.1$var)["Pr(>F)"][[1]][2],3)
    htests[4,3] <- roundsig(linearHypothesis(fm.rob.2, c(1,-1,1,-1), vcov=fm.rob.2$var)["Pr(>F)"][[1]][2],3)
    htests[4,4] <- roundsig(linearHypothesis(fm.rob.3, c(1,-1,1,-1), vcov=fm.rob.3$var)["Pr(>F)"][[1]][2],3)
    
    # create the base regression table
    varlabels.mu <- pairlist("lpostodds"           = "Log Posterior Odds",
                              "lpreodds"            = "$\\delta$",
                              "llr1"                = "$\\beta_H$",
                              "llr0"                = "$\\beta_L$",
                              "llr1_mu"             = "$\\beta_H^{Robot}$",
                              "llr0_mu"             = "$\\beta_L^{Robot}$")
    vars.mu <- c("llr1","llr0","llr1_mu","llr0_mu")
    table.mu <- multiregtable(vars.mu, varlabels.mu, list(fm.rob.1, fm.rob.2, fm.rob.3), 3, addrows=htests)
    
    # output
    result <- hacktex(table.mu,
                    file=paste(tabdir, "updating_mu.tex", sep="/"),
                    label="tab:updating_mu",
                    table.env=FALSE,
                    caption.loc="top",
                    center="none",
                    rowname =NULL,
                    colheads=c("Regressor","I","II","III"),
                    collabel.just=c("l","c","c","c"))
    
### Section 4.2: costs of deviations from Bayes' rule (relative to benchmark of no updating)

    # in main experiment, benchmarked against no updating
    data$profit_noupdating <- (0.5 - 0.5 * (data$preconfidence1/100)^2 + (data$preconfidence1/100) * data$temptophalf
                            + 0.5 - 0.5 * (data$preconfidence1/100)^2 + (data$preconfidence1/100) * data$temptophalf
                            + 0.5 - 0.5 * (data$preconfidence1/100)^2 + (data$preconfidence1/100) * data$temptophalf
                            + 0.5 - 0.5 * (data$preconfidence1/100)^2 + (data$preconfidence1/100) * data$temptophalf)
    mean(data$profit_bayes, na.rm=T)
    mean(data$profit_actual, na.rm=T)
    mean(data$profit_noupdating, na.rm=T)
    (mean(data$profit_actual, na.rm=T) - mean(data$profit_noupdating, na.rm=T))/(mean(data$profit_bayes, na.rm=T) - mean(data$profit_noupdating, na.rm=T))    
    
    (mean(data$profit_bayes, na.rm=T) - mean(data$profit_actual, na.rm=T))/(mean(data$profit_bayes, na.rm=T) - mean(data$profit_noupdating, na.rm=T))    
    
  
### Table S-2: implied valuations for information

    # fraction of observations censored
    data$censorings <- (as.numeric(data$talkingbob == 0 | data$talkingbob==400)
                        + as.numeric(data$precisebob == 0 | data$precisebob==400))
    sum(data$censorings[data$restrict2==1], na.rm=T) / (2 * nrow(data[data$restrict2==1,]))
    
    # define helper function to calculate the stats
    sumstats <- function(dvec){
        
        retmat <- array(0, dim=c(1,4))
        retmat[1,1] <- NROW(dvec[!is.na(dvec)])
        retmat[1,2] <- roundsig(mean(dvec, na.rm = TRUE), 1)
        retmat[1,3] <- roundsig(sd(dvec, na.rm = TRUE), 1)
        retmat[1,4] <- roundsig(NROW(dvec[dvec<0]) / NROW(dvec), 2)
        
        return(retmat)
    }
    
    # for reference, how prominent is aversion in the full sample
    sumstats(data$coarse)
    
    # how many underbids for $2 were there?
    plot(density(data$silentbob[data$restrict2==1 & !is.na(data$silentbob)]))
    table(data$silentbob[data$restrict2==1])
    table(as.numeric(data$silentbob[data$restrict2==1] < 199))/NROW(data$silentbob[data$restrict2==1])
    
    # focus on restrict2 guys
    usedata <- data[data$restrict2==1,]
    
    # assemble the table
    output <- array(data="", dim=c(9,4))
    output[2,1:4] <- sumstats(usedata$coarse)
    output[3,1:4] <- sumstats(usedata$precise)
    output[5,1:4] <- sumstats(usedata$coarse[usedata$male==0])
    output[6,1:4] <- sumstats(usedata$precise[usedata$male==0])
    output[8,1:4] <- sumstats(usedata$coarse[usedata$male==1])
    output[9,1:4] <- sumstats(usedata$precise[usedata$male==1])
    
    nrow(usedata)
    nrow(usedata[usedata$male==1,])
    nrow(usedata[usedata$male==0,])
    
    rnames <- rep(c("","\\hspace{1em}Learning top/bottom half","\\hspace{1em}Learning percentile"), 3)
    rnames[1] <- "\\textbf{Estimation Sample}"
    rnames[4] <- "\\textbf{Women}"
    rnames[7] <- "\\textbf{Men}"
    
    # output
    result <- hacktex(output, 
        file=paste(tabdir, "value_summary_stats.tex", sep="/"),
        title="",
        label="tab:value_summary_stats",
        table.env=FALSE,
        colheads=c("$N$", "Mean","Std. Dev.","$P(v< 0)$"),
        rowname=rnames,
        collabel.just=c("c","c","c","c"))
    
### Table S-3: confidence and positive information value      

    # helper function for computing heteroskedasticity-robust SEs in IV models
    ivrobcov <- function(fm){
        X <- fm$x[["regressors"]]
        Z <- fm$x[["instruments"]]
        Sigma <- diag(fm$residuals * fm$residuals)
        fm$var <- solve((t(X) %*% Z) %*% solve(t(Z) %*% Sigma %*% Z) %*% (t(Z) %*% X))
        return(fm)
    }

    # prep data
    usedata <- data[data$restrict2==1 & !is.na(data$coarse),]
    usedata$lpostodds4 <- log(usedata$postconfidence4 / (100 - usedata$postconfidence4))
    usedata <- usedata[is.finite(usedata$lpostodds4),]
    usedata$lpostodds4.2 <- usedata$lpostodds4 * usedata$lpostodds4
    
    # estimate models
    fm.1 <- robcov(ols(as.numeric(coarse>=0) ~ lpostodds4, data=usedata,x=T,y=T))
    fm.2 <- robcov(ols(as.numeric(coarse>=0) ~ lpostodds4 + temptophalf, data=usedata,x=T,y=T))
    fm.3 <- robcov(ols(as.numeric(coarse>=0) ~ lpostodds4 + temptophalf + male + year, data=usedata,x=T,y=T))
    fm.4 <- ivrobcov(ivreg(as.numeric(coarse>=0) ~ lpostodds4 + temptophalf, ~ sigsum + othermeanscore + temptophalf, data=usedata,x=T,y=T))
    fm.5 <- ivrobcov(ivreg(as.numeric(coarse>=0) ~ lpostodds4 + temptophalf + male + year, ~ sigsum + othermeanscore + temptophalf + male + year, data=usedata,x=T,y=T))

    # calculate first-stage F-statistics (note kludge as linearHypothesis() doesn't play nicely with ols())
    fs.4 <- robcov(ols(lpostodds4 ~ sigsum + othermeanscore + temptophalf, data=usedata, x=T,y=T))
    fs.4.lm <- lm(lpostodds4 ~ sigsum + othermeanscore + temptophalf, data=usedata, x=T,y=T)
    fstat.4 <- roundsig(linearHypothesis(fs.4.lm, rbind(c(0,1,0,0), c(0,0,1,0)), vcov=fs.4$var, test="F")$F[2], 2)
    fs.5 <- robcov(ols(lpostodds4 ~ sigsum + othermeanscore + temptophalf + male + year, data=usedata, x=T,y=T))
    fs.5.lm <- lm(lpostodds4 ~ sigsum + othermeanscore + temptophalf + male + year, data=usedata, x=T,y=T)
    fstat.5 <- roundsig(linearHypothesis(fs.5.lm, rbind(c(0,1,0,0,0,0), c(0,0,1,0,0,0)), vcov=fs.5$var, test="F")$F[2], 2)

    varlabels <- pairlist("lpostodds4"           = "$\\logit(\\mu)$",
                           "male"                = "Male",
                           "year"                = "YOG",
                           "temptophalf"         = "Top Half")
    vars.info <- c("lpostodds4","temptophalf","male","year")
    addrows.info <- rbind(c("First-Stage $F$-Statistic","-","-","-",fstat.4,fstat.5))
    table.info <- multiregtable(vars.info, varlabels, list(fm.1,fm.2,fm.3,fm.4,fm.5), 3, addrows=addrows.info)
    
    # output combined table
    result <- hacktex(table.info, 
                    file=paste(tabdir, "infovalue_regs_posterior.tex", sep="/"),
                    label="tab:infovalue_regs_posterior",
                    table.env=FALSE,
                    caption.loc="top",
                    rowname=NULL,
                    center="none",
                    cgroup=c("","OLS","IV"),
                    n.cgroup=c(1,3,2),
                    colheads=c("Regressor","I","II","III","IV","V"),
                    collabel.just=c("l","c","c","c","c","c"))

### Section S-3.1: info value noise test
    
    usedata <- data[data$restrict2==1,]
    
    # intuitive check
    table(as.numeric(usedata$coarse < 0)) / nrow(usedata[!is.na(usedata$coarse),])
    table(as.numeric(usedata$precise < 0)) / nrow(usedata[!is.na(usedata$precise),])
    table(as.numeric(usedata$coarse < 0), as.numeric(usedata$precise < 0))
    cor(as.numeric(usedata$coarse < 0), as.numeric(usedata$precise < 0), use="complete.obs")
    
    # calculate sample second moments
    hatvs <- var(usedata$silentbob, na.rm=TRUE)
    hatvc <- var(usedata$talkingbob, na.rm=TRUE)
    hatvp <- var(usedata$precisebob, na.rm=TRUE)
    hatccs <- cov(usedata$talkingbob, usedata$silentbob, use="complete.obs")
    hatccp <- cov(usedata$talkingbob, usedata$precisebob, use="complete.obs")
    hatcps <- cov(usedata$precisebob, usedata$silentbob, use="complete.obs")
    k1 <- hatccp + hatvs - hatcps - hatccs
    k2 <- hatvc + hatvs - 2*hatccs
    k3 <- hatvp + hatvs - 2*hatcps
    
    # parameters must satisfy veps <= min{k2/2,k3/2} and posfun(veps) >=0
    # a conservative value for sdeps is sqrt(699) using restrict2 sample
    posfun <- function(veps){
        return(sqrt((k2 - 2 *veps)*(k3 - 2*veps)) + veps - k1)
    }
    plot(1:10000, posfun(1:10000), type="l")
    posfun(699)
    posfun(700)
    k2/2
    k3/2
    sdeps <- sqrt(699)
    
    # check the second inequality restriction on the correlation of C and P
    posfun2 <- function(veps){
        return(sqrt((k2 - 2 *veps)*(k3 - 2*veps)) + k1 - veps)
    }
    plot(1:10000, posfun2(1:10000), type="l")
    
    # find p-values for a range of thresholds
    ps.c <- rep(NA, 15)
    ps.j <- rep(NA, 15)
    for (t in 1:15){
        
        x <- - t * 10
        print(x)
        ### test using only coarse bids
        n <- nrow(usedata[!is.na(usedata$talkingbob),])
        nx <- nrow(usedata[!is.na(usedata$talkingbob) & usedata$coarse < x,])
        Fx <- pnorm(x, sd=sdeps*sqrt(2))
        
        # easier to estimate the conjugate probability -- fewer calculcations
        pval <- 0
        for (m in 0:(nx)){
            pval <- pval + ((Fx)^m)*((1 - Fx)^(n-m)) * choose(n, m)
            #print(pval)
        }
        ps.c[t] <- 1 - pval
        
        ### joint coarse-precise tests
        n <- nrow(usedata[!is.na(usedata$talkingbob),])
        nx <- nrow(usedata[!is.na(usedata$talkingbob) & usedata$coarse < x & usedata$precise < x,])
        
        # need to calculate the probabability of one individual having both P < x and C < x numerically
        sims <- 100000
        res <- rep(NA, sims)
        for (s in 1:sims){
            epsvec <- rnorm(3, mean=0, sd=sdeps)
            res[s] <- as.numeric(epsvec[1] > max(epsvec[2:3]) - x)
        }
        Fx <- mean(res)
        
        # now calculate the probability of getting at most nx observations below x and get its conjugate
        pval <- 0
        for (m in 0:(nx)){
            pval <- pval + ((Fx)^m)*((1 - Fx)^(n-m)) * choose(n, m)
            #print(pval)
        }
        ps.j[t] <- 1 - pval
    }
    
    # plot max info value for comparison
    usedata$maxinfoval <- pmax(usedata$coarse, usedata$precise)
    maxvalcdf <- ecdf(usedata$maxinfoval)
    
    postscript(file=paste(figdir,"noisetest_rejections.eps",sep=""))
    plot(-(1:15)*10, ps.j, type="l", lty="solid", lwd=2, xlab="x", ylab="p")
    points(-(1:15)*10, ps.j)
    #lines(-(1:15)*10, maxvalcdf(-(1:15)*10), lty="dotted", lwd=2)
    dev.off()
