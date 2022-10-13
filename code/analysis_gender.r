#################################################################################################
#
#   analysis_gender.r
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
    library(rms)
    library(gplots)
    
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

    # log posterior odds
    data$lpreodds <- log(data$top50 / (100 - data$top50))
    data$lpostodds <- log(data$top50afters4 / (100 - data$top50afters4))
    data$lpostodds2 <- data$lpostodds * data$lpostodds
    
    # prep vars for plots
    data$qt_qscore <- interaction(data$quiztype, data$quizscore)
    data$llrafterquiz <- log(data$top50afterquiz / (100 - data$top50afterquiz))
    data$lpostodds <- log(data$top50afters4 / (100 - data$top50afters4))
    data$lpostodds2 <- data$lpostodds * data$lpostodds
    usedata <- data[data$restrict2==1,]

### Section S-1: statistics on initial confidence gap
    
    # differences in initial confidence, conditional and unconditional on performance
    summary(lm(usedata$top50 ~ usedata$male))
    summary(lm(usedata$quizscore ~ usedata$male))
    fm <- lm(top50 ~ as.factor(qt_qscore), data=usedata)
    summary(lm(fm$residuals ~ usedata$male[!is.na(usedata$top50)]))
    
    # differences in demand for information
    robcov(ols(as.numeric(coarse >=0) ~ male, data=usedata, x=T, y=T))
    robcov(ols(as.numeric(precise >=0) ~ male, data=usedata, x=T, y=T))

### Figure S-1: barplot of gender differences in demand/belief relationship

    # find quantiles of posterior distribution
    breakpoints <- 100 * (0:4)/4
    usedata$postquantile <- cut(usedata$top50afters4, breaks=breakpoints, include.lowest=TRUE)
    
    # cell means
    meanvalues.m <- tapply(usedata$coarse[usedata$male==1], usedata$postquantile[usedata$male==1], FUN=mean, na.rm=TRUE)
    meanvalues.w <- tapply(usedata$coarse[usedata$male==0], usedata$postquantile[usedata$male==0], FUN=mean, na.rm=TRUE)
    
    # cell standard deviations and confidence intervals
    meansd <- function(usedata, subset){
        sd <- tapply(usedata$coarse[subset], usedata$postquantile[subset], FUN=sd, na.rm=TRUE)
        n <- tapply(usedata$coarse[subset], usedata$postquantile[subset], FUN=NROW)
        return(sd / sqrt(n))
    }
    meansd.m <- meansd(usedata, usedata$male==1)
    meansd.w <- meansd(usedata, usedata$male==0)
    
    # plot
    postscript(file=paste(figdir,"infovalue_bars.eps",sep=""))
    barplot(rbind(meanvalues.m, meanvalues.w), beside=TRUE, ylim=c(0,75), 
            legend.text=c("Men","Women"),
            args.legend=list(x="topleft"),
            col=grey(c(0.8,0.3)),
            names.arg=c("0%-25%","26%-50%","51%-75%","76%-100%"))
    plotCI(x = (1:(NROW(breakpoints)-1))*3 - 1.5, y=meanvalues.m, uiw = 2 * meansd.m, add=TRUE, pch=NA, gap=0)
    plotCI(x = (1:(NROW(breakpoints)-1))*3 - 0.5, y=meanvalues.w, uiw = 2 * meansd.w, add=TRUE, pch=NA, gap=0)
    dev.off()

