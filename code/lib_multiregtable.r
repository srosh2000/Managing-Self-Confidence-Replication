####################################################################################
# multiregtable(varnames, varlabels, models, roundparam, addrows)
#
# Description
# -----------
# Take a list of fitted models and organize the results into a data frame
# with each model in a separate column.
#
# Arguments
# ---------
#       varnames        List of the variables to be included as rows
#       varlabels       Hash mapping variable names into labels for printing
#       models          List of fitted models
#       roundparam      Number of decimal places to round to
#       addrows         Any extra rows to insert in between the usual coef
#                       estimates and the N/R^2; for example, a row of Y/Ns
#                       indicating whether time FEs were included.
# Notes
# -----
# Calls helper functions multiregtable_addons(), multiregtable_sigsymbol().
#
####################################################################################

multiregtable <- function(varnames, varlabels, models, roundparam, addrows = NULL){

    rows <- NROW(varnames)
    cols <- NROW(models) + 1
    rettab <- as.data.frame(matrix(NA, nrow=rows*2, ncol=cols))

    # print row names
    for (v in 1:rows){
        rettab[v*2 - 1, 1] <- varlabels[varnames[v]]
        rettab[v*2, 1] <- c("")
    }

    # fill in numeric quantities
    for (m in 1:(cols-1)){
        for (v in 1:rows){
            if (varnames[v] %in% names(models[[m]]$coefficients)) {
            
                # find variance matrix -- named differently by ols(), tsls(), gmm()
                if(!is.null(models[[m]]$var)){
                    thevar <- models[[m]]$var
                } else if (!is.null(models[[m]]$V)){
                    thevar <- models[[m]]$V
                } else if (!is.null(models[[m]]$vcov)){
                    thevar <- models[[m]]$vcov
                }
            
                rettab[v*2 - 1, m+1] <- roundsig(models[[m]]$coefficients[varnames[v]], roundparam)
                t.ratio <- models[[m]]$coefficients[varnames[v]]/sqrt(thevar[varnames[v],varnames[v]])
                rettab[v*2, m+1] <- paste("{\\scriptsize (", roundsig(sqrt(thevar[varnames[v],varnames[v]]), roundparam), ")", "$^{", multiregtable_sigsymbol(t.ratio) ,"}$}", sep="")
            }
            else {
                rettab[v*2 - 1, m+1] <- c("")
                rettab[v*2, m+1] <- c("")
            }
        }
    }

    # tack on any extra rows handed to us
    if(!is.null(addrows)){
        startrow <- rows*2+1
        endrow <- rows*2+nrow(addrows)
        rettab[startrow:endrow,] <- addrows
        nadded <- nrow(addrows)
    } else {
        nadded <- 0
    }

    # add final summary stats
    startrow <- rows*2+nadded+1
    endrow <- rows*2+nadded+2
    rettab[startrow:endrow,] <- multiregtable_addons(models, roundparam)

    # return
    return(rettab)
}

multiregtable_sigsymbol <- function(t.ratio){

    if (pnorm(abs(t.ratio)) > 0.995){ return(c("***")) }
    else if (pnorm(abs(t.ratio)) > 0.975){ return(c("**")) }
    else if (pnorm(abs(t.ratio)) > 0.95){ return(c("*")) }
    else { return(c("")) }

}

multiregtable_addons <- function(models, roundparam){

    cols <- NROW(models) + 1
    rettab <- as.data.frame(matrix(NA, nrow=2, ncol=cols))
    rettab[1,1] <- c("N")
    rettab[2,1] <- c("$R^2$")

    for (m in 1:(cols-1)){
        if (models[[m]]$call[1] == "lm()"){
            rettab[1,m+1] <- NROW(models[[m]]$residuals)
            rettab[2,m+1] <- roundsig(summary.lm(models[[m]])$r.squared, roundparam)
        } else if (!is.null(models[[m]]$n)) {
            rettab[1,m+1] <- models[[m]]$n
            rettab[2,m+1] <- c("-")
        } else if (!is.null(models[[m]]$x)) {
            rettab[1,m+1] <- NROW(models[[m]]$x)
            rettab[2,m+1] <- roundsig(summary.lm(models[[m]])$r.squared, roundparam)
        } else {
            rettab[1,m+1] <- c("-")
            rettab[2,m+1] <- c("-")
        }
    }

    return(rettab)

}

roundsig <- function(vector, sigfigs){
    retvec <- rep("", NROW(vector))
    for (i in 1:NROW(vector)){
    
        number <- vector[i]
        prelim <- round(number, sigfigs)
        out <- prelim
        if ((sigfigs > 0) & (prelim %% 1 == 0)){out <- paste(out, ".", sep="")}
        addzeros <- ""
        az <- 1
        while ((az <= sigfigs) & ((prelim * 10^sigfigs) %% (10^az) == 0)){
            addzeros <- paste(addzeros, "0", sep="")
            az <- az + 1
        }
        retvec[i] <- paste(out, addzeros, sep="")
    }
    return(retvec)
}
