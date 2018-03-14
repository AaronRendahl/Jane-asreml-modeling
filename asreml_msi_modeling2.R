library(asreml)
library(aod)

## get coefficients and p-values
getcoefp <- function(model) {
  x <- as.data.frame(summary(model, all=TRUE)$coef.fi)
  #x <- x[rev(seq_len(nrow(x))),]
  pval <- as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
  out <- cbind(x, pval)
  out[!is.na(out[,2]),]
}

## overall Wald chi square test of a given effect on entire trajectory
WaldChi2 <- function(model, var) {
  covariance_matrix <- diag(length(model$coefficients$fixed))
  covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)] <- model$Cfixed
  covariance_matrix[lower.tri(covariance_matrix)] <- t(covariance_matrix)[lower.tri(t(covariance_matrix))]
  colnames(covariance_matrix) <- rownames(covariance_matrix) <- names(model$coefficients$fixed)
  use <- grep(var, rownames(covariance_matrix))
  terms <- use[rowSums(covariance_matrix[use,])!=0]
  out <- wald.test(b = model$coefficients$fixed[terms], 
            Sigma = covariance_matrix[terms, terms], 
            Terms = 1:length(terms))$result$chi2
  out <- rbind(out)
  rownames(out) <- var
  out
}

## get fitted values for a given model
getfit <- function(model, var, data=OST_R_asreml_NS_hr) {
  if(missing(var)) {
    predict(model, classify="time")$predictions$pvals
  } else {
    if(is.numeric(data[[var]])) {
      m <- mean(data[[var]])
      s <- sd(data[[var]])
      levs <- setNames(list(c(m, m+s)), var)
      predict(model, classify=paste0("time:", var), levels=levs)$predictions$pvals
    } else {
      predict(model, classify=paste0("time:", var))$predictions$pvals
    }
  }
}


#########################################################################################
message("Using ", R.version.string, " and asreml ", packageVersion("asreml"), ".")

OST_R_asreml_NS_hr <- read.csv("OST_R_asreml_NS_hr.csv", check.names=FALSE)
OST_R_asreml_NS_hr$Number_ID <- factor(OST_R_asreml_NS_hr$Number_ID)
OST_R_asreml_NS_hr$logINS <- log(OST_R_asreml_NS_hr$insulin)
OST_R_asreml_NS_hr$time <- OST_R_asreml_NS_hr$time*60
OST_R_asreml_NS_hr$time2 <- NULL
OST_R_asreml_NS_hr$EMS <- substring(OST_R_asreml_NS_hr$EMS, 1, 1)

## these are the covariates we'll add
vars <- c("SEX", "AGE", "Breed", "Adiponectin", "TG", "NEFA", "Leptin", "TBFM", "SI", "AIRg", "EMS")

## primary models
glu <- asreml(glucose ~ pol(time,-2), data=OST_R_asreml_NS_hr, 
               random=~str(~Number_ID:pol(time,2), ~us(3):id(Number_ID)), 
               rcov=~units, maxit=400, Cfixed=TRUE)

ins <- asreml(logINS ~ pol(time,-2), data=OST_R_asreml_NS_hr, 
              random=~str(~Number_ID:pol(time,2), ~us(3):id(Number_ID)), 
              rcov=~units, maxit=400, Cfixed=TRUE)

## get coefficients and p-values
glu1 <- list(model=glu, coef=getcoefp(glu), fit=getfit(glu))
ins1 <- list(model=ins, coef=getcoefp(ins), fit=getfit(ins))

## add covariates one by one, and save all results
glu2 <- lapply(vars, function(v) {
  mm <- update(glu, as.formula(sprintf("~ . + %s * pol(time, 2)", v)))
  list(model=mm, coef=getcoefp(mm), Wald=WaldChi2(mm, v), fit=getfit(mm, v))
})
ins2 <- lapply(vars, function(v) {
  mm <- update(ins, as.formula(sprintf("~ . + %s * pol(time, 2)", v)))
  list(model=mm, coef=getcoefp(mm), Wald=WaldChi2(mm, v), fit=getfit(mm, v))
})

