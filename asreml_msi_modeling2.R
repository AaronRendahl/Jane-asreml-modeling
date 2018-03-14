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
      m <- mean(data[[var]], na.rm=TRUE)
      s <- sd(data[[var]], na.rm=TRUE)
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
OST_R_asreml_NS_hr$TBFM <- as.numeric(as.character(OST_R_asreml_NS_hr$TBFM))

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
  mm <- update(glu, as.formula(sprintf("~ . + %s * pol(time, 2)", v)), 
               data=OST_R_asreml_NS_hr[!is.na(OST_R_asreml_NS_hr[[v]])])
  list(model=mm, coef=getcoefp(mm), Wald=WaldChi2(mm, v), fit=getfit(mm, v))
})
ins2 <- lapply(vars, function(v) {
  mm <- update(ins, as.formula(sprintf("~ . + %s * pol(time, 2)", v)),
               data=OST_R_asreml_NS_hr[!is.na(OST_R_asreml_NS_hr[[v]])])
  list(model=mm, coef=getcoefp(mm), Wald=WaldChi2(mm, v), fit=getfit(mm, v))
})

## combine together for output
g <- setNames(c(list(glu1), glu2), c("time_only", vars))
i <- setNames(c(list(ins1), ins2), c("time_only", vars))
saveRDS(g, "results_glucose.RDS")
saveRDS(i, "results_insulin.RDS")

## and make some csv files
g.wald <- do.call(rbind, lapply(g, function(x) data.frame(var=rownames(x$Wald), x$Wald, stringsAsFactors=FALSE)))
i.wald <- do.call(rbind, lapply(i, function(x) data.frame(var=rownames(x$Wald), x$Wald, stringsAsFactors=FALSE)))

g.coef <- do.call(rbind, lapply(names(g), function(n) {
  out <- cbind(var=n, term=rownames(g[[n]]$coef), as.data.frame(g[[n]]$coef), stringsAsFactors=FALSE)
  rownames(out) <- NULL
  out
}))

i.coef <- do.call(rbind, lapply(names(i), function(n) {
  out <- cbind(var=n, term=rownames(i[[n]]$coef), as.data.frame(i[[n]]$coef), stringsAsFactors=FALSE)
  rownames(out) <- NULL
  out
}))

g.fit <- do.call(rbind, lapply(names(g), function(n) {
  tmp <- as.data.frame(g[[n]]$fit)
  if(ncol(tmp)==4) { tmp <- data.frame(tmp[,1,drop=FALSE], X=NA, tmp[,2:4]) }
  stopifnot(ncol(tmp)==5)
  names(tmp)[2] <- "level"
  cbind(var=n, tmp, stringsAsFactors=FALSE)
}))

i.fit <- do.call(rbind, lapply(names(i), function(n) {
  tmp <- as.data.frame(i[[n]]$fit)
  if(ncol(tmp)==4) { tmp <- data.frame(tmp[,1,drop=FALSE], X=NA, tmp[,2:4]) }
  stopifnot(ncol(tmp)==5)
  names(tmp)[2] <- "level"
  cbind(var=n, tmp, stringsAsFactors=FALSE)
}))

## combine together glucose and insulin into one file
cc <- function(x) {
  out <- do.call(rbind, lapply(names(x), function(n) {
    cbind(response=n, x[[n]], stringsAsFactors=FALSE)
  }))
  rownames(out) <- NULL
  out
}

all.wald <- cc(list(glucose=g.wald, insulin=i.wald))
all.coef <- cc(list(glucose=g.coef, insulin=i.coef))
all.fit <- cc(list(glucose=g.fit, insulin=i.fit))

write.csv(all.wald, "results-wald.csv", na="", row.names=FALSE)
write.csv(all.coef, "results-coef.csv", na="", row.names=FALSE)
write.csv(all.fit, "results-fit.csv", na="", row.names=FALSE)