two.bt <- function (sample1, sample2, FUN, R, ratio=FALSE, student = FALSE, M, weights = NULL, 
    ...) 
{
    func.name <- ifelse(is.character(FUN), FUN, deparse(substitute(FUN)))
    func <- match.fun(FUN)
    ind <- c(rep(1, length(sample1)), rep(2, length(sample2)))
    nobsgrp <- as.numeric(table(ind))
    extra <- list(...)
    if (func.name == "quantile") {
        if (is.na(match("probs", names(extra)))) 
            stop("'probs' argument must be specified")
        if (length(extra$probs) > 1) 
            stop("can only bootstrap a single quantile")
    }
    boot.func <- function(x, idx) {
        d1 <- x[idx[ind == 1]]
        d2 <- x[idx[ind == 2]]
        if(ratio) fval <- func(d1, ...) / func(d2, ...) else fval <- func(d1, ...) - func(d2, ...)
        if (student) {
            b <- two.boot(d1, d2, FUN, R = M, student = FALSE, 
                M = NULL, weights = NULL, ...)
            fval <- c(fval, var(b$t))
        }
        fval
    }
    if (!is.null(weights)) 
        weights <- unlist(weights)
    b <- boot(c(sample1, sample2), statistic = boot.func, R = R, 
        weights = weights, strata = ind)
    b$student <- student
    b$student <- student
    structure(b, class = "simpleboot")
}

sum.bt <- function(b) unlist(list(Observed=b$t0,SE=sd(b$t),Bias=mean(b$t)-b$t0))

plot.bt <- function(b){
  par(mfrow=c(1,2))
  nc <- min(ceiling(length(b$t)/25),100)
  nc <- max(nc,10)
  h <- hist(b$t,breaks=nc,prob=TRUE,main="Histogram of t")
  mx <- max(h$density) 
  lines(c(b$t0,b$t0),c(-mx/5,1.1*mx),lty=2,col='red')
  qqnorm(b$t,xlab="Quantiles of standard normal",cex=.7)
  qqline(b$t,lty=2)
  par(mfrow=c(1,1))
}

btt.ci <- function(boot.out,conf=.95,df=NULL) {
  # computes bootstrap t CI's
  # can either specify df (df=0 will use z) or default is to
  # use n-1 for one.boot and min sample size -1 for two.boot
  if (!inherits(boot.out, "simpleboot")) 
    stop("only use this function on 'simpleboot' objects")
  if (conf.level < 0 || conf.level > 1) 
    stop("conf.level must be between 0 and 1")
  if(!is.null(df)) if(df < 0) stop("df must be > 0")
  if (!is.null(boot.out$student) && boot.out$student) 
    x <- boot.out$t[, 1, drop = FALSE]
  else x <- boot.out$t
  m <- qnorm((1+conf)/2)
  if(is.null(df)) {
    m <- qt((1+conf.level)/2,min(table(boot.out$strata))-1)
  } else if(df!=0) m <- qt((1+conf)/2,df)
  out <- c(boot.out$t0,boot.out$t0 + c(-1,1)*m*sd(x))
  names(out) <- c("estimate","lower","upper")
  out
}