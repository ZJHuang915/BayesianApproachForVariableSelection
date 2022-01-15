#' This is the data set to be included in this package
#' @format  100 simulation data set with n=100 and p=10
"SimdataN100P10"


#' This is the data set to be included in this package
#' @format  100 simulation data set with n=50 and p=100
"SimdataN50P100"


#' AssignValue
#' @name AssignValue
#' @param x vector to be modified
#' @param idx index of elements to be modified
#' @param value values to be assigned
#' @return the vector after modified x
#' @export
AssignValue <- function(x,idx,value){
  x[idx] <- value ; return(x)
}


#' PriorProbability
#' @name PriorProbability
#' @param data.X predictors
#' @param data.Y responses
#' @return prior probability with calculate correlation between predictors and responses
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @export
PriorProbability <- function(data.X, data.Y){
  p.temp <- apply(data.X, 2, function(x){stats::cor(x,data.Y) %>% abs})
  p.prior <- p.temp/sum(p.temp)
  return(p.prior)
}


#' RidgeLambda
#' @name RidgeLambda
#' @param n.size sample size
#' @param p.size number of predictor
#' @param kt model size at time t
#' @return lambda value by formula 2.4
#' @export
RidgeLambda <- function(n.size, p.size, kt){
  if(n.size <= (p.size-1)){return(max(1/kt,1/300))}else{return(0)}
}


#' qBeta
#' @name qBeta
#' @param At1 random set, A at time t + 1
#' @param gt1 g-prior, g at time t + 1
#' @param sig2t1 \eqn{\sigma^2} at time t + 1
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of predictor
#' @return \eqn{\beta_{t+1}}
#' @importFrom mvtnorm rmvnorm
#' @export
qBeta <- function(At1, gt1, sig2t1, data.X, data.Y, n.size, p.size){
  kt1 <- length(At1)
  lambda <- RidgeLambda(n.size, p.size, kt=kt1)
  Ft1 <- (sig2t1^(-1) * t(data.X[,At1]) %*% data.X[,At1]) * (1+gt1^(-1)) + lambda * diag(1,kt1)
  mut1 <- sig2t1^(-1) * solve(Ft1) %*% t(data.X[,At1]) %*% data.Y
  betat1 <- AssignValue(x=rep(0,p.size), idx=At1, value=mvtnorm::rmvnorm(1, mean=mut1, sigma=solve(Ft1)))
  return(betat1)
}


#' probBeta
#' @name probBeta
#' @param betat1 \eqn{\beta} at time t + 1
#' @param At1 random set, A at time t + 1
#' @param gt1 g-prior, g at time t + 1
#' @param sig2t1 \eqn{\sigma^2} at time t + 1
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of predictor
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return proposal density of \eqn{\beta}
#' @importFrom mvtnorm dmvnorm
#' @export
probBeta <- function(betat1, At1, gt1, sig2t1, data.X, data.Y, n.size, p.size, log.p=F){
  kt1 <- length(At1)
  lamb <- RidgeLambda(n.size, p.size, kt=kt1)
  Ft1 <- (sig2t1^(-1) * t(data.X[,At1,drop=F]) %*% data.X[,At1,drop=F]) * (1+gt1^(-1)) +
    lamb * diag(1,kt1)
  mut1 <- sig2t1^(-1) * solve(Ft1) %*% t(data.X[,At1,drop=F]) %*% data.Y
  return(mvtnorm::dmvnorm(betat1[At1], mean=mut1, sigma=solve(Ft1), log=log.p))
}


#' qSigma2
#' @name qSigma2
#' @param sig2t \eqn{\sigma^2} at time t
#' @param eps_sig parameter for proposal function of \eqn{\sigma^2}
#' @return \eqn{\sigma^2_{t+1}}
#' @export
qSigma2 <- function(sig2t, eps_sig){
  sig2t1 <- stats::runif(1, min=max(10^(-8),sig2t-eps_sig), max=(sig2t+eps_sig))
  return(sig2t1)
}


#' probSigma2
#' @name probSigma2
#' @param sig2t1 \eqn{\sigma^2} at time t + 1
#' @param sig2t \eqn{\sigma^2} at time t
#' @param eps_sig parameter for proposal function of \eqn{\sigma^2}
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return proposal density of \eqn{\sigma^2}
#' @importFrom stats dunif
#' @export
probSigma2 <- function(sig2t1, sig2t, eps_sig, log.p=F){
  dunif(sig2t1, min=max(10^(-8),sig2t-eps_sig), max=(sig2t+eps_sig), log=log.p)
}


#' qGprior
#' @name qGprior
#' @param gt g-prior, g at time t
#' @param eps_g parameter for proposal function of g
#' @return \eqn{g_{t+1}}
#' @export
qGprior <- function(gt, eps_g){
  gt1 <- stats::runif(1, min=max(10^(-8),gt-eps_g), max=(gt+eps_g))
  return(gt1)
}


#' probGprior
#' @name probGprior
#' @param gt1 g-prior, g at time t + 1
#' @param gt g-prior, g at time t
#' @param eps_g parameter for proposal function of g
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return proposal density of g
#' @export
probGprior <- function(gt1, gt, eps_g, log.p=F){
  stats::dunif(gt1, min=max(10^(-8),gt-eps_g), max=(gt+eps_g), log=log.p)
}


#' qAlpha
#' @name qAlpha
#' @param alp element in A
#' @param At random set, A at time t
#' @param p.prio prior probability
#' @return \eqn{\alpha} for \eqn{A_{t}}
#' @export
qAlpha <- function(alp, At, p.prio){
  kt <- length(At)
  if(kt>1 & (alp%in%At)==F){
    return(p.prio[alp])
  }else if(kt>1 & (alp%in%At)==T){
    return(sum(p.prio[At])/(p.prio[alp] * sum(1/p.prio[At])))
  }else if(kt==1 & (alp%in%At)==T){
    return(0)
  }else if(kt==1 & (alp%in%At)==F){
    return(p.prio[alp]/sum(p.prio[-At]))
  }
}


#' qRandomsetA
#' @name qRandomsetA
#' @param p.size number of parameter
#' @param At random set, A at time t
#' @param ch \eqn{c_h} indicator of model change
#' @param p.prio prior probability
#' @return \eqn{A_{t+1}}
#' @export
qRandomsetA <- function(p.size, At, ch, p.prio){
  A.alpha.i <- sapply(1:p.size, qAlpha, At=At, p.prio=p.prio)
  A.alpha <- sample(1:p.size, size=1, prob=A.alpha.i)
  if(ch==0){
    At1 <- At
  }else{
    ifelse((A.alpha %in% At), At1 <- subset(At, At!=A.alpha), At1 <- c(At,A.alpha))
  }
  return(At1)
}


#' probRandomsetA
#' @name probRandomsetA
#' @param At1 random set, A at time t + 1
#' @param At random set, A at time t
#' @param ch \eqn{c_h} indicator of model change
#' @param p.prio prior probability
#' @return proposal density of A
#' @importFrom magrittr %>%
#' @export
probRandomsetA <- function(At1, At, ch, p.prio){
  lngAt1 <- length(At1)
  lngAt <- length(At)
  if(ch==0){
    ifelse((lngAt1==lngAt & sum(At1 %in% At)==lngAt), 1, 0)
  }else{
    if(lngAt1 > lngAt){
      At1[!(At1 %in% At)] %>% qAlpha(At=At, p.prio=p.prio)
    }else{
      At[!(At %in% At1)] %>% qAlpha(At=At, p.prio=p.prio)
    }
  }
}


#' fconditionGprior
#' @name fconditionGprior
#' @param g g-prior
#' @param beta \eqn{\beta}
#' @param sig2 \eqn{\sigma^2}
#' @param A random set, A
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of parameter
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return full conditional density of g
#' @export
fconditionGprior <- function(g, beta, sig2, A, data.X, data.Y, n.size, p.size, log.p=F){
  k <- length(A)
  lamb <- RidgeLambda(n.size=n.size, p.size=p.size, kt=k)
  Cov <- g^(-1) * sig2^(-1) * t(data.X[,A,drop=F]) %*% data.X[,A,drop=F] + lamb * diag(1,k)
  p.g <- det(Cov)^(1/2) * exp(-(t(beta[A]) %*% Cov %*% beta[A])/2) *
    g^(-1/2-1) * exp(-n.size/(2*g))
  if(log.p==T){
    p.g <- (1/2)*log(det(Cov)) -(t(beta[A]) %*% Cov %*% beta[A])/2 -3/2*log(g) -n.size/(2*g)
  }
  return(p.g)
}


#' fconditionSigma2
#' @name fconditionSigma2
#' @param sig2 \eqn{\sigma^2}
#' @param beta \eqn{\beta}
#' @param g g-prior
#' @param A random set, A
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of parameter
#' @param sig.a parameter a for full conditional of \eqn{\sigma^2}
#' @param sig.b parameter b for full conditional of \eqn{\sigma^2}
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return full conditional density of \eqn{\sigma^2}
#' @export
fconditionSigma2 <- function(sig2, beta, g, A, data.X, data.Y, n.size, p.size,
                             sig.a, sig.b, log.p=F){
  k <- length(A)
  lamb <- RidgeLambda(n.size=n.size, p.size=p.size, kt=k)
  Cov <- g^(-1) * sig2^(-1) * t(data.X[,A,drop=F]) %*% data.X[,A,drop=F] + lamb * diag(1,k)
  p.sig2 <- sig2^(-n.size/2) * exp(-(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
                                       (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2)) *
    det(Cov)^(1/2) * exp(-(t(beta[A]) %*% Cov %*% beta[A])/2) *
    sig2^(-sig.a-1) * exp(-sig.b/sig2)
  if(log.p==T){
    p.sig2 <- (-n.size/2)*log(sig2) -(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
                                        (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2) +
      (1/2)*log(det(Cov)) -(t(beta[A]) %*% Cov %*% beta[A])/2 -(sig.a+1)*log(sig2) -sig.b/sig2
  }
  return(p.sig2)
}


#' fconditionBeta
#' @name fconditionBeta
#' @param sig2 \eqn{\sigma^2}
#' @param beta \eqn{\beta}
#' @param g g-prior
#' @param A random set, A
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of parameter
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return full conditional density of \eqn{\beta}
#' @export
fconditionBeta <- function(beta, sig2, g, A, data.X, data.Y, n.size, p.size, log.p=F){
  k <- length(A)
  lamb <- RidgeLambda(n.size, p.size, kt=k)
  Cov <- g^(-1) * sig2^(-1) * t(data.X[,A,drop=F]) %*% data.X[,A,drop=F] + lamb * diag(1,k)
  p.beta <- exp(-(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
                    (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2)) *
    exp(-(t(beta[A]) %*% Cov %*% beta[A])/2)
  if(log.p==T){
    p.beta <- -(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
                  (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2) -(t(beta[A]) %*% Cov %*% beta[A])/2
  }
  return(p.beta)
}


#' fconditionRandomset
#' @name fconditionRandomset
#' @param sig2 \eqn{\sigma^2}
#' @param beta \eqn{\beta}
#' @param g g-prior
#' @param A random set, A
#' @param data.X predictors
#' @param data.Y responses
#' @param n.size sample size
#' @param p.size number of parameter
#' @param p.prio prior probability
#' @param power.l power transformation of truncated binomial
#' @param log.p logical. If TRUE, probability p is given as log(p)
#' @return full conditional density of A
#' @importFrom extraDistr dtbinom
#' @export
fconditionRandomset <- function(A, beta, sig2, g, data.X, data.Y, p.prio,
                                n.size, p.size, power.l, log.p=F){
  k <- length(A)
  lamb <- RidgeLambda(n.size, p.size, kt=k)
  Cov <- g^(-1) * sig2^(-1) * t(data.X[,A,drop=F]) %*% data.X[,A,drop=F] + lamb * diag(1,k)
  p.A <- exp(-(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
                 (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2)) *
    det(Cov)^(1/2) * exp(-(t(beta[A]) %*% Cov %*% beta[A])/2) *
    sum(p.prio[A]) * extraDistr::dtbinom(k, size=p.size, prob=2.5/p.size, a=0)^power.l / k
  if(log.p==T){
    p.A <- -(t(data.Y-data.X[,A,drop=F] %*% beta[A]) %*%
               (data.Y-data.X[,A,drop=F] %*% beta[A]))/(2*sig2) + (1/2)*log(det(Cov)) -
      (t(beta[A]) %*% Cov %*% beta[A])/2 + log(sum(p.prio[A])) +
      power.l*log(extraDistr::dtbinom(k, size=p.size, prob=2.5/p.size, a=0)) - log(k)
  }
  return(p.A)
}

#' MHAlgorithm
#' @name MHAlgorithm
#' @param sim.size simulation sample size
#' @param data observed data
#' @param p.prio prior of p can be specified
#' @param tune.p_h tuning parameter
#' @param tune.eps_sig tuning parameter
#' @param tune.eps_g tuning parameter
#' @param tune.mu_truncbinom tuning parameter
#' @param tune.sig_a tuning parameter
#' @param tune.sig_b tuning parameter
#' @param ini.g initial value g
#' @param ini.sig2 initial value \eqn{\sigma^2}
#' @param ini.A initial random set, A
#' @param ini.beta initial \eqn{\beta}
#' @return simulation samples and acceptance rate
#' @examples
#' \dontrun{
#' MHAlgorithm(sim.size = 50000, data = cbind(diabetes.data$y,diabetes.data$x),
#' tune.p_h = 0.5, tune.eps_sig = 0.1, tune.eps_g = 60,
#' tune.mu_truncbinom = 2.5, tune.sig_a = 0.001, tune.sig_b = 0.001)
#' }
#' @export
MHAlgorithm <- function(sim.size, data, p.prio=NULL, tune.p_h, tune.eps_sig, tune.eps_g,
                        tune.mu_truncbinom, tune.sig_a, tune.sig_b,
                        ini.g=NULL, ini.sig2=NULL, ini.A=NULL, ini.beta=NULL){
  data.Y <- data[,1]
  data.X <- data[,-1]
  if(is.null(p.prio)){p.prio <- PriorProbability(data.X, data.Y)}
  n.size <- nrow(data.X)
  p.size <- ncol(data.X)
  power.l <- ifelse(p.size < 30, 1, 3)

  if(is.null(ini.g)){ini.g <- 1}
  if(is.null(ini.sig2)){ini.sig2 <- 1}
  if(is.null(ini.A)){ini.A <- sample(1:p.size, ceiling(tune.mu_truncbinom))}
  if(is.null(ini.beta)){
    ini.beta <- qBeta(At1=ini.A, gt1=ini.g, sig2t1=ini.sig2, data.X, data.Y,
                      n.size, p.size)
  }

  g.sample <- c()
  sig2.sample <- c()
  A.sample <- list()
  A.sample_ch <- c()  # length < length(A.sample_ch) - 1
  A.sample_k <- c()
  beta.sample <- list()
  sim.sample <- list()

  g.sample[1] <- ini.g
  sig2.sample[1] <- ini.sig2
  A.sample[[1]] <- ini.A
  A.sample_k[1] <- length(A.sample[[1]])
  beta.sample[[1]] <- ini.beta

  acceptance.time <- c(acpt.g=0, acpt.sig2=0, acpt.A=0, acpt.beta=0)

  for(i in 2:sim.size){
    # draw g
    g.cand <- qGprior(gt=g.sample[i-1], eps_g=tune.eps_g)
    alp <- min(1,exp(fconditionGprior(g=g.cand, beta=beta.sample[[i-1]], sig2=sig2.sample[i-1],
                                      A=A.sample[[i-1]], data.X, data.Y, n.size, p.size, log.p=T) -
                       fconditionGprior(g=g.sample[i-1], beta=beta.sample[[i-1]], sig2=sig2.sample[i-1],
                                        A=A.sample[[i-1]], data.X, data.Y, n.size, p.size, log.p=T) +
                       probGprior(gt1=g.sample[i-1], gt=g.cand, eps_g=tune.eps_g, log.p=T) -
                       probGprior(gt1=g.cand, gt=g.sample[i-1], eps_g=tune.eps_g, log.p=T)))
    if(stats::runif(1) <= alp){
      g.sample[i] <- g.cand
      acceptance.time[1] <- acceptance.time[1] + 1
    }else{
      g.sample[i] <- g.sample[i-1]
    }

    # draw sigma2
    sig2.cand <- qSigma2(sig2t=sig2.sample[i-1], eps_sig=tune.eps_sig)
    alp <- min(1,exp(fconditionSigma2(sig2=sig2.cand, beta=beta.sample[[i-1]], g=g.sample[i],
                                      A=A.sample[[i-1]], data.X, data.Y, n.size, p.size,
                                      sig.a=tune.sig_a, sig.b=tune.sig_b, log.p=T) -
                       fconditionSigma2(sig2=sig2.sample[i-1], beta=beta.sample[[i-1]], g=g.sample[i],
                                        A=A.sample[[i-1]], data.X, data.Y, n.size, p.size,
                                        sig.a=tune.sig_a, sig.b=tune.sig_b, log.p=T) +
                       probSigma2(sig2t1=sig2.sample[i-1], sig2t=sig2.cand, eps_sig=tune.eps_sig,
                                  log.p=T) -
                       probSigma2(sig2t1=sig2.cand, sig2t=sig2.sample[i-1], eps_sig=tune.eps_sig,
                                  log.p=T)))
    if(stats::runif(1) <= alp){
      sig2.sample[i] <- sig2.cand
      acceptance.time[2] <- acceptance.time[2] + 1
    }else{
      sig2.sample[i] <- sig2.sample[i-1]
    }

    # draw A
    A.sample_ch[i-1] <- stats::rbinom(1, size=1, prob=tune.p_h)
    A.cand <- qRandomsetA(p.size=p.size, At=A.sample[[i-1]], ch=A.sample_ch[i-1], p.prio=p.prio)
    if(A.sample_ch[i-1]==0){
      alp <- 0
    }else{
      alp <- min(1,exp(fconditionRandomset(A=A.cand, beta=beta.sample[[i-1]],sig2=sig2.sample[i],
                                           g=g.sample[i], data.X, data.Y, p.prio, n.size,
                                           p.size, power.l, log.p=T) -
                         fconditionRandomset(A=A.sample[[i-1]], beta=beta.sample[[i-1]],sig2=sig2.sample[i],
                                             g=g.sample[i], data.X, data.Y, p.prio, n.size,
                                             p.size, power.l, log.p=T)) *
                   probRandomsetA(At1=A.sample[[i-1]], At=A.cand, ch=A.sample_ch[i-1], p.prio) /
                   probRandomsetA(At1=A.cand, At=A.sample[[i-1]], ch=A.sample_ch[i-1], p.prio))
    }
    if(stats::runif(1) <= alp){
      A.sample[[i]] <- A.cand
      acceptance.time[3] <- acceptance.time[3] + 1
    }else{
      A.sample[[i]] <- A.sample[[i-1]]
    }
    A.sample_k[i] <- length(A.sample[[i]])

    # draw beta
    beta.cand <- qBeta(At1=A.sample[[i]], gt1=g.sample[i], sig2t1=sig2.sample[i],
                       data.X, data.Y, n.size, p.size)
    alp <- min(1,exp(fconditionBeta(beta=beta.cand, sig2=sig2.sample[i], g=g.sample[i],
                                    A=A.sample[[i]], data.X, data.Y, n.size, p.size, log.p=T) -
                       fconditionBeta(beta=beta.sample[[i-1]], sig2=sig2.sample[i], g=g.sample[i],
                                      A=A.sample[[i]], data.X, data.Y, n.size, p.size, log.p=T) +
                       probBeta(betat1=beta.sample[[i-1]], At1=A.sample[[i]],
                                gt1=g.sample[i], sig2t1=sig2.sample[i],
                                data.X, data.Y, n.size, p.size, log.p=T) -
                       probBeta(betat1=beta.cand, At1=A.sample[[i]],
                                gt1=g.sample[i], sig2t1=sig2.sample[i],
                                data.X, data.Y, n.size, p.size, log.p=T)))
    if(stats::runif(1) <= alp){
      beta.sample[[i]] <- beta.cand
      acceptance.time[4] <- acceptance.time[4] + 1
    }else{
      beta.sample[[i]] <- beta.sample[[i-1]]
    }
  }

  sim.sample <- list(beta.sample, g.sample, sig2.sample,
                     A.sample)
  acceptance.rate <- acceptance.time / (sim.size-1)

  return(list(sim.sample=sim.sample, acceptance.rate=acceptance.rate))
}
