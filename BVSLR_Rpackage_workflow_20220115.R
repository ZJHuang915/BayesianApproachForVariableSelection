library(devtools)
library(magrittr)
library(roxygen2)
##### Bayesian approach variable selection in linear regression
##### BAVSLR
setwd("C:\\Users\\88692\\Desktop\\")
create("BVSLR")


##### data construction #####
### simdata_n50_p100: "sim.beta"; "sim.X_n50"; "sim.Y_n50"; "sim.cov"
### simdata_n100_p10: "sim.beta"; "sim.X_n100"; "sim.Y_n100"; "sim.cov"

# function for assigning value 
AssignValue <- function(x,indx,value){x[indx] <- value ; return(x)}

# number of predictor
sim1_P <- 10
sim2_P <- 100

# covariance matrix
sim1_Cov <- matrix(0.6, nrow=sim1_P, ncol=sim1_P) %>% `diag<-`(1)
sim2_Cov <- matrix(0.6, nrow=sim2_P, ncol=sim2_P) %>% `diag<-`(1)

# true betas 
sim1_beta <- rep(0, sim1_P) %>% AssignValue(x=.,
                                            indx=c(1,2,3,4),
                                            value=c(-2.5,-1.5,1.5,2.5))
sim2_beta <- rep(0, sim2_P) %>% AssignValue(x=.,
                                            indx=c(1,11,21,51,71,81),
                                            value=c(-2.5,-2,-1.5,1.5,2,2.5))

# generating X with n=50, n=100
set.seed(3072)
sim1_X <- lapply(1:100,function(i) mvtnorm::rmvnorm(n=100, mean=rep(0,sim1_P), sigma=sim1_Cov)) # 100 groups: 100x10(nxp)
sim2_X <- lapply(1:100,function(i) mvtnorm::rmvnorm(n=50, mean=rep(0,sim2_P), sigma=sim2_Cov)) # 100 groups: 50x100(nxp)

# generating Y with n=50, n=100
set.seed(3072)
sim1_Y <- lapply(1:100,function(i) 
  sim1_X[[i]] %*% sim1_beta + rnorm(n=100, mean=0, sd=1)) # 100 groups: 100x1
sim2_Y <- lapply(1:100,function(i) 
  sim2_X[[i]] %*% sim2_beta + rnorm(n=50, mean=0, sd=1)) # 100 groups: 50x1


SimdataN100P10 <- list(beta=sim1_beta, Xdata=sim1_X, Ydata=sim1_Y, Cov=sim1_Cov)
SimdataN50P100 <- list(beta=sim2_beta, Xdata=sim2_X, Ydata=sim2_Y, Cov=sim2_Cov)

rm('AssignValue')
setwd("C:\\Users\\88692\\Desktop\\BVSLR")

usethis::use_data(SimdataN100P10, SimdataN50P100, internal=F,overwrite = TRUE)
usethis::use_package("magrittr")
usethis::use_package("mvtnorm")
usethis::use_package("extraDistr")
usethis::use_package("stats")
#  DESCRIPTION file

document()
check()
build()

#### terminal running
# cd C:/Users/88692/Desktop
# R CMD build BVSLR
# R CMD check BVSLR_1.0.tar.gz
# R CMD check BVSLR_1.0.tar.gz --as-cran
