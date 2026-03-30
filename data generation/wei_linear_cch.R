## Prerequisites ##

library(MASS)

## Simulation setting ##

path <- "/nas/longleaf/home/haolin/dissertation2/round4/data/"

gamma.weibull = 0.7
p = 10
p.full = 2
n = 1500
q = 0.1
nsim = 500
corr = 0.6
beta = matrix(c(rep(0,p.full-2), c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), rep(0,p-p.full-8)), nrow=p)

cenc = rep(0, nsim)

## Data generation (training) ##

for (j in 1:nsim){
  cat(j)
  
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  X = as.matrix(mvrnorm(n = 2*n, matrix(0, nrow = p, ncol = 1), ar1_cor(p, corr)))
  
  true.lin.pred <- X %*% beta
  v <- runif(n=2*n)
  Tlat <- (- log(v) / (1* exp(true.lin.pred)))^(1 / gamma.weibull)
 
  C = runif(2*n, min = 0, max = 0.03) # censoring times

  end <- rep(0.0162, 2*n)
  time <- pmin(Tlat, C, end) # follow-up times
  
  status <- as.numeric(Tlat == time) #event indicators
  
  dat.pre = data.frame(X, status, time)
  
  dat = dat.pre[1:n,]
  test = dat.pre[(n+1):(2*n),]

  dat$subid = 1:n
  test$subid = 1:n
  
  dat$sc = sample(c(1:0), n, replace=T, prob = c(q,1-q))
  q.comp = 0.145
  dat$r = sample(c(1:0), n, replace=T, prob = c(q.comp,1-q.comp))
 
  cenc[j] = table(dat$status)[1]/n
  
  ntilde = sum(dat$sc)
  
  cc <- subset(dat, select = -c(r))
  for (m in 1:n){
    if ((cc$sc[m]==1)|(dat$status[m]==1)){
      cc[m,c((p.full+1):p)] = cc[m,c((p.full+1):p)]
    }else{
      cc[m,c((p.full+1):p)] = NA
    }
  }
  cc$weight = 1*(cc$status==1)+n/ntilde*(cc$status==0)
  write.csv(cc,file=paste0(path, 'train', j, '.csv'), row.names = F)
  
  write.csv(test,file=paste0(path, 'test', j, '.csv'), row.names = F)
  
}

mean(cenc)

