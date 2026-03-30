#path.test <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/"
#path.cc <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/"

path.results <- "/nas/longleaf/home/haolin/dissertation2/G_round4/results/"
#path.c <- "/nas/longleaf/home/haolin/dissertation1/round2/data/c/"
path.cc <- "/nas/longleaf/home/haolin/dissertation2/G_round4/data/"

library(rpart)
library(partykit)
library(survival)
library(randomForestSRC)
library(tidyr)
library(intsurv)

nsim = 100

tau = 0.0162
#p=4
#true = c(0.1, 0.3, -0.3, -0.1)
p=10
p.full = 2
gamma.weibull = 0.7
beta = matrix(c(rep(0,p.full-2), c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), rep(0,p-p.full-8)), nrow=p)

rf.cc.pa = matrix(NA, nrow = 1, ncol = nsim)
C.cc.pa = matrix(NA, nrow = 1, ncol = nsim)
rf.cc.time = matrix(NA, nrow = 1, ncol = nsim)



#m=2
for (m in 1:nsim){
  try({
  cat(m)
  
  # read datasets
  cc = read.csv(paste0(path.cc, 'train', m, '.csv'))
  cc = subset(cc, select = -sc)
  test = read.csv(paste0(path.cc, 'test', m, '.csv'))
  
  # true survival
  dat.pa = test
  
  time.pa = data.frame(time = sort(dat.pa$time), ind = 1)
  
  true.lin.pred <- as.matrix(dat.pa[,1:p]) %*% beta +1.5
  #true.lin.pred <- 0.3*as.matrix(dat.pa[,1]) - 0.2*as.matrix(dat.pa[,2]) +0.1*as.matrix(dat.pa[,3]) - 0.9*as.matrix(dat.pa[,4]) -0.3 + as.matrix(dat.pa[,3])*as.matrix(dat.pa[,2]) - 1.6*as.matrix(dat.pa[,4])*as.matrix(dat.pa[,4])
  survival.true <- function(t){
    survival =  exp(-exp(true.lin.pred)*t^gamma.weibull)
    return(survival)
  }
  true.survival = data.frame(subid = c(1:nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                             t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                             t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                             t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                             t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  true.survival$t1 = survival.true(1*tau/15)
  true.survival$t2 = survival.true(2*tau/15)
  true.survival$t3 = survival.true(3*tau/15)
  true.survival$t4 = survival.true(4*tau/15)
  true.survival$t5 = survival.true(5*tau/15)
  true.survival$t6 = survival.true(6*tau/15)
  true.survival$t7 = survival.true(7*tau/15)
  true.survival$t8 = survival.true(8*tau/15)
  true.survival$t9 = survival.true(9*tau/15)
  true.survival$t10 = survival.true(10*tau/15)
  true.survival$t11 = survival.true(11*tau/15)
  true.survival$t12 = survival.true(12*tau/15)
  true.survival$t13 = survival.true(13*tau/15)
  true.survival$t14 = survival.true(14*tau/15)
  
  # prediction
  srf.survival.pred = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                                 t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                                 t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                                 t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                                 t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
  cc = cc[is.na(cc[,p])==0,]
  weight = cc$weight
  cc = subset(cc, select = -c(weight, subid))
  
  start_time <- Sys.time()
  number = tune(Surv(time, status) ~ ., case.wt = weight, data = cc, ntree=500)
  srf.fit = rfsrc(Surv(time, status) ~ ., case.wt = weight, data = cc, ntree=500,mtry=number$optimal[1], nodesize=number$optimal[2])
  
  CHF.test = matrix(NA, nrow = nrow(dat.pa), ncol = 1)
  for (i in 1:nrow(dat.pa)){
    srf.pred = data.frame(predict(srf.fit, newdata = dat.pa[i,1:p])$survival[1,], predict(srf.fit, newdata = dat.pa[i,1:p])$time.interest)
    chf.pred = data.frame(chf=predict(srf.fit, newdata = dat.pa[i,1:p])$chf[1,], time=predict(srf.fit, newdata = dat.pa[i,1:p])$time.interest)
    chf.test.temp = merge(chf.pred, time.pa, by="time", all=T)
    chf.test.temp = chf.test.temp %>% fill(chf)
    chf.test.temp = chf.test.temp[is.na(chf.test.temp$ind)==F,]
    CHF.test[i] = sum(chf.test.temp$chf, na.rm = T)
    subid = i
    srf.pred$time1 = as.numeric(srf.pred[,2]-1*tau/15 <0)
    srf.pred$time2 = as.numeric(srf.pred[,2]-2*tau/15 <0)
    srf.pred$time3 = as.numeric(srf.pred[,2]-3*tau/15 <0)
    srf.pred$time4 = as.numeric(srf.pred[,2]-4*tau/15 <0)
    srf.pred$time5 = as.numeric(srf.pred[,2]-5*tau/15 <0)
    srf.pred$time6 = as.numeric(srf.pred[,2]-6*tau/15 <0)
    srf.pred$time7 = as.numeric(srf.pred[,2]-7*tau/15 <0)
    srf.pred$time8 = as.numeric(srf.pred[,2]-8*tau/15 <0)
    srf.pred$time9 = as.numeric(srf.pred[,2]-9*tau/15 <0)
    srf.pred$time10 = as.numeric(srf.pred[,2]-10*tau/15 <0)
    srf.pred$time11 = as.numeric(srf.pred[,2]-11*tau/15 <0)
    srf.pred$time12 = as.numeric(srf.pred[,2]-12*tau/15 <0)
    srf.pred$time13 = as.numeric(srf.pred[,2]-13*tau/15 <0)
    srf.pred$time14 = as.numeric(srf.pred[,2]-14*tau/15 <0)
    srf.survival.pred$subid[subid] = subid
    srf.survival.pred$t1[subid] = min(srf.pred[srf.pred$time1>0, 1])
    srf.survival.pred$t2[subid] = min(srf.pred[srf.pred$time2>0, 1])
    srf.survival.pred$t3[subid] = min(srf.pred[srf.pred$time3>0, 1])
    srf.survival.pred$t4[subid] = min(srf.pred[srf.pred$time4>0, 1])
    srf.survival.pred$t5[subid] = min(srf.pred[srf.pred$time5>0, 1])
    srf.survival.pred$t6[subid] = min(srf.pred[srf.pred$time6>0, 1])
    srf.survival.pred$t7[subid] = min(srf.pred[srf.pred$time7>0, 1])
    srf.survival.pred$t8[subid] = min(srf.pred[srf.pred$time8>0, 1])
    srf.survival.pred$t9[subid] = min(srf.pred[srf.pred$time9>0, 1])
    srf.survival.pred$t10[subid] = min(srf.pred[srf.pred$time10>0, 1])
    srf.survival.pred$t11[subid] = min(srf.pred[srf.pred$time11>0, 1])
    srf.survival.pred$t12[subid] = min(srf.pred[srf.pred$time12>0, 1])
    srf.survival.pred$t13[subid] = min(srf.pred[srf.pred$time13>0, 1])
    srf.survival.pred$t14[subid] = min(srf.pred[srf.pred$time14>0, 1])
  }
  srf.survival.pred[srf.survival.pred==Inf] = 1
  end_time <- Sys.time()
  rf.cc.time[m] = end_time - start_time
  
  C.test = cIndex(dat.pa$time, dat.pa$status, CHF.test)
  C.cc.pa[m] = 1-C.test[1]
  
  rf.cc.pa[m] = sum((as.matrix(srf.survival.pred[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
  
  })
}

print("weighted")

summary(c(rf.cc.pa))
summary(c(C.cc.pa))
summary(c(rf.cc.time))

result = data.frame(nsim = c(1:nsim), rf.cc.pa=c(rf.cc.pa),C.cc.pa=c(C.cc.pa))

write.csv(result,file=paste0(path.results, 'cc_w.csv'), row.names = F)
