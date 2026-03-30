path <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/07_weighted/github/"
path.results <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/07_weighted/github/"


library(rpart)
library(partykit)
library(survival)
library(tidyr)
library(intsurv)

#detach("package:randomForestSRC", unload = TRUE)
install.packages("C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/07_weighted/github/randomForestSRC", repos = NULL, type = "source")
library(randomForestSRC)

nsim = 100

ntree=500
nodesize = 750
mtry = 2
tau = 0.0162
#p=4
#true = c(0.1, 0.3, -0.3, -0.1)
p=10
p.full = 2
gamma.weibull = 0.7
beta = matrix(c(rep(0,p.full-2), c(0.1, 0.1, 0.3, -0.3, -0.2, -0.3, 0.01, -0.01, 0.25, -0.25), rep(0,p-p.full-8)), nrow=p)

rf.cc.pa = matrix(NA, nrow = 1, ncol = nsim)
OOB.cc.pa = matrix(NA, nrow = 1, ncol = nsim)
C.cc.pa = matrix(NA, nrow = 1, ncol = nsim)
C2.cc.pa = matrix(NA, nrow = 1, ncol = nsim)


#m=2
for (m in 1:nsim){
  try({
    cat(m)
    
    # read datasets
    cc = read.csv(paste0(path, 'train', m, '.csv'))
    cc = subset(cc, select = -sc)
    test = read.csv(paste0(path, 'test', m, '.csv'))
    #test = test[1:10,]
    
    # true survival
    dat.pa = test
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
    
    chf.survival.pred = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                                   t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                                   t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                                   t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                                   t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
    
    cc = subset(cc, select = -c(subid))
    cc = cc[is.na(cc[,p]) == 0,]
    
    ind.mat = matrix(NA, nrow = nrow(cc), ncol = ntree)
    inbag.mat = matrix(NA, nrow = nrow(cc), ncol = ntree)
    OOB.mat = matrix(1, nrow = nrow(cc), ncol = ntree)
    
    ind.mat = as.matrix(rmultinom(n=ntree, size = 1500, prob = rep(1/1500, 1500)))
    ind.mat = ind.mat[1:nrow(cc),]
    
    for (l in 1:nrow(cc)){
      for (k in 1:ntree){
        inbag.mat[l, k] = as.numeric(ind.mat[l,k]>0)
      }
    }
    
    OOB.mat = OOB.mat - inbag.mat
    
    time = data.frame(time = cc$time, weight = cc$weight)
    time = time[order(time$time),]
    
    time.pa = data.frame(time = sort(dat.pa$time), ind = 1)
    
    OOB.CHF.cum = matrix(0, nrow = nrow(cc), ncol = 1)
    CHF.test.cum = matrix(0, nrow = nrow(dat.pa), ncol = 1)
    
    for (l in 1:ntree){
      cat(l)
      #set.seed(l*10)
      ind = ind.mat[,l]
      boot.sample = as.data.frame(lapply(cc, rep, ind))
      boot.sample$weight = round(boot.sample$weight)
      boot.sample <- as.data.frame(lapply(boot.sample, rep, boot.sample$weight))
      boot.sample = subset(boot.sample, select = -c(weight))
      srf.fit = rfsrc(Surv(time, status) ~ ., data = boot.sample,splitrule = "custom", ntree = 1, mtry = mtry, sampsize = nrow(boot.sample), nodesize = nodesize)
      srf.survival.pred = data.frame(subid = rep(0, nrow(dat.pa)), t1 = rep(0, nrow(dat.pa)), t2 = rep(0, nrow(dat.pa)), 
                                     t3 = rep(0, nrow(dat.pa)), t4 = rep(0, nrow(dat.pa)), t5 = rep(0, nrow(dat.pa)), 
                                     t6 = rep(0, nrow(dat.pa)), t7 = rep(0, nrow(dat.pa)), t8 = rep(0, nrow(dat.pa)), 
                                     t9 = rep(0, nrow(dat.pa)), t10 = rep(0, nrow(dat.pa)), t11 = rep(0, nrow(dat.pa)),
                                     t12 = rep(0, nrow(dat.pa)), t13 = rep(0, nrow(dat.pa)), t14 = rep(0, nrow(dat.pa)))
      
      OOB.CHF = matrix(NA, nrow = nrow(cc), ncol = 1)
      for (i in 1:nrow(cc)){
        srf.pred = data.frame(chf = predict(srf.fit, newdata = cc[i,1:p])$chf[1,], time= predict(srf.fit, newdata = cc[i,1:p])$time.interest)
        chf.temp = merge(srf.pred, time, by="time", all.y=T)
        chf.temp = chf.temp %>% fill(chf)
        OOB.CHF[i] = sum(chf.temp$chf*chf.temp$weight, na.rm = T)*OOB.mat[i,l]
      }
      OOB.CHF.cum = OOB.CHF.cum+OOB.CHF
      
      CHF.test = matrix(NA, nrow = nrow(dat.pa), ncol = 1)
      for (i in 1:nrow(dat.pa)){
        srf.pred = data.frame(chf = predict(srf.fit, newdata = dat.pa[i,1:p])$chf[1,], time= predict(srf.fit, newdata = dat.pa[i,1:p])$time.interest)
        chf.test.temp = merge(srf.pred, time.pa, by="time", all=T)
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
        srf.survival.pred$t1[subid] = max(srf.pred[srf.pred$time1>0, 1])
        srf.survival.pred$t2[subid] = max(srf.pred[srf.pred$time2>0, 1])
        srf.survival.pred$t3[subid] = max(srf.pred[srf.pred$time3>0, 1])
        srf.survival.pred$t4[subid] = max(srf.pred[srf.pred$time4>0, 1])
        srf.survival.pred$t5[subid] = max(srf.pred[srf.pred$time5>0, 1])
        srf.survival.pred$t6[subid] = max(srf.pred[srf.pred$time6>0, 1])
        srf.survival.pred$t7[subid] = max(srf.pred[srf.pred$time7>0, 1])
        srf.survival.pred$t8[subid] = max(srf.pred[srf.pred$time8>0, 1])
        srf.survival.pred$t9[subid] = max(srf.pred[srf.pred$time9>0, 1])
        srf.survival.pred$t10[subid] = max(srf.pred[srf.pred$time10>0, 1])
        srf.survival.pred$t11[subid] = max(srf.pred[srf.pred$time11>0, 1])
        srf.survival.pred$t12[subid] = max(srf.pred[srf.pred$time12>0, 1])
        srf.survival.pred$t13[subid] = max(srf.pred[srf.pred$time13>0, 1])
        srf.survival.pred$t14[subid] = max(srf.pred[srf.pred$time14>0, 1])
      }
      CHF.test.cum = CHF.test.cum+CHF.test/ntree
      srf.survival.pred[srf.survival.pred==Inf] = 1
      chf.survival.pred = chf.survival.pred+srf.survival.pred/ntree
    }
    
    srf.survival.pred = exp(-chf.survival.pred)
    
    OOB.CHF.cum = OOB.CHF.cum/rowSums(OOB.mat)
    
    
    rf.cc.pa[m] = sum((as.matrix(srf.survival.pred[,-1]) - as.matrix(true.survival[,-1]))^2)/nrow(dat.pa)/14
    C = cIndex(cc$time, cc$status, OOB.CHF.cum, cc$weight)
    OOB.cc.pa[m] = 1-C[1]
    
    C.test = cIndex(dat.pa$time, dat.pa$status, CHF.test.cum)
    C.cc.pa[m] = 1-C.test[1]
    
    C2.test = cIndex(dat.pa$time, dat.pa$status, 1-srf.survival.pred$t14)
    C2.cc.pa[m] = 1-C2.test[1]
    
  })
}

print(paste0("cc, mtry = ", mtry, ", nodesize = ", nodesize, ", ntree = ", ntree))

summary(c(rf.cc.pa))
summary(c(OOB.cc.pa))
summary(c(C.cc.pa))
summary(c(C2.cc.pa))

result = data.frame(nsim = c(1:nsim), rf.cc.pa=c(rf.cc.pa), OOB.cc.pa=c(OOB.cc.pa), C.cc.pa=c(C.cc.pa), C2.cc.pa=c(C2.cc.pa))

write.csv(result,file=paste0(path.results, 'cc_m_', mtry, '_', nodesize, '.csv'), row.names = F)
