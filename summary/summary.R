path <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/07_weighted/results/"

#mtry = c(2, 4)
mtry = c(2, 4, 6, 8)
#nodesize = c(500, 600, 700, 800)
nodesize = c(188, 375, 563, 750)
npoints = length(mtry)*length(nodesize)

# cc tuned
cc = data.frame()

for (i in 1:length(mtry)){
  for (j in 1:length(nodesize)){
    try({
    mtry.work = mtry[i]
    nodesize.work = nodesize[j]
    cc2 = read.csv(paste0(path, 'cc_m_', mtry.work, '_', nodesize.work, '.csv'))
    cc = rbind(cc, cc2)
    })
  }
}

cc.sort = cc[order(cc$nsim, cc$OOB.cc.pa),]
cc.new = cc.sort[seq(1, nrow(cc), npoints), ]
# MIE
summary(cc.new$rf.cc.pa)
summary(cc.new$C.cc.pa)

# weighted manually tuned

cc = data.frame()

for (i in 1:length(mtry)){
  for (j in 1:length(nodesize)){
    mtry.work = mtry[i]
    nodesize.work = nodesize[j]
    cc2 = read.csv(paste0(path, 'cc_w', mtry.work, '_', nodesize.work, '.csv'))
    cc = rbind(cc, cc2)
  }
}

cc.sort = cc[order(cc$nsim, -cc$OOB.cc.pa),]
cc.new = cc.sort[seq(1, nrow(cc), 32), ]
summary(cc.new$rf.cc.pa)
summary(cc.new$C.cc.pa)

# weighted self-tuned

w = read.csv(paste0(path, 'cc_w.csv'))
summary(w$rf.cc.pa)

# naive self-tuned

w = read.csv(paste0(path, 'cc_uw.csv'))
summary(w$rf.cc.pa)
