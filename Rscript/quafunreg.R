library("glmnet")
library("MASS")
library("R.matlab")

library("pracma")
library("wavethresh")

# install.packages('gpuR')

source("/mnt/md0/zlyrebecca/sp/MOSJ-CT/script/Rscript/QFM-code/getquantlets.R")
source("/mnt/md0/zlyrebecca/sp/MOSJ-CT/script/Rscript/QFM-code/qfreg.R")
source("/mnt/md0/zlyrebecca/sp/MOSJ-CT/script/Rscript/QFM-code/inference.R")
source("/mnt/md0/zlyrebecca/sp/MOSJ-CT/script/Rscript/QFM-code/PcrQuant.R")

# Load data ----

# load("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg.RData")

# save(raw.dataset, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.RData")

load("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.RData")

workd = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints"
setwd(workd)
path1 = "/NT-t/"
list_NT = list.files(path = paste0(workd, path1), pattern = "_VTK_IO.csv")
list_NT = list_NT[-1]
path2 = "/DRB-t/"
list_DRB = list.files(path = paste0(workd, path2), pattern = "_VTK_IO.csv")
path3 = "/DkkMo-t/"
list_DK = list.files(path = paste0(workd, path3), pattern = "_VTK_IO.csv")
path4 = "/DkkMoDRB-t/"
list_DKD = list.files(path = paste0(workd, path4), pattern = "_VTK_IO.csv")
path5 = "/H-t/"
list_H = list.files(path = paste0(workd, path5), pattern = "_VTK_IO.csv")


list_NT1 = sprintf("nt%02d", 1:length(list_NT))
list_DRB1 = sprintf("drb%02d", 1:length(list_DRB))
list_DK1 = sprintf("dk%02d", 1:length(list_DK))
list_DKD1 = sprintf("dkd%02d", 1:length(list_DKD))
list_H1 = sprintf("hh%02d", 1:length(list_H))

for (i in 1:length(list_NT)){
  assign(list_NT1[i], 
         drop(as.matrix(read.csv(paste0(workd, path1, list_NT[i]), skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F))))
}
for (i in 1:length(list_DRB)){
  assign(list_DRB1[i], 
         drop(as.matrix(read.csv(paste0(workd, path2, list_DRB[i]), skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F))))
}
for (i in 1:length(list_DK)){
  assign(list_DK1[i], 
         drop(as.matrix(read.csv(paste0(workd, path3, list_DK[i]), skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F))))
}
for (i in 1:length(list_DKD)){
  assign(list_DKD1[i], 
         drop(as.matrix(read.csv(paste0(workd, path4, list_DKD[i]), skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F))))
}
for (i in 1:length(list_H)){
  assign(list_H1[i], 
         drop(as.matrix(read.csv(paste0(workd, path5, list_H[i]), skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F))))
}

## Organize into list for convenience.
x <- factor(c(rep("DkkMo", length(list_DK)), 
              rep("DkkMoDRB", length(list_DKD)), 
              rep("DRB", length(list_DRB)), 
              rep("NT", length(list_NT)))
)
raw.dataset <- mget(eval(c(list_DK1, list_DKD1, list_DRB1, list_NT1)))

raw.length <- rep(0, length(raw.dataset))
for (i in 1:length(raw.dataset)) {
  raw.length[[i]] <- length(raw.dataset[[i]])
}

# save.image(file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg.RData")

# Thinning the data ----

load("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.RData")

seq(1, length(raw.dataset[[i]]), 10)
raw.dataset.thin = raw.dataset

for (i in 1:length(raw.dataset)) {
  raw.dataset.thin[[i]] = raw.dataset[[i]][order(unlist(raw.dataset[[i]]))[seq(1, length(raw.dataset[[i]]), 10)]]
}

save(raw.dataset.thin, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.thin.RData")

# load thin data ----
load("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.thin.RData")
raw.dataset = raw.dataset.thin
rm(raw.dataset.thin)
########################  lasso selection! 30min########################

Sys.time() -> start
lasso.list <- vector("list", length(raw.dataset))
for (i in 24:length(raw.dataset)) {
  cat("\n",i, "start")
  y <- raw.dataset[[i]]
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
  set.seed(12345 + i)
  lasso_fit <- glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE)
  cvfit.lasso <- cv.glmnet(BETA_BASE_TOTAL_2, y, intercept = TRUE, nfolds = 3)
  zeros <- as.vector(coef(lasso_fit, s = cvfit.lasso$lambda.1se) == 0)
  selects <- seq(0, dim(BETA_BASE_TOTAL_2)[2], 1)[  zeros == FALSE ]
  lasso.list[[i]] <- selects
  save(lasso.list, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_lasso.thin.list.RData")
  rm(y,y.long,grid.p,CDFBETA,NQ,BNQ,BETA_BASE_TOTAL_2,lasso_fit,cvfit.lasso,zeros,selects)
  gc()
  cat("\n",i, "finished")
}

Times.over.lasso <- (Sys.time() - start)

# save(lasso.list, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_lasso.thin.list.RData")

######################## Compute common basis  ########################

# library("Matrix") # crossprod
library("bigmemory") # as.big.matrix
library("bigalgebra")
library("biganalytics") # apply,big.matrix-method

lasso.list1 <- lapply(lasso.list, CatNorm) ## 1 term fix
lasso.nonzero.obs1 <- lapply(lasso.list1, rep.list) ## D_i, generate 1 vector
lasso.counts.fit <- CountBasis(lasso.list1, lasso.nonzero.obs1)


n <- length(raw.dataset)

lasso_IncidenceVec_i_ <- vector("list", n) ### generate leave-one-out
for (i in 1:n) {
  lasso_fitIncd <- IncidenceVec(lasso.list1[-i], lasso.nonzero.obs1[-i])
  lasso_IncidenceVec_i_[[i]] <- lasso_fitIncd
}



p1024 <- signif(seq(0.001, 0.999, length = 1024), 4)
quantiles_p <- function(x, probs = p1024) {
  quantile(x, probs, type=6)
}


Qy = matrix( round(unlist( lapply( raw.dataset,  quantiles_p )  ),3) , 1024 )
write.csv(Qy, "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_Qy_1024.csv", row.names = FALSE )
Qy_3 <- read.csv("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_Qy_1024.csv", header = TRUE)
Qy <- Qy_3






leaveout.list <- lasso_IncidenceVec_i_
remain.counts <- lasso.counts.fit[[3]]
remain.basis <- lasso.counts.fit[[2]]
Y.list <- raw.dataset

active.set <- (remain.basis < 1000)
feasible.long <- sum(active.set)
n <- length(leaveout.list)

max.long <- max(unlist(lapply(Y.list, length)))
checks <- matrix(NA, nrow = n, ncol = feasible.long)
Values <- array(NA, c(max.long, n, feasible.long))


for (i in 1:n) {
  y <- raw.dataset[[i]]
  y.long <- length(y)
  grid.p <- seq(1 / (y.long + 1), y.long / (y.long + 1), 1 / (y.long + 1))
  CDFBETA <- GENERATE_BETA_CDF(a1, a2, grid.p)
  NQ <- qnorm(grid.p, 0, 1)
  BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
  BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))
  Psi <- cbind(rep(1, length(grid.p)), BETA_BASE_TOTAL_2)
  gc()
  for (j in 1:feasible.long) {
    colum_i_ <- leaveout.list[[i]][[1]]
    obs_i_ <- leaveout.list[[i]][[2]]
    
    SET1_i_ <- sort(colum_i_[  obs_i_ > remain.counts[active.set ][j]   ])
    smPsi_i_ <- Psi[, (SET1_i_) + 1 ]
    # https://pj.freefaculty.org/blog/?p=122
    gitsmp = ginv(crossprod(smPsi_i_), tol = sqrt(.Machine$double.eps))
    try1 = try(crossprod(smPsi_i_,y))
    try2 = try(gitsmp %*% try1)
    Values[(1:y.long), i, j] <- try(smPsi_i_ %*% try2)
    # Values[(1:y.long), i, j] <- try(smPsi_i_ %*% ginv(t(smPsi_i_) %*% smPsi_i_, tol = sqrt(.Machine$double.eps)) %*% t(smPsi_i_) %*% y)
    checks[i, j] <- length(SET1_i_)
    save(Values, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/loop_Values.thin.RData")
    save(checks, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/loop_checks.thin.RData")
    rm(colum_i_,obs_i_,SET1_i_,smPsi_i_,gitsmp,try1,try2)
    gc()
  }
  rm(y,y.long,grid.p,CDFBETA,NQ,BNQ,BETA_BASE_TOTAL_2,Psi)
  gc()
}

save.image(file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_2.RData")

##################### leave one-out cci ! for beta basis ###################

Y.mat <- matrix(NA, nrow = max.long, ncol = n)
for (i in 1:length(raw.dataset)) {
  Y.mat[(1:length(raw.dataset[[i]])), i] <- raw.dataset[[i]]
}


lasso.long <- sum(lasso.counts.fit[[2]] < 1000)
lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)

for (j in 1:lasso.long) {
  lasso.values[, j] <- agree.cci(Y.mat, Values[, , j])
}

lasso.Chary_i_ <- apply(lasso.values, 2, mean)
lasso.Chary1_i_ <- apply(lasso.values, 2, min)

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x1 <- lasso.counts.fit[[3]][ (lasso.counts.fit[[2]] < 1000) ]

##################### leave one-out cci ! for PCA ##################################################################################

## LCCC.FIT  --  computes the leave-one-out concordance correlation index based on PC regression
## TruncPCA  --  computes PCs basis with its maximum number (max.K)

fit.pca <- TruncPCA(Qy, length.train = NULL, max.K = length(raw.dataset))
lccc.pca <- LCCC.PCR.FIT(Y = Qy, B = fit.pca$v, remain.basis = lasso.counts.fit[[2]])



CDFBETA <- GENERATE_BETA_CDF(a1, a2, p1024)
NQ <- qnorm(p1024, 0, 1)
BNQ <- (NQ - mean(NQ)) / sqrt(sum((NQ - mean(NQ))^2))
BETA_BASE_TOTAL_2 <- cbind(BNQ, t(CDFBETA))


raw.dataset1 <- vector("list", length(raw.dataset))
for (i in 1:length(raw.dataset)) {
  raw.dataset1[[i]] <- Qy[, i]
}

##################### leave one-out cci ! for quanvelts ##################################################################################


lcccwave.lasso <- LCCC.WAVE.FIT(
  global.list = lasso.counts.fit[[1]], B = BETA_BASE_TOTAL_2,
  remain.counts = lasso.counts.fit[[3]],
  remain.basis = lasso.counts.fit[[2]], Y.list = raw.dataset1
)

Y.mat <- matrix(NA, nrow = 1024, ncol = n)
for (i in 1:length(raw.dataset)) {
  Y.mat[, i] <- raw.dataset1[[i]]
}


lcccwave.lasso.values <- matrix(NA, nrow = n, ncol = lasso.long)
for (j in 1:lasso.long) {
  lcccwave.lasso.values[, j] <- agree.cci(Y.mat, lcccwave.lasso[[1]][, , j])
}

wave.lasso.Chary_i_ <- apply(lcccwave.lasso.values, 2, mean)
wave.lasso.Chary1_i_ <- apply(lcccwave.lasso.values, 2, min)



texts1 <- expression(paste(bar(rho), " vs K"))
texts2 <- expression(paste(rho^0, " vs K"))


cex.lab <- 0.9 ### -> Size of labels for both x and y axis!
cex.axis <- 0.65 ### -> Size of coordinates for both x and y axis!
cex.main <- 0.7 ### -> Size of main topic ( USELESS since we will cut title )
cex.in <- 0.5 ### -> Size of all in boxplot (Ex: outlier )

yaxis.at <- c(-2, seq(0.85, 1.0, by = 0.05), 1.5)
yaxis.lab <- c("", "0.85", "0.90", "0.95", "1.00", "")
ylim <- c(0.84, 1.0)



itbasis <- sort(unique(c(lasso.x, lccc.pca[[2]])), decreasing = TRUE)

itxaxis <- seq(itbasis)
itidx1 <- rep(NA, length(itbasis))
itidx3 <- rep(NA, length(itbasis))
for (i in 1:length(itbasis)) {
  if (sum(lasso.x == itbasis[i]) != 0) {
    itidx1[i] <- i
  }
  if (sum(lccc.pca[[2]] == itbasis[i]) != 0) {
    itidx3[i] <- i
  }
}

xaxis.at <- c(-5, itxaxis[seq(3, length(itbasis), 3)], 35)
xaxis.lab <- c("", itbasis[ itxaxis[seq(3, length(itbasis), 3)] ], "")

xaxis.lab1 <- c("", lasso.x1[ itxaxis[seq(3, length(itbasis), 3)] ], "")


lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x1 <- lasso.counts.fit[[3]][ (lasso.counts.fit[[2]] < 1000) ]

########### start Reproduce Figure 3 ################################################################################


tiff("Figure3.tiff", width = 8, height = 4, units = "in", res = 120, bg = "transparent")

par(mfrow = c(1, 2), mar = c(4.5, 2.5, 3, 2))

plot(0, type = "o", xlab = "", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim, xlim = c(0, length(itbasis) + 1))
points(itidx1[is.na(itidx1) != TRUE], lasso.Chary1_i_, col = "red", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], lasso.Chary1_i_, col = "red")
points(itidx1[is.na(itidx1) != TRUE], lasso.Chary_i_, col = "blue", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], lasso.Chary_i_, col = "blue")
axis(side = 2, at = yaxis.at, label = yaxis.lab, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at, label = rep("", length(yaxis.at)), line = -0.5, , tick = TRUE, las = 1)

axis(side = 1, at = xaxis.at, label = xaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis)

axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at, label = xaxis.lab, line = 0.5, tick = FALSE, cex.axis = cex.axis)
mtext("C", side = 1, line = 0.5, at = 0.5, col = "blue", cex = 0.9)
mtext(expression(K[C]), side = 1, line = 1.7, at = 0.5, col = "blue", cex = 0.9)

title("A", cex.main = 1.2)


legend(4, 0.95, c(
  expression(paste(paste(paste(D^C, " ("), rho^0), ")")),
  expression(paste(paste(paste(D^C, " ("), bar(rho)), " )"))
),
lty = rep(1, 2), pch = rep(19, 2), col = c("red", "blue"),
cex = 1.1, bty = "n", ncol = 1
)


yaxis.at1 <- c(-2, seq(0.80, 1.0, by = 0.05), 1.5)
yaxis.lab1 <- c("", "0.80", "0.85", "0.90", "0.95", "1.00", "")
ylim1 <- c(0.80, 1.0)


plot(0, type = "o", xlab = "K", ylab = "", cex.lab = 1.1, cex.axis = 0.7, axes = FALSE, ylim = ylim1, xlim = c(0, length(itbasis) + 1))
points(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, mean)), col = "gray", pch = 19)
lines(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, mean)), col = "gray", lty = 2)
points(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, min)), col = "gray", pch = 19)
lines(seq(length(lccc.pca[[2]])), rev(apply(lccc.pca[[1]], 2, min)), col = "gray")

points(itidx1, rev(wave.lasso.Chary1_i_), col = "red", pch = 19)
lines(itidx1, rev(wave.lasso.Chary1_i_), col = "red")
points(itidx1[is.na(itidx1) != TRUE], rev(wave.lasso.Chary_i_), col = "blue", pch = 19)
lines(itidx1[is.na(itidx1) != TRUE], rev(wave.lasso.Chary_i_), col = "blue")

axis(side = 2, at = yaxis.at1, label = yaxis.lab1, line = -0.5, tick = FALSE, cex.axis = cex.axis, las = 1)
axis(side = 2, at = yaxis.at1, label = rep("", length(yaxis.at1)), line = -0.5, , tick = TRUE, las = 1)

xaxis.at22 <- c(-5, itxaxis[ seq(1, length(itbasis), 3)], 35)
xaxis.lab22 <- c("", rev(itbasis)[  seq(1, length(itbasis), 3) ], "")

axis(side = 1, at = xaxis.at22, label = rep("", length(xaxis.at22)), line = -0.5, tick = TRUE)
axis(side = 1, at = xaxis.at22, label = (xaxis.lab22), line = -0.5, tick = FALSE, cex.axis = cex.axis)

legend(6, 0.94, c(
  expression(paste("Quantlets (", rho^0, ")")),
  expression(paste("Quantlets (", bar(rho), " )")),
  expression(paste("PCA (", rho^0, ")")),
  expression(paste("PCA (", bar(rho), " )"))
),
lty = c(rep(1, 3), 2), pch = c(rep(19, 4)), col = c("red", "blue", "gray", "gray"),
cex = 1.1, bty = "n", ncol = 1
)
title("B", cex.main = 1.2)


dev.off()

########### end Reproduce Figure 4 ################################################################################

############# We choose # of basis ###############################################################################

lasso.x <- lasso.counts.fit[[2]][ (lasso.counts.fit[[2]] < 1000) ]
lasso.x.idx <- seq(length(lasso.list1))[ (lasso.counts.fit[[2]] < 1000) ]
list.order <- lasso.counts.fit[[1]]
unlist(lapply(list.order, length))



be <- 3
REDUCED_BASE9 <- BETA_BASE_TOTAL_2[, list.order[[be + 8 ] ]] # our choice in paper
REDUCED_BASE21 <- BETA_BASE_TOTAL_2[, list.order[[be + 20 ] ]] # Normal case in paper



save.image(file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_3.RData")
