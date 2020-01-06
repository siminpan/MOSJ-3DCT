####
#### Analysis of bone image data.
####
#### Analysis goal: Test for whether there are differences between distance distributions 
#### across comparison groups.
####

##
## Load data on distances between images. The images of two legs were superimposed on 
## each other, and pixel-wise distances between one image and its nearest point were 
## computed. The idea is that, if there is a tumor on one leg, it will be evident by 
## non-zero distances, where bone has either grown or degraded. 
##
setwd("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints")


y_DkkMo_1 <- read.csv("DkkMo-t/27_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)
y_DkkMo_2 <- drop(as.matrix(read.csv("02_Distances/DkkMo-t/dist_csvs/dist_2.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMo_3 <- drop(as.matrix(read.csv("02_Distances/DkkMo-t/dist_csvs/dist_3.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_DkkMoDRB_1 <- drop(as.matrix(read.csv("02_Distances/DkkMoDRB-t/dist_csvs/dist_1.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMoDRB_2 <- drop(as.matrix(read.csv("02_Distances/DkkMoDRB-t/dist_csvs/dist_2.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMoDRB_3 <- drop(as.matrix(read.csv("02_Distances/DkkMoDRB-t/dist_csvs/dist_3.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_DRB_1 <- drop(as.matrix(read.csv("02_Distances/DRB-t/dist_csvs/dist_1.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DRB_2 <- drop(as.matrix(read.csv("02_Distances/DRB-t/dist_csvs/dist_2.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DRB_3 <- drop(as.matrix(read.csv("02_Distances/DRB-t/dist_csvs/dist_3.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_NT_1 <- drop(as.matrix(read.csv("02_Distances/NT-t/dist_csvs/dist_1.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_2 <- drop(as.matrix(read.csv("02_Distances/NT-t/dist_csvs/dist_2.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_3 <- drop(as.matrix(read.csv("02_Distances/NT-t/dist_csvs/dist_3.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

## Organize into list for convenience.
x <- factor(c(rep("DkkMo", 3), rep("DkkMoDRB", 3), rep("DRB", 3), rep("NT", 3)))
Y <- list(y_DkkMo_1, y_DkkMo_2, y_DkkMo_3, y_DkkMoDRB_1, y_DkkMoDRB_2, y_DkkMoDRB_3, 
  y_DRB_1, y_DRB_2, y_DRB_3, y_NT_1, y_NT_2, y_NT_3)
  
##
## Compare "treatment" and "control" groups. We pool the distances for the control 
## samples, take bootstrap samples from them, then compare distribution-test statistics 
## to those observed.
##

## Pool the control samples.
y_0 <- c(Y[[10]], Y[[11]], Y[[12]])

## Compute KS test statistics comparing treatment samples to the pooled control samples.
KS_DkkMo_1 <- ks.test(Y[[1]], y_0)$stat
KS_DkkMo_2 <- ks.test(Y[[2]], y_0)$stat
KS_DkkMo_3 <- ks.test(Y[[3]], y_0)$stat
KS_DkkMo <- max(KS_DkkMo_1, KS_DkkMo_2, KS_DkkMo_3)

KS_DkkMoDRB_1 <- ks.test(Y[[4]], y_0)$stat
KS_DkkMoDRB_2 <- ks.test(Y[[5]], y_0)$stat
KS_DkkMoDRB_3 <- ks.test(Y[[6]], y_0)$stat
KS_DkkMoDRB <- max(KS_DkkMoDRB_1, KS_DkkMoDRB_2, KS_DkkMoDRB_3)

KS_DRB_1 <- ks.test(Y[[7]], y_0)$stat
KS_DRB_2 <- ks.test(Y[[8]], y_0)$stat
KS_DRB_3 <- ks.test(Y[[9]], y_0)$stat
KS_DRB <- max(KS_DRB_1, KS_DRB_2, KS_DRB_3)

## Compute the null sampling distribution of the test statistic using bootstrap. Take 
## bootstrap samples from pooled control distances and compute the same statistic each 
## time.
B <- 1000
KS_0 <- rep(NA, B)
for(b in 1:B) {
  cat(b)

  y_0_1 <- sample(y_0, replace = TRUE)
  y_0_2 <- sample(y_0, replace = TRUE)
  y_0_3 <- sample(y_0, replace = TRUE)
  KS_0[b] <- max(ks.test(y_0_1, y_0)$st, ks.test(y_0_2, y_0)$st, ks.test(y_0_3, y_0)$st)
}

## Compute p-values.
p_val_DkkMo <- mean(KS_DkkMo <= KS_0)
p_val_DkkMoDRB <- mean(KS_DkkMoDRB <= KS_0)
p_val_DRB <- mean(KS_DRB <= KS_0)



