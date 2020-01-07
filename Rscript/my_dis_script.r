####
#### Analysis of bone image data.
####
#### Analysis goal: Test for whether there are differences between distance distributions 
#### across comparison groups.
####

library(e1071)
library(ggplot2)

##
## Load data on distances between images. The images of two legs were superimposed on 
## each other, and pixel-wise distances between one image and its nearest point were 
## computed. The idea is that, if there is a tumor on one leg, it will be evident by 
## non-zero distances, where bone has either grown or degraded. 
##
setwd("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints")


y_DkkMo_1 <- drop(as.matrix(read.csv("DkkMo-t/27_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMo_2 <- drop(as.matrix(read.csv("DkkMo-t/28_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMo_3 <- drop(as.matrix(read.csv("DkkMo-t/30_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMo_4 <- drop(as.matrix(read.csv("DkkMo-t/31_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMo_5 <- drop(as.matrix(read.csv("DkkMo-t/32_VTK_IO.csv", 
                                     skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_DkkMoDRB_1 <- drop(as.matrix(read.csv("DkkMoDRB-t/36_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMoDRB_2 <- drop(as.matrix(read.csv("DkkMoDRB-t/37_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMoDRB_3 <- drop(as.matrix(read.csv("DkkMoDRB-t/38_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DkkMoDRB_4 <- drop(as.matrix(read.csv("DkkMoDRB-t/39_VTK_IO.csv", 
                                        skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_DRB_1 <- drop(as.matrix(read.csv("DRB-t/13_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DRB_2 <- drop(as.matrix(read.csv("DRB-t/14_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DRB_3 <- drop(as.matrix(read.csv("DRB-t/19_VTK_IO.csv", 
  skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_DRB_4 <- drop(as.matrix(read.csv("DRB-t/21_VTK_IO.csv", 
                                   skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

y_NT_1 <- drop(as.matrix(read.csv("NT-t/04_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_2 <- drop(as.matrix(read.csv("NT-t/05_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_3 <- drop(as.matrix(read.csv("NT-t/06_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_4 <- drop(as.matrix(read.csv("NT-t/07_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_5 <- drop(as.matrix(read.csv("NT-t/08_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_6 <- drop(as.matrix(read.csv("NT-t/09_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_7 <- drop(as.matrix(read.csv("NT-t/10_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_8 <- drop(as.matrix(read.csv("NT-t/11_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_9 <- drop(as.matrix(read.csv("NT-t/12_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))
y_NT_10 <- drop(as.matrix(read.csv("NT-t/99_VTK_IO.csv", skip = 1, colClasses = c(rep("numeric", 1), rep("NULL", 13)), header = F)))

## Organize into list for convenience.
x <- factor(c(rep("DkkMo", 5), rep("DkkMoDRB", 4), rep("DRB", 4), rep("NT", 10)))
Y <- list(y_DkkMo_1, y_DkkMo_2, y_DkkMo_3, y_DkkMo_4, y_DkkMo_5, 
          y_DkkMoDRB_1, y_DkkMoDRB_2, y_DkkMoDRB_3, y_DkkMoDRB_4,
          y_DRB_1, y_DRB_2, y_DRB_3, y_DRB_4, 
          y_NT_1, y_NT_2, y_NT_3, y_NT_4, y_NT_5, y_NT_6, y_NT_7, y_NT_8, y_NT_9, y_NT_10)
  
##
## Compare "treatment" and "control" groups. We pool the distances for the control 
## samples, take bootstrap samples from them, then compare distribution-test statistics 
## to those observed.
##

## Pool the control samples.
y_0 <- c(Y[[14]], Y[[15]], Y[[16]], Y[[17]], Y[[18]], Y[[19]], Y[[20]], Y[[21]], Y[[22]], Y[[23]])

## Compute KS test statistics comparing treatment samples to the pooled control samples.
KS_DkkMo_1 <- ks.test(Y[[1]], y_0)$stat
KS_DkkMo_2 <- ks.test(Y[[2]], y_0)$stat
KS_DkkMo_3 <- ks.test(Y[[3]], y_0)$stat
KS_DkkMo_4 <- ks.test(Y[[4]], y_0)$stat
KS_DkkMo_5 <- ks.test(Y[[5]], y_0)$stat
KS_DkkMo <- max(KS_DkkMo_1, KS_DkkMo_2, KS_DkkMo_3, KS_DkkMo_4, KS_DkkMo_5)

KS_DkkMoDRB_1 <- ks.test(Y[[6]], y_0)$stat
KS_DkkMoDRB_2 <- ks.test(Y[[7]], y_0)$stat
KS_DkkMoDRB_3 <- ks.test(Y[[8]], y_0)$stat
KS_DkkMoDRB_4 <- ks.test(Y[[9]], y_0)$stat
KS_DkkMoDRB <- max(KS_DkkMoDRB_1, KS_DkkMoDRB_2, KS_DkkMoDRB_3, KS_DkkMoDRB_4)

KS_DRB_1 <- ks.test(Y[[10]], y_0)$stat
KS_DRB_2 <- ks.test(Y[[11]], y_0)$stat
KS_DRB_3 <- ks.test(Y[[12]], y_0)$stat
KS_DRB_4 <- ks.test(Y[[13]], y_0)$stat
KS_DRB <- max(KS_DRB_1, KS_DRB_2, KS_DRB_3, KS_DRB_4)

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

KS_00 <- as.data.frame(x=KS_0)
write.csv(KS_00, "KS_00.csv")
## Compute p-values.
p_val_DkkMo <- mean(KS_DkkMo <= KS_0)
p_val_DkkMoDRB <- mean(KS_DkkMoDRB <= KS_0)
p_val_DRB <- mean(KS_DRB <= KS_0)


# plot ----

# y_01 <- data.frame(data = -y_NT_1, id = rep("y_01", length(y_NT_1)))
# y_02 <- data.frame(data = -y_NT_2, id = rep("y_02", length(y_NT_2)))
# y_03 <- data.frame(data = -y_NT_3, id = rep("y_03", length(y_NT_3)))
# y_04 <- data.frame(data = -y_NT_4, id = rep("y_04", length(y_NT_4)))
# y_05 <- data.frame(data = -y_NT_5, id = rep("y_05", length(y_NT_5)))
# y_06 <- data.frame(data = -y_NT_6, id = rep("y_06", length(y_NT_6)))

y_00 <- data.frame(data = numeric(), id = character())
skewness_00 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_00 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:10){
  nam <- paste("y_0", i, sep = "")
  assign(nam, 
         data.frame(data = -get(paste("y_NT_", i, sep = "")), 
                    id = rep(paste("y_0", i, sep = ""), 
                             length(get(paste("y_NT_", i, sep = ""))))
                    )
         )
  y_00 <- rbind(y_00, get(nam))
  skewness_00[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_00[i,2] <- paste("y_0", i, sep = "")
  median_00[i,1] <- median(get(nam)$data) 
  median_00[i,2] <- paste("y_0", i, sep = "")
}


h_NT <- ggplot(y_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(y_00,id == 'y_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(y_00,id == 'y_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_05'), fill = colors()[41], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_07'), fill = colors()[61], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_08'), fill = colors()[71], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_09'), fill = colors()[81], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_010'), fill = colors()[91], alpha = 0.3)

h_NT

h_NT2 <- ggplot(y_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(y_00,id == 'y_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(y_00,id == 'y_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_03'), fill = colors()[21], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_04'), fill = colors()[31], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_05'), fill = colors()[41], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_07'), fill = colors()[31], alpha = 0.3) # +
  # geom_histogram(data = subset(y_00,id == 'y_08'), fill = colors()[71], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_09'), fill = colors()[81], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_010'), fill = colors()[91], alpha = 0.3)

h_NT2

dk_00 <- data.frame(data = numeric(), id = character())
skewness_01 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_01 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:5){
  nam <- paste("dk_0", i, sep = "")
  assign(nam, 
         data.frame(data = -get(paste("y_DkkMo_", i, sep = "")), 
                    id = rep(paste("dk_0", i, sep = ""), 
                             length(get(paste("y_DkkMo_", i, sep = ""))))
         )
  )
  dk_00 <- rbind(dk_00, get(nam))
  skewness_01[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_01[i,2] <- paste("dk_0", i, sep = "")
  median_01[i,1] <- median(get(nam)$data) 
  median_01[i,2] <- paste("dk_0", i, sep = "")
}
h_DK <- ggplot(dk_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(dk_00,id == 'dk_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(dk_00,id == 'dk_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_05'), fill = colors()[41], alpha = 0.3)

h_DK 

dkdrb_00 <- data.frame(data = numeric(), id = character())
skewness_02 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_02 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:4){
  nam <- paste("dkdrb_0", i, sep = "")
  assign(nam, 
         data.frame(data = -get(paste("y_DkkMoDRB_", i, sep = "")), 
                    id = rep(paste("dkdrb_0", i, sep = ""), 
                             length(get(paste("y_DkkMoDRB_", i, sep = ""))))
         )
  )
  dkdrb_00 <- rbind(dkdrb_00, get(nam))
  skewness_02[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_02[i,2] <- paste("dkdrb_0", i, sep = "")
  median_02[i,1] <- median(get(nam)$data) 
  median_02[i,2] <- paste("dkdrb_0", i, sep = "")
}
h_DKDRB <- ggplot(dkdrb_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkdrb_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(dkdrb_00,id == 'dkdrb_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkdrb_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkdrb_04'), fill = colors()[31], alpha = 0.3)

h_DKDRB

drb_00 <- data.frame(data = numeric(), id = character())
skewness_03 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_03 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:4){
  nam <- paste("drb_0", i, sep = "")
  assign(nam, 
         data.frame(data = -get(paste("y_DRB_", i, sep = "")), 
                    id = rep(paste("drb_0", i, sep = ""), 
                             length(get(paste("y_DRB_", i, sep = ""))))
         )
  )
  drb_00 <- rbind(drb_00, get(nam))
  skewness_03[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_03[i,2] <- paste("drb_0", i, sep = "")
  median_03[i,1] <- median(get(nam)$data) 
  median_03[i,2] <- paste("drb_0", i, sep = "")
}
h_DRB <- ggplot(drb_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(drb_00,id == 'drb_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(drb_00,id == 'drb_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_04'), fill = colors()[31], alpha = 0.3)

h_DRB


