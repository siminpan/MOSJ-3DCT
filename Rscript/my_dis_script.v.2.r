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

list_NT1 = sprintf("nt%02d", 1:length(list_NT))
list_DRB1 = sprintf("drb%02d", 1:length(list_DRB))
list_DK1 = sprintf("dk%02d", 1:length(list_DK))
list_DKD1 = sprintf("dkd%02d", 1:length(list_DKD))

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

## Organize into list for convenience.
x <- factor(c(rep("DkkMo", length(list_DK)), 
              rep("DkkMoDRB", length(list_DKD)), 
              rep("DRB", length(list_DRB)), 
              rep("NT", length(list_NT)))
            )
Y <- mget(eval(c(list_DK1, list_DKD1, list_DRB1, list_NT1)))

##
## Compare "treatment" and "control" groups. We pool the distances for the control 
## samples, take bootstrap samples from them, then compare distribution-test statistics 
## to those observed.
##

## Pool the control samples.
y_0 <- c(Y[[30]],Y[[31]], Y[[32]], Y[[33]], Y[[34]], Y[[35]], Y[[36]], Y[[37]], Y[[38]], Y[[39]])

## Compute KS test statistics comparing treatment samples to the pooled control samples.

list_NT2 = sprintf("ks_nt%02d", 1:length(list_NT))
list_DRB2 = sprintf("ks_drb%02d", 1:length(list_DRB))
list_DK2 = sprintf("ks_dk%02d", 1:length(list_DK))
list_DKD2 = sprintf("ks_dkd%02d", 1:length(list_DKD))

# for (i in 1:length(list_NT)){
#   assign(list_NT2[i], 
#          ks.test(get(list_NT1[i]), y_0)$stat)
#   
# }
for (i in 1:length(list_DRB)){
  assign(list_DRB2[i], 
         ks.test(get(list_DRB1[i]), y_0)$stat)
  
}
for (i in 1:length(list_DK)){
  assign(list_DK2[i], 
         ks.test(get(list_DK1[i]), y_0)$stat)
  
}
for (i in 1:length(list_DKD)){
  assign(list_DKD2[i], 
         ks.test(get(list_DKD1[i]), y_0)$stat)
  
}

KS_DkkMo <- max(unlist(mget(eval(list_DK2))))

KS_DkkMoDRB <- max(unlist(mget(eval(list_DKD2))))

KS_DRB <- max(unlist(mget(eval(list_DRB2))))

## Compute the null sampling distribution of the test statistic using bootstrap. Take 
## bootstrap samples from pooled control distances and compute the same statistic each 
## time.
B <- 10000
KS_0 <- rep(NA, B)
for(b in 1:B) {
  cat(b)

  y_0_1 <- sample(y_0, replace = TRUE)
  y_0_2 <- sample(y_0, replace = TRUE)
  y_0_3 <- sample(y_0, replace = TRUE)
  KS_0[b] <- max(ks.test(y_0_1, y_0)$st, 
                 ks.test(y_0_2, y_0)$st, 
                 ks.test(y_0_3, y_0)$st)
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
  xlim(-0.5, 0.5) +
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
  xlim(-0.5, 0.5) +
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
  xlim(-0.5, 0.5) +
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
  xlim(-0.5, 0.5) +
  geom_histogram(data = subset(drb_00,id == 'drb_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(drb_00,id == 'drb_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_04'), fill = colors()[31], alpha = 0.3)

h_DRB


