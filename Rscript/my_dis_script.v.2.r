####
#### Analysis of bone image data.
####
#### Analysis goal: Test for whether there are differences between distance distributions 
#### across comparison groups.
####

library(e1071)
library(ggplot2)

# Load data ----
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
Y <- mget(eval(c(list_DK1, list_DKD1, list_DRB1, list_NT1)))

##
## Compare "treatment" and "control" groups. We pool the distances for the control 
## samples, take bootstrap samples from them, then compare distribution-test statistics 
## to those observed.
##

# Pool ----
## Pool the control samples.
# y_0 <- c(Y[["nt01"]],Y[["nt02"]], Y[["nt03"]], Y[["nt04"]], Y[["nt05"]], Y[["nt06"]], Y[["nt07"]], Y[["nt08"]], Y[["nt09"]], Y[["nt10"]])
# y_1 <- c(nt01, nt02, nt03, nt04, nt05, nt06, nt07, nt08, nt09, nt10)

## compare with DK group
y_0 <- unlist(mget(eval(list_DK1)))
## Compute KS test statistics comparing treatment samples to the pooled control samples.

# get KS stats ----
list_NT2 = sprintf("ks_nt%02d", 1:length(list_NT))
list_DRB2 = sprintf("ks_drb%02d", 1:length(list_DRB))
list_DK2 = sprintf("ks_dk%02d", 1:length(list_DK))
list_DKD2 = sprintf("ks_dkd%02d", 1:length(list_DKD))

for (i in c(1:4, 6:length(list_NT))){
  assign(list_NT2[i],
         ks.test(get(list_NT1[i]), y_0)$stat)

}
ks.test(nt03, nt02)
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

KS_NT <- max(unlist(mget(eval(list_NT2[c(1:4,6:10)]))))
# KS_NT <- min(unlist(mget(eval(list_NT2))))

KS_DkkMo <- max(unlist(mget(eval(list_DK2))))
# KS_DkkMo <- min(unlist(mget(eval(list_DK2))))

KS_DkkMoDRB <- max(unlist(mget(eval(list_DKD2))))
# KS_DkkMoDRB <- min(unlist(mget(eval(list_DKD2))))

KS_DRB <- max(unlist(mget(eval(list_DRB2))))
# KS_DkkMoDRB <- min(unlist(mget(eval(list_DKD2))))

# Bootstrap for null KS D ----
## Compute the null sampling distribution of the test statistic using bootstrap. Take 
## bootstrap samples from pooled control distances and compute the same statistic each 
## time.
l1 = ceiling(mean(lengths(Y)))
B <- 100
KS_0 <- rep(NA, B)
for(b in 1:B) {
  cat(b)

  y_0_1 <- sample(y_0, l1, replace = TRUE)
  y_0_2 <- sample(y_0, l1, replace = TRUE)
  y_0_3 <- sample(y_0, l1, replace = TRUE)
  KS_0[b] <- max(ks.test(y_0_1, y_0)$st, 
                 ks.test(y_0_2, y_0)$st, 
                 ks.test(y_0_3, y_0)$st)
}

B <- 100
KS_0 <- rep(NA, B)
list_loop = sprintf("y_0_%02d", 1:length(list_DK))
for(b in 1:B) {
  cat(b)
  for (i in 1:length(list_DK)){
    assign(list_loop[i],
           sample(get(list_DK1[i]), replace = TRUE)
           )
  }
    KS_0[b] <- max(ks.test(y_0_01, y_0)$st, 
                   ks.test(y_0_02, y_0)$st, 
                   ks.test(y_0_03, y_0)$st, 
                   ks.test(y_0_04, y_0)$st, 
                   ks.test(y_0_05, y_0)$st, 
                   ks.test(y_0_06, y_0)$st, 
                   ks.test(y_0_07, y_0)$st, 
                   ks.test(y_0_08, y_0)$st, 
                   ks.test(y_0_09, y_0)$st, 
                   ks.test(y_0_10, y_0)$st
                   )
}

KS_00 <- as.data.frame(x=KS_0)
write.csv(KS_00, "KS_00_Bsize100.csv")
## Compute p-values.
KS_0 <- read.csv("KS_dk_smalsize100.csv")
KS_0 = KS_0[2]
p_val_DkkMo <- mean(KS_DkkMo <= KS_0)
p_val_DkkMoDRB <- mean(KS_DkkMoDRB <= KS_0)
p_val_DRB <- mean(KS_DRB <= KS_0)
p_val_NT <- mean(KS_NT <= KS_0)


# plot ----

# y_01 <- data.frame(data = -y_NT_1, id = rep("y_01", length(y_NT_1)))
# y_02 <- data.frame(data = -y_NT_2, id = rep("y_02", length(y_NT_2)))
# y_03 <- data.frame(data = -y_NT_3, id = rep("y_03", length(y_NT_3)))
# y_04 <- data.frame(data = -y_NT_4, id = rep("y_04", length(y_NT_4)))
# y_05 <- data.frame(data = -y_NT_5, id = rep("y_05", length(y_NT_5)))
# y_06 <- data.frame(data = -y_NT_6, id = rep("y_06", length(y_NT_6)))

# |__nt ----
y_00 <- data.frame(data = numeric(), id = character())
skewness_00 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_00 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)

for(i in 1:10){
  nam <- paste("y_", sprintf("%02d", as.numeric(i)), sep = "")
  assign(nam, 
         data.frame(data = -get(paste("nt", sprintf("%02d", as.numeric(i)), sep = "")), 
                    id = rep(paste("y_", sprintf("%02d", as.numeric(i)), sep = ""), 
                             length(get(paste("nt", sprintf("%02d", as.numeric(i)), sep = ""))))
                    )
         )
  y_00 <- rbind(y_00, get(nam))
  skewness_00[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_00[i,2] <- paste("y_", sprintf("%02d", as.numeric(i)), sep = "")
  median_00[i,1] <- median(get(nam)$data) 
  median_00[i,2] <- paste("y_", sprintf("%02d", as.numeric(i)), sep = "")
}


h_NT <- ggplot(y_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  geom_histogram(data = subset(y_00,id == 'y_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(y_00,id == 'y_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(y_00,id == 'y_04'), fill = colors()[31], alpha = 0.3) +
  # geom_histogram(data = subset(y_00,id == 'y_05'), fill = colors()[41], alpha = 0.3) +
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

# |__dk ----
dk_00 <- data.frame(data = numeric(), id = character())
skewness_01 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_01 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)

for(i in 1:10){
  nam <- paste("dk_", sprintf("%02d", as.numeric(i)), sep = "")
  assign(nam, 
         data.frame(data = -get(paste("dk", sprintf("%02d", as.numeric(i)), sep = "")), 
                    id = rep(paste("dk_", sprintf("%02d", as.numeric(i)), sep = ""), 
                             length(get(paste("dk", sprintf("%02d", as.numeric(i)), sep = ""))))
         )
  )
  dk_00 <- rbind(dk_00, get(nam))
  skewness_01[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_01[i,2] <- paste("dk_", sprintf("%02d", as.numeric(i)), sep = "")
  median_01[i,1] <- median(get(nam)$data) 
  median_01[i,2] <- paste("dk_", sprintf("%02d", as.numeric(i)), sep = "")
}

h_DK <- ggplot(dk_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  xlim(-0.5, 0.5) +
  geom_histogram(data = subset(dk_00,id == 'dk_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(dk_00,id == 'dk_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_05'), fill = colors()[41], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_07'), fill = colors()[61], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_08'), fill = colors()[71], alpha = 0.3) +
  geom_histogram(data = subset(dk_00,id == 'dk_09'), fill = colors()[81], alpha = 0.3)+
  geom_histogram(data = subset(dk_00,id == 'dk_10'), fill = colors()[91], alpha = 0.3)

h_DK 

# |__dkdrb ----
dkdrb_00 <- data.frame(data = numeric(), id = character())
skewness_02 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_02 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:9){
  nam <- paste("dkd_", sprintf("%02d", as.numeric(i)), sep = "")
  assign(nam, 
         data.frame(data = -get(paste("dkd", sprintf("%02d", as.numeric(i)), sep = "")), 
                    id = rep(paste("dkd_", sprintf("%02d", as.numeric(i)), sep = ""), 
                             length(get(paste("dkd", sprintf("%02d", as.numeric(i)), sep = ""))))
         )
  )
  dkdrb_00 <- rbind(dkdrb_00, get(nam))
  skewness_02[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_02[i,2] <- paste("dkd_", sprintf("%02d", as.numeric(i)), sep = "")
  median_02[i,1] <- median(get(nam)$data) 
  median_02[i,2] <- paste("dkd_", sprintf("%02d", as.numeric(i)), sep = "")
}

h_DKDRB <- ggplot(dkdrb_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  xlim(-0.5, 0.5) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_05'), fill = colors()[41], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_07'), fill = colors()[61], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_08'), fill = colors()[71], alpha = 0.3) +
  geom_histogram(data = subset(dkdrb_00,id == 'dkd_09'), fill = colors()[81], alpha = 0.3)

h_DKDRB

# |__drb ----
drb_00 <- data.frame(data = numeric(), id = character())
skewness_03 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_03 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:10){
  nam <- paste("drb_", sprintf("%02d", as.numeric(i)), sep = "")
  assign(nam, 
         data.frame(data = -get(paste("drb", sprintf("%02d", as.numeric(i)), sep = "")), 
                    id = rep(paste("drb_", sprintf("%02d", as.numeric(i)), sep = ""), 
                             length(get(paste("drb", sprintf("%02d", as.numeric(i)), sep = ""))))))
  drb_00 <- rbind(drb_00, get(nam))
  skewness_03[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_03[i,2] <- paste("drb_", sprintf("%02d", as.numeric(i)), sep = "")
  median_03[i,1] <- median(get(nam)$data) 
  median_03[i,2] <- paste("drb_", sprintf("%02d", as.numeric(i)), sep = "")
}

h_DRB <- ggplot(drb_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  xlim(-0.5, 0.5) +
  geom_histogram(data = subset(drb_00,id == 'drb_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(drb_00,id == 'drb_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_05'), fill = colors()[41], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_07'), fill = colors()[61], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_08'), fill = colors()[71], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_09'), fill = colors()[81], alpha = 0.3) +
  geom_histogram(data = subset(drb_00,id == 'drb_10'), fill = colors()[91], alpha = 0.3)

h_DRB

# |__h ----
hh_00 <- data.frame(data = numeric(), id = character())
skewness_04 <- data.frame(skewness = numeric(), id = character(), stringsAsFactors = FALSE)
median_04 <- data.frame(median = numeric(), id = character(), stringsAsFactors = FALSE)
for(i in 1:7){
  nam <- paste("hh_", sprintf("%02d", as.numeric(i)), sep = "")
  assign(nam, 
         data.frame(data = -get(paste("hh", sprintf("%02d", as.numeric(i)), sep = "")), 
                    id = rep(paste("hh_", sprintf("%02d", as.numeric(i)), sep = ""), 
                             length(get(paste("hh", sprintf("%02d", as.numeric(i)), sep = ""))))))
  hh_00 <- rbind(hh_00, get(nam))
  skewness_04[i,1] <- skewness(get(nam)$data, type = 2)
  skewness_04[i,2] <- paste("hh_", sprintf("%02d", as.numeric(i)), sep = "")
  median_04[i,1] <- median(get(nam)$data) 
  median_04[i,2] <- paste("hh_", sprintf("%02d", as.numeric(i)), sep = "")
}

h_hh <- ggplot(hh_00, aes(x=data)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  xlim(-0.5, 0.5) +
  geom_histogram(data = subset(hh_00,id == 'hh_01'), fill = colors()[1], alpha = 0.3) + 
  geom_histogram(data = subset(hh_00,id == 'hh_02'), fill = colors()[11], alpha = 0.3) +
  geom_histogram(data = subset(hh_00,id == 'hh_03'), fill = colors()[21], alpha = 0.3) +
  geom_histogram(data = subset(hh_00,id == 'hh_04'), fill = colors()[31], alpha = 0.3) +
  geom_histogram(data = subset(hh_00,id == 'hh_05'), fill = colors()[41], alpha = 0.3) +
  geom_histogram(data = subset(hh_00,id == 'hh_06'), fill = colors()[51], alpha = 0.3) +
  geom_histogram(data = subset(hh_00,id == 'hh_07'), fill = colors()[61], alpha = 0.3)
  # geom_histogram(data = subset(hh_00,id == 'hh_08'), fill = colors()[71], alpha = 0.3) +
  # geom_histogram(data = subset(hh_00,id == 'hh_09'), fill = colors()[81], alpha = 0.3) +
  # geom_histogram(data = subset(hh_00,id == 'hh_10'), fill = colors()[91], alpha = 0.3)

h_hh
