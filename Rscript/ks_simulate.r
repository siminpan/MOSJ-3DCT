library(e1071)
library(ggplot2)

library(foreach)
library(doParallel)

list_DK3 = sprintf("mean_dk%02d", 1:length(list_DRB))
list_DK4 = sprintf("sd_dk%02d", 1:length(list_DRB))

for (i in 1:length(list_DK)){
  assign(list_DK3[i], 
         mean(get(list_DK1[i])))
  
}

for (i in 1:length(list_DK)){
  assign(list_DK4[i], 
         sd(get(list_DK1[i])))
         
}

# skewness(get(nam)$data, type = 2)

# unlist(mget(eval(list_DK2)))

DK_0 <- unlist(mget(eval(list_DK1)))
DK_001 <- as.data.frame(x=DK_0)

h0_DK <- ggplot(DK_001, aes(x=DK_0)) + 
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", size=1.0) +
  xlim(-0.5, 0.5) +
  geom_histogram(data = DK_001, fill = colors()[11], alpha = 0.3)

h0_DK


# parallel-computing example ----
# https://www.r-bloggers.com/r-with-parallel-computing-from-user-perspectives/
cores <- detectCores(logical = FALSE)
cores = cores-7
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
chunk.size <- nrow(bar1)/cores

res2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
  {
    res <- matrix(0, nrow=chunk.size, ncol=2)
    for(x in ((i-1)*chunk.size+1):(i*chunk.size)){
      res[x - (i-1)*chunk.size,] <- c(x, i)
    }
    res
  }

# ks simulate ----
mean(nt03)
sd(nt03)
c1 = rnorm(1*10^6, mean(nt01), sd(nt01))
c2 = rnorm(1*10^6, mean(nt02), sd(nt02))
c3 = rnorm(1*10^6, mean(nt03), sd(nt03))

c1 = rnorm(1*10^6, 0, sd(nt01))
c2 = rnorm(1*10^6, 0, sd(nt02))
c3 = rnorm(1*10^6, 0, sd(nt03))


c0 = c(c1, c2, c3)

KS_1 <- ks.test(c1, c0)$stat
KS_2 <- ks.test(c2, c0)$stat
KS_3 <- ks.test(c3, c0)$stat
KS_c <- max(KS_1, KS_2, KS_3)


B <- 100
KS_0 <- data.frame(ks = rep(NA, B),
                   i = rep(NA, B)
)

cores <- detectCores(logical = FALSE)
cores = cores-7
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
chunk.size <- nrow(KS_0)/cores

res2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
  {
    res <- matrix(0, nrow=chunk.size, ncol=2)
    for(x in ((i-1)*chunk.size+1):(i*chunk.size)){
      y_0_1 <- sample(c0, replace = TRUE)
      y_0_2 <- sample(c0, replace = TRUE)
      y_0_3 <- sample(c0, replace = TRUE)
      res[x - (i-1)*chunk.size,] <- c(max(ks.test(y_0_1, c0)$st, ks.test(y_0_2, c0)$st, ks.test(y_0_3, c0)$st), i)
    }
    res
  }

## Compute p-values.
p_val_c <- mean(KS_c <= res2.p[,1])
