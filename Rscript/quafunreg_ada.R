# module load R/4.0.3-foss-2020a-recommended-mt

## module load rstudio/1.3.1093-foss-2020a-Java-11-R-4.0.0
## rserver --server-daemonize=0 --www-port 8787 --rsession-which-r=$(which R)

# /scratch/user/simonpan/rlibrary
.libPaths( c( "/scratch/user/simonpan/rlibrary" , .libPaths() ) )

library("MASS")

library("glmnet")
library("R.matlab")
library("pracma")
library("wavethresh")

source("/scratch/user/simonpan/ct/QFM-code/getquantlets.R")
source("/scratch/user/simonpan/ct/QFM-code/qfreg.R")
source("/scratch/user/simonpan/ct/QFM-code/inference.R")
source("/scratch/user/simonpan/ct/QFM-code/PcrQuant.R")

# Load data ----

# load("/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg.RData")

# save(raw.dataset, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_raw.data.RData")

load("/scratch/user/simonpan/ct/quafunreg_raw.data.RData")

########################  lasso selection! 30min########################

Sys.time() -> start
lasso.list <- vector("list", length(raw.dataset))
for (i in 1:length(raw.dataset)) {
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
}

Times.over.lasso <- (Sys.time() - start)

save.image(file = "/scratch/user/simonpan/ct/quafunreg_lasso.RData")

save(lasso.list, file = "/mnt/md0/zlyrebecca/sp/MOSJ-CT/05.6_3Dpoints/quafunreg_lasso.list.RData")