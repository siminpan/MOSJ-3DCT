library(e1071)
library(ggplot2)

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
