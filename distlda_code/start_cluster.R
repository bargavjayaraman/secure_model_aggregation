### Cluster at server
library(doParallel)

cl <- makeCluster(24)

registerDoParallel(cl)
###
