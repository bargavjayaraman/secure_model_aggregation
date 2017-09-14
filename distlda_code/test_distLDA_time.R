# test the computing time of distLDA (in one machine)
library(flare)
load("~/syndata_200_1.Rdata")

source("distLDA.R")
source("centrLDA.R")
source("SLDA.R")
source("my_sugm.R")
source("my_dantzig.R")
source("my_dantzig_inner.R")
n = N/m

C1 = 3
C2 = 1
C3 = 3

beta = Theta %*% (mu1 - mu2)
expr_num = length(experiments)
#t = C4 * 9 * sqrt(log(d)/N) * 6 + C5 * 11 * 9^2 * m * log(d) * 6 / N
t = 0.1

ptime <- system.time({
tmps = lapply( experiments, distLDA, C1*sqrt(log(d)/n), C2*sqrt(log(d)/n), t, doTruncateForNaive = TRUE)
})[3]
