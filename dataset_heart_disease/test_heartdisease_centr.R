# train and test: DistLDA
library(flare)
toolbox_path = '~/Dropbox/Lu Tian/nips2016/Code/'
source(paste(toolbox_path, "distLDA_asym.R",sep=''))
source(paste(toolbox_path, "SLDA.R",sep=''))
source(paste(toolbox_path, "my_sugm.R",sep=''))
source(paste(toolbox_path, "my_dantzig.R",sep=''))
source(paste(toolbox_path, "my_dantzig_inner.R",sep=''))
source(paste(toolbox_path, "centrLDA_withmu.R",sep=''))

load("train_test_idx.Rdata")
fid_output = file("centralized_swipe.txt", "w")
for (C3 in 2^seq(-5,5)){
Nexpr = length(train_idx) # 10 repetitive experiments
result_table = matrix(0, 3, Nexpr) # container for all resutls
for (expr_no in 1:Nexpr) {
  
  # the container for train data in this expr
  data_train = vector("list",2*m)
  
  data_test = dataset[ test_idx[[expr_no]], 1:d ] # total test data in current experiment
  test_label = label[ test_idx[[expr_no]] ] # the label of current test data
  
  #current train data index
  current_train_idx = train_idx[[expr_no]]
  #the hospital id of current train data
  current_train_hsptid = hspt_id[ current_train_idx ]
  #the label of current train data
  current_train_label = label[ current_train_idx ]
  
  #the lambdas for different machines
  #lambda1 = c(0,0,0,0)
  #lambda2 = c(0,0,0,0)
  #C1 = 1
  #C2 = 1
  #C3 = 2^-9
  
  N = 0
  
  for (i in 1:m){
    # load train data for current hospital
    data_train[[i]] = t(dataset[ current_train_idx[ current_train_hsptid == i 
                                                  & current_train_label == T], 1:d])
    data_train[[i+m]] = t(dataset[ current_train_idx[ current_train_hsptid == i 
                                                  & current_train_label == F], 1:d])
    # the number for current hospital
    n = sum( current_train_hsptid == i)
    
    #lambda1[i] = C1 * sqrt(log(d) / n)
    #lambda2[i] = C2 * sqrt(log(d) / n)
    N = N + n
  }
  #C4 = 0.5
  #C5 = 0.5
  
  # truncate threshold of this experiment
  #t = C4 * sqrt(log(d)/N) + C5 * m * log(d) / N
  
  # lambda for centralized expr
  lambda_c = C3 * sqrt(log(d) / N)
  
  #res = distLDA_asym(data_train, lambda1, lambda2, t)
  res2 = centrLDA_withmu(data_train, lambda_c)
  
  #mubar = (res[[3]] + res[[4]]) / 2
  mubar2 = (res2[[2]]+res2[[3]]) / 2 #average, mubar and mubar2 should be the same
  
  difference = apply(data_test, 1, '-', mubar2)
  #predict1 = (crossprod(res[[1]],difference) > 0)
  #predict2 = (crossprod(res[[2]],difference) > 0)
  predict3 = (crossprod(res2[[1]],difference) > 0)
  #result_table[1, expr_no] = sum(predict1 != test_label)/length(test_label)
  #result_table[2, expr_no] = sum(predict2 != test_label)/length(test_label)
  result_table[3, expr_no] = sum(predict3 != test_label)/length(test_label)
}
cat(paste(as.character(C3),as.character(mean(result_table[3,])),
          as.character(sd(result_table[3,])),sep='\t'), file = fid_output)
cat('\n',file=fid_output)
}
close(fid_output)