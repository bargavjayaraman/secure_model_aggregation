# generate experiments
Nexpr = 10
train_idx = vector("list", Nexpr)
test_idx = vector("list", Nexpr)
Ntotal = nrow(dataset)
for(i in 1:Nexpr){
  for(j in 1:4){
    available_idx = which(hspt_id == j & label == T)
    train_idx[[i]] = c(train_idx[[i]],
                       sample(available_idx, floor(length(available_idx) / 2)))
    test_idx[[i]] = c(test_idx[[i]],setdiff(available_idx, train_idx[[i]]))
    available_idx = which(hspt_id == j & label == F)
    train_idx[[i]] = c(train_idx[[i]],
                       sample(available_idx, floor(length(available_idx) / 2)))
    test_idx[[i]] = c(test_idx[[i]],setdiff(available_idx, train_idx[[i]]))
  }
}
save(train_idx,test_idx,file="train_test_idx.Rdata")