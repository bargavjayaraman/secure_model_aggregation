#preprocessing heart disease data
data1 = read.table('processed.cleveland.data.txt',sep=',',na.strings='?')
data2 = read.table('processed.hungarian.data.txt',sep=',',na.strings='?')
data3 = read.table('processed.switzerland.data.txt',sep=',',na.strings='?')
data4 = read.table('processed.va.data.txt',sep=',',na.strings='?')

dataset = rbind(data1,data2,data3,data4)
hspt_id = c(rep(1,nrow(data1)),rep(2,nrow(data2)),
            rep(3,nrow(data3)),rep(4,nrow(data4)))
m = max(hspt_id)
# a function to extend some variables to dummy variables
myDummy <- function(x, lvls){
  n = length(x)
  d = length(lvls)
  output = matrix(0,n,d)
  for (i in 1:n){
    if(is.na(x[i])){
      output[i,] = NA
    }else{
    output[i,which(lvls == x[i])] = 1
    }
  }
  return(output)
}


extend_dummy <- function(X, dummy_list, lvls_list){
  newX = NULL
  for (i in 1:ncol(X)){
    if(any(dummy_list == i)){
      cur = myDummy(X[,i],lvls_list[[which(dummy_list == i)]])
    }
    else{
      cur = X[,i]
    }
    newX = cbind(newX, cur)
  }
  return(newX)
}

lvls = list(c(1,2,3,4),c(0,1,2),c(1,2,3),c(3,6,7))
dataset = extend_dummy(dataset, c(3,7,11,13), lvls_list = lvls)
## replace NA with averages
d = ncol(dataset)-1

for (i in 1:d){
  dataset[,i] = replace(dataset[,i], is.na(dataset[,i]), mean(dataset[,i], na.rm = T))
}

label = (dataset[,d+1] > 0)