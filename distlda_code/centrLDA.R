#centralized LDA
centrLDA = function(data,lambda){
  m = length(data) / 2
  data1 = data[[1]]
  data2 = data[[m + 1]]
  for(i in 2:m){
    data1 = cbind(data1, data[[i]])
    data2 = cbind(data2, data[[m + i]])
  }
  tmp = SLDA(data1, data2,lambda)
  return (tmp[[1]])
}