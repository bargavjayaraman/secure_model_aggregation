# centralized LDA, used for n fixed and m growing, N growing
# todo
centrLDA2 = function(data,lambda,m_s){ # lambda is a vector of same length with m_s
  m = length(data) / 2
  data1 = data[[1]]
  data2 = data[[m + 1]]
  counter = 1
  result = vector("list", length(m_s))
  for(j in 1:length(m_s)) {
    while(counter < m_s[j]){
      counter = counter + 1
      data1 = cbind(data1, data[[counter]])
      data2 = cbind(data2, data[[m + counter]])
    }
    SLDA(data1, data2,lambda[j])
    result[[j]] = SLDA(data1, data2,lambda[j])[[1]]
  }
  return (result)
}