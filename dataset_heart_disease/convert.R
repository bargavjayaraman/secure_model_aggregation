#convert multi precision matrix to 4 digits
X = read.table("distributed_swipe.txt",header = F)
#X = as.numeric(X)
fid_output = file("d2.txt",open='w')
for(i in 1:nrow(X)){
  cat(sprintf("%f\t%.4f\t%.4f\n", X[i,1],X[i,2],X[i,3]), file = fid_output )
}
close(fid_output)