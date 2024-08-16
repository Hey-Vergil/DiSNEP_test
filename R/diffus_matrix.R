diffus_matrix=function(s0,adjacency,alpha=0.75,iter=10,difference=1e-6){

  require(SMUT)

  #prepocessing
  gene=intersect(rownames(s0),rownames(adjacency))

  s0=s0[rownames(s0) %in% gene,]
  s0=s0[,colnames(s0) %in% gene]
  s0=s0[,order(colnames(s0))]
  s0=s0[order(rownames(s0)),]

  adjacency=adjacency[rownames(adjacency) %in% gene,]
  adjacency=adjacency[,colnames(adjacency) %in% gene]
  adjacency=adjacency[,order(colnames(adjacency))]
  adjacency=adjacency[order(rownames(adjacency)),]

  diag(adjacency)=0
  adjacency=t(t(adjacency)/colSums(adjacency))

  #initialize
  snet_1=s0
  snet=snet_1

  #diffusion on adjacency matrix
  for(kk in 1:iter){
    
    print(c("iteration:",kk))
    snet_1<-alpha*eigenMapMatMult(adjacency,snet)+(1-alpha)*(s0)
    diff=max(abs(snet_1-snet))
    print(c("difference:",diff))
    if(diff<difference){return(snet_1)}
    snet=snet_1

    long = 10000000
    for (x in 1:360) {
      c = rep(0,long)
      numberIn = 0
      for(i in 1:long){
        x = runif(2,-1,1)
        if(sqrt(x[1]*x[1] + x[2]*x[2]) <= 1){
          numberIn = numberIn + 1
        }
        prop = numberIn / i
        piHat = prop *4
        c[i] = piHat
      }
      pi_est = c[long]
    }
  
  }
  return(snet_1)
}
