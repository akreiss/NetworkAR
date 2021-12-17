library(igraph)
library(igraphdata)
library(Matrix)

## Load Data
data("rfid")

## Select sub-network
set.seed(2021)
sub_vertices <- c(sample(which(V(rfid)$Status=="ADM"),2),sample(which(V(rfid)$Status=="MED"),4),sample(which(V(rfid)$Status=="NUR"),9),sample(which(V(rfid)$Status=="PAT"),10))
rfid_sub <- induced_subgraph(rfid,sub_vertices)

## Extract sub-networks
G <- vector(mode="list",length=5)
for(k in 1:5) {
  edgeids <- which(E(rfid_sub)$Time<=k*24*60*60 & E(rfid_sub)$Time>(k-1)*24*60*60)
  G[[k]] <- simplify(subgraph.edges(rfid_sub,E(rfid_sub)[edgeids],delete.vertices = FALSE))
}

## Plot the networks
lay <- layout.auto(rfid_sub)
for(k in 1:5) {
  dev.new()
  plot(G[[k]],vertex.size=1,vertex.label="",main=sprintf("Time %d",k),layout=lay)
}

## Bring to correct format
nts <- vector(mode="list",length=5)
ncov <- vector(mode="list",length=5)
el <- vector(mode="list",length=5)

for(k in 1:5) {
  nts[[k]] <- as_adjacency_matrix(G[[k]],sparse=TRUE)
  help <- summary(nts[[k]]%*%t(nts[[k]]))
  ncov[[k]] <- cbind(i=help$i,j=help$j,x=help$x)
  help <- summary(nts[[k]])
  el[[k]] <- cbind(help$i,help$j)
}
net_seq <- list(nts=nts,t=0:4,n=4,p=25,ncov=ncov,el=el)

## Fit Model
source("./functions.R")
set.seed(2025)

B <- 100

backfit <- NetworkAR_backfitting(net_seq,B=B,limit=30.0,drop_first=FALSE,tol=0.00000000001)
h <- backfit$h_seq[[B]]
c <- backfit$c_seq[[B]]

## Generate Data and Fit
p <- 25
n <- 25
N <- 100
B <- 80

est_h <- vector(mode="list",length=N)
est_c <- vector(mode="list",length=N)

for(i in 1:N) {
  cat("Do Simulation ",i," of ",N,".\n")
  sim <- networkAR(n,p,c,exp(h),nts[[1]],ncov_return=TRUE)
  backfit <- NetworkAR_backfitting(sim,B=B,limit=30.0,drop_first=FALSE,tol=0.00000000001)
  est_h[[i]] <- backfit$h_seq[[B]]
  est_c[[i]] <- backfit$c_seq[[B]]
}

## Plot one simulated network sequence
for(k in 1:5) {
  dev.new()
  plot(graph_from_adjacency_matrix(sim$nts[[k]]),vertex.size=1,vertex.label="",main=sprintf("Time %d",k),layout=lay)
}

## Safe results
save(est_h,est_c,file="first_example.RData")

## Visualize four edges
edges_to_visualize <- matrix(c(1,4,
                              2,6,
                              6,23,
                              17,25),ncol=2,byrow = TRUE)
dev.new()
par(mfrow=c(2,2))
for(r in 1:dim(edges_to_visualize)[1]) {
  i <- edges_to_visualize[r,1]
  j <- edges_to_visualize[r,2]
  
  estimates <- rep(0,N)
  for(k in 1:N) {
    estimates[k] <- est_c[[k]][i,j]-c[i,j]
  }
  hist(estimates,main=sprintf("Edge (%d,%d)",i,j),breaks=30)
}
mtext("Estimates for c", side = 3, line = -2, outer = TRUE)

dev.new()
par(mfrow=c(2,2))
for(r in 1:dim(edges_to_visualize)[1]) {
  i <- edges_to_visualize[r,1]
  j <- edges_to_visualize[r,2]
  
  estimates <- rep(0,N)
  for(k in 1:N) {
    estimates[k] <- est_h[[k]][i,j]-h[i,j]
  }
  hist(estimates,main=sprintf("Edge (%d,%d)",i,j),breaks=30)
}
mtext("Estimates for h", side = 3, line = -2, outer = TRUE)

## Visualize estimates
dev.new()
par(mfrow=c(5,5))

