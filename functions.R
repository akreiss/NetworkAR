library(Matrix)
library(nloptr)

dyn.load("optimize.dll")


## Simulates a time series of networks generated by the Network AR model with
## clustering and assortative mixing.
networkAR<- function(n,p,c,eh,X0,ncov_return=FALSE) {
  network_sequence <- vector(mode="list",length=n+1)
  times <- 0:n
  
  ## Allocate output Data if required
  if(ncov_return==TRUE) {
    ncov_out <- vector(mode="list",length=n+1)
    el_out <- vector(mode="list",length=n+1)
  }

  ## Set initial network
  network_sequence[[1]] <- 1*X0
  
  ## Simulate the network time series
  for(k in 1:n) {
    ## Compute number of common outgoing vertices
    ncov <- network_sequence[[k]]%*%t(network_sequence[[k]])
    if(ncov_return==TRUE) {
      help <- summary(ncov)
      ncov_out[[k]] <- cbind(i=help$i,j=help$j,x=help$x)
      help <- summary(network_sequence[[k]])
      el_out[[k]] <- cbind(help$i,help$j)
    }
    
    ## Compute alpha and beta
    ec <- exp(c*ncov)
    dn <- (1+ec+eh)
    alpha <- ec/dn
    beta  <- eh/dn
    
    ## Generate randomness
    eps0 <- matrix(runif(p*p),ncol=p,nrow=p)
    one_ind <- which(eps0<=alpha,arr.ind=TRUE)
    eps  <- sparseMatrix(i=one_ind[,1],j=one_ind[,2],x=1,dims=c(p,p))
    eps[eps0<=alpha+beta & eps0>alpha] <- -1
    
    ## Rewrite the matrix in list form
    eps_list_h <- summary(eps)

    ## Compute new adjacency matrix
    network_sequence[[k+1]] <- network_sequence[[k]]
    network_sequence[[k+1]][cbind(eps_list_h$i[eps_list_h$x== 1],eps_list_h$j[eps_list_h$x== 1])] <- 1
    network_sequence[[k+1]][cbind(eps_list_h$i[eps_list_h$x==-1],eps_list_h$j[eps_list_h$x==-1])] <- 0
    network_sequence[[k+1]] <- drop0(network_sequence[[k+1]])
  }
  ## Compute number of common outgoing vertices for the last network
  if(ncov_return==TRUE) {
    help <- summary(network_sequence[[n+1]]%*%t(network_sequence[[n+1]]))
    ncov_out[[n+1]] <- cbind(i=help$i,j=help$j,x=help$x)
    help <- summary(network_sequence[[n+1]])
    el_out[[n+1]] <- cbind(help$i,help$j)
  }

  if(ncov_return==TRUE) {
    return(list(nts=network_sequence,t=times,n=n,p=p,ncov=ncov_out,el=el_out))
  } else {
    return(list(nts=network_sequence,t=times,n=n,p=p)) 
  }
}

## Computes the negative log-likelihood
log_likelihood <- function(c,h,net_seq,vectorize=TRUE,combine=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  
  if(combine) {
    combineR <- 1
  } else {
    combineR <- 0
  }
  
  out <- .Call("log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p),as.integer(combineR))
  
  if(combine) {
    return(out)
  } else {
    return(matrix(out,ncol=p,nrow=p))
  }
  
  return()
}

log_likelihoodh <- function(h,c,net_seq,vectorize=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  return(.Call("log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p),as.integer(1)))
}

ddc_log_likelihood <- function(c,h,net_seq,vectorize=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  
  o <- .Call("ddc_log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p))
  if(vectorize) {
    o <- matrix(o,ncol=net_seq$p,nrow=net_seq$p)
  }
  
  return(o)
}

ddh_log_likelihood <- function(c,h,net_seq,vectorize=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  
  o <- .Call("ddh_log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p))
  if(vectorize) {
    o <- matrix(o,ncol=net_seq$p,nrow=net_seq$p)
  }
  
  return(o)
}

d2dc2_log_likelihood <- function(c,h,net_seq,vectorize=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  
  o <- .Call("d2dc2_log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p))
  if(vectorize) {
    o <- matrix(o,ncol=net_seq$p,nrow=net_seq$p)
  }
  
  return(o)
}

d2dh2_log_likelihood <- function(c,h,net_seq,vectorize=TRUE) {
  if(vectorize) {
    cmat <- matrix(c,ncol=net_seq$p,nrow=net_seq$p)
    hmat <- matrix(h,ncol=net_seq$p,nrow=net_seq$p)
  } else {
    cmat <- c
    hmat <- h
  }
  
  o <- .Call("d2dh2_log_likelihood",net_seq$el,net_seq$ncov,cmat,hmat,as.integer(net_seq$n+1),as.integer(net_seq$p))
  if(vectorize) {
    o <- matrix(o,ncol=net_seq$p,nrow=net_seq$p)
  }
  
  return(o)
}

## Implements the backfitting algorithm
NetworkAR_backfitting <- function(net_seq,B,limit,drop_first=FALSE,tol=0.00000000001) {
  ## Read data for easier writing
  p <- net_seq$p
  n <- net_seq$n
  
  ## Drop initial step if asked
  if(drop_first) {
    cat("Drop first observation.\n")
    net_seq$nts <- net_seq$nts[2:(n+1)]
    net_seq$t <- net_seq$t[2:(n+1)]
    net_seq$ncov <- net_seq$ncov[2:(n+1)]
    net_seq$el <- net_seq$el[2:(n+1)]
    net_seq$n <- n-1
  }
  
  ## Prepare output
  c_seq <- vector(mode="list",length=B)
  h_seq <- vector(mode="list",length=B+1)
  
  ## Compute initial value for h
  cat("Compute initial value for h.\n")
  h_seq[[1]] <- matrix(.Call("initialh",net_seq$el,as.integer(net_seq$n),as.integer(net_seq$p)),ncol=net_seq$p,nrow=net_seq$p)
  
  ## Do iterations
  for(i in 1:B) {
    ## Optimize with respect to c
    c_seq[[i]] <- matrix(.Call("optimize_c",net_seq$el,net_seq$ncov,h_seq[[i]],as.integer(net_seq$n+1),as.integer(net_seq$p),limit,tol),ncol=p,nrow=p)
    
    ## Optimize with respect to h
    h_seq[[i+1]] <- matrix(.Call("optimize_h",net_seq$el,net_seq$ncov,c_seq[[i]],as.integer(net_seq$n+1),as.integer(net_seq$p),limit,tol),ncol=p,nrow=p)
    
    if(i>=2) {
      cat("In interation ",i,"/",B,": Max change in h: ",max(abs(h_seq[[i+1]]-h_seq[[i]]))," max chage in c: ",max(abs(c_seq[[i]]-c_seq[[i-1]])),".\n")
    }
  }
  
  return(list(h_seq=h_seq,c_seq=c_seq))
}