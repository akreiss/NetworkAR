#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>


double max(double x1,double x2)
{
  if(x1>=x2)
    return(x1);
  else
    return(x2);
}
double min(double x1,double x2)
{
  if(x1>=x2)
    return(x2);
  else
    return(x1);
}
double myabs(double x)
{
  if(x>=0)
    return(x);
  else
    return(-x);
}

int get_index_from_edge_list(SEXP X,int i,int j,int t)
{
  int ub,lb,K;
  int found_it;
  int k;
  
  // Get number of edges
  K=LENGTH(VECTOR_ELT(X,t))/2;
  
  // Return if there are no edges
  if(K==0)
    return(-1);
  
  ub=K-1;
  lb=0;
  
  // Check if boundaries are what we look for
  if(INTEGER(VECTOR_ELT(X,t))[lb]==i && INTEGER(VECTOR_ELT(X,t))[lb+K]==j)
    return(lb);
  if(INTEGER(VECTOR_ELT(X,t))[ub]==i && INTEGER(VECTOR_ELT(X,t))[ub+K]==j)
    return(ub);
  
  // If not, search
  found_it=0;
  while(found_it==0) {
    if((ub-lb)%2==0)
      k=lb+(ub-lb)/2;
    else
      k=lb+(ub-lb-1)/2;
    
    if(ub-lb<=1)
      found_it=-1;
    
    if(INTEGER(VECTOR_ELT(X,t))[k+K]==j)
    {
      if(INTEGER(VECTOR_ELT(X,t))[k]==i)
        found_it=1;
      else if(INTEGER(VECTOR_ELT(X,t))[k]<i)
        lb=k;
      else
        ub=k;
    }
    else if(INTEGER(VECTOR_ELT(X,t))[k+K]<j)
      lb=k;
    else
      ub=k;
  }
  
  if(found_it==1)
    return(k);
  else
    return(-1);
}

int get_index_from_edge_list_ncov(SEXP ncov,int i,int j,int t)
{
  int ub,lb,K;
  int found_it;
  int k;
  
  // Get number of edges
  K=LENGTH(VECTOR_ELT(ncov,t))/3;
  
  // Return if there are no edges
  if(K==0)
    return(-1);
  
  ub=K-1;
  lb=0;
  
  // Check if boundaries are what we look for
  if((int)REAL(VECTOR_ELT(ncov,t))[lb]==i && (int)REAL(VECTOR_ELT(ncov,t))[lb+K]==j)
    return(lb);
  if((int)REAL(VECTOR_ELT(ncov,t))[ub]==i && (int)REAL(VECTOR_ELT(ncov,t))[ub+K]==j)
    return(ub);
  
  // If not, search
  found_it=0;
  while(found_it==0) {
    if((ub-lb)%2==0)
      k=lb+(ub-lb)/2;
    else
      k=lb+(ub-lb-1)/2;
    
    if(ub-lb<=1)
      found_it=-1;
    
    if((int)REAL(VECTOR_ELT(ncov,t))[k+K]==j)
    {
      if((int)REAL(VECTOR_ELT(ncov,t))[k]==i)
        found_it=1;
      else if((int)REAL(VECTOR_ELT(ncov,t))[k]<i)
        lb=k;
      else
        ub=k;
    }
    else if((int)REAL(VECTOR_ELT(ncov,t))[k+K]<j)
      lb=k;
    else
      ub=k;
  }
  
  if(found_it==1)
    return(k);
  else
    return(-1);
}