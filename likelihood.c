#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>

// Functions from support.c
double max(double x1,double x2);
double min(double x1,double x2);
double myabs(double x);
int get_index_from_edge_list(SEXP X,int i,int j,int t);
int get_index_from_edge_list_ncov(SEXP ncov,int i,int j,int t);

// Functions contained in this file
SEXP log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP combineR);
double* log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int combine);

SEXP ddc_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);
double* ddc_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);

SEXP ddc_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j);
double ddc_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j);

SEXP d2dc2_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j);
double* d2dc2_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);
double d2dc2_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j);

SEXP ddh_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);
double* ddh_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);

SEXP ddh_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j);
double ddh_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j);

SEXP d2dh2_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);
double* d2dh2_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR);

SEXP d2dh2_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j);
double d2dh2_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j);

SEXP log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP combineR)
{
  double* log_lik;
  int i,p;
  int combine;
  SEXP out;
  
  combine=INTEGER(combineR)[0];
  
  // Compute Likelihood
  log_lik=log_likelihood_internal(X,ncov,c,h,nR,pR,combine);
  
  // Output
  if(INTEGER(combineR)[0]==0)
  {
    p=INTEGER(pR)[0];
    out=PROTECT(allocVector(REALSXP,p*p));
    for(i=0;i<=p*p-1;i++)
      REAL(out)[i]=log_lik[i];
  }
  else
  {
    out=PROTECT(allocVector(REALSXP,1));
    REAL(out)[0]=*log_lik;
  }
  
  free(log_lik);
  UNPROTECT(1);
  
  return(out);
}

double* log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int combine)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2,ri3;
  int s1,s2,s3;
  int ind;
  double* log_lik;
  double log_lik_temp;
  double ca;
  int curr_one,prev_one;
  
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for output (combine=0 means that the likelihood is not added)
  if(combine==0)
  {
    log_lik=malloc(sizeof(double)*p*p);
    for(i=0;i<=p*p-1;i++)
      log_lik[i]=0;
  }
  else
  {
    log_lik=malloc(sizeof(double));
    *log_lik=0;
  }
  
  // Compute Likelihood
  ri1=0;
  ri2=0;
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri1=0;
    ri2=0;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
        // Check if edge is presently there.
        if(s1>0 && ri1<=s1-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t))[ri1]==i && INTEGER(VECTOR_ELT(X,t))[ri1+s1]==j)
          {
            // There is a one in the current time step
            curr_one=1;
            ri1++;
          }
          else
          {
            // No there is no one in the current time step
            curr_one=0;
          }
        }
        else
        {
          curr_one=0;
        }
        
        // Check if edge was previously there.
        if(s2>0 && ri2<=s2-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t-1))[ri2]==i && INTEGER(VECTOR_ELT(X,t-1))[ri2+s2]==j)
          {
            // There is a one in the previous time step
            prev_one=1;
            ri2++;
          }
          else
          {
            // No there is no one in the previous time step
            prev_one=0;
          }
        }
        else
        {
          prev_one=0;
        }
        
        // Find the number of common successors.
        if(s3>0 && ri3<=s3-1)
        {
          // There are pairs with common successors left
          if((int)REAL(VECTOR_ELT(ncov,t-1))[ri3]==i && (int)REAL(VECTOR_ELT(ncov,t-1))[ri3+s3]==j)
          {
            // There are common ancestors
            ca=REAL(VECTOR_ELT(ncov,t-1))[ri3+2*s3];
            ri3++;
          }
          else
          {
            ca=0;
          }
        }
        else
        {
          ca=0;
        }
        
        // Compute contribution to likelihood
        ind=i-1+(j-1)*p;
        log_lik_temp=0;
        if(curr_one==1)
        {
          if(prev_one==1)
          {
            // There is an edge now, there was one before.
            if(REAL(c)[ind]*ca<200)
              log_lik_temp=log_lik_temp+log(1+exp(REAL(c)[ind]*ca));
            else
              log_lik_temp=log_lik_temp+REAL(c)[ind]*ca;
          } else
          {
            // There is an edge now, there was none before.
            log_lik_temp=log_lik_temp+REAL(c)[ind]*ca;
          }
        }
        else
        {
          if(prev_one==1)
          {
            // There is no edge now, there was one before.
            log_lik_temp=log_lik_temp+REAL(h)[ind];
          } else
          {
            // There is no edge now, there was none before.
            log_lik_temp=log_lik_temp+log(1+exp(REAL(h)[ind]));
          }
        }
        if(REAL(c)[ind]*ca<200)
          log_lik_temp=log_lik_temp-log(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]));
        else
          log_lik_temp=log_lik_temp-REAL(c)[ind]*ca;
        
        // Save the contrbution appropriately
        if(combine==0)
          log_lik[ind]=log_lik[ind]+log_lik_temp;
        else
          *log_lik=*log_lik+log_lik_temp;
      }
    }
  }
  
  return(log_lik);
}


SEXP ddc_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  double* ddc_log_lik;
  SEXP out;
  int p;
  int i;
  
  p=INTEGER(pR)[0];
  
  // Compute derivative
  ddc_log_lik=ddc_log_likelihood_internal(X,ncov,c,h,nR,pR);
  
  // Output
  out=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(out)[i]=ddc_log_lik[i];
  
  free(ddc_log_lik);
  
  UNPROTECT(1);
  return(out);
  
}


double* ddc_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2,ri3;
  int s1,s2,s3;
  int ind;
  double* ddc_log_lik;
  double ca;
  int curr_one,prev_one;
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for output and initialize
  ddc_log_lik=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
    ddc_log_lik[i]=0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri1=0;
    ri2=0;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
        // Check if edge is presently there.
        if(s1>0 && ri1<=s1-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t))[ri1]==i && INTEGER(VECTOR_ELT(X,t))[ri1+s1]==j)
          {
            // There is a one in the current time step
            curr_one=1;
            ri1++;
          }
          else
          {
            // No there is no one in the current time step
            curr_one=0;
          }
        }
        else
        {
          curr_one=0;
        }
        
        // Check if edge was previously there.
        if(s2>0 && ri2<=s2-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t-1))[ri2]==i && INTEGER(VECTOR_ELT(X,t-1))[ri2+s2]==j)
          {
            // There is a one in the previous time step
            prev_one=1;
            ri2++;
          }
          else
          {
            // No there is no one in the previous time step
            prev_one=0;
          }
        }
        else
        {
          prev_one=0;
        }
        
        // Find the number of common successors.
        if(s3>0 && ri3<=s3-1)
        {
          // There are pairs with common successors left
          if((int)REAL(VECTOR_ELT(ncov,t-1))[ri3]==i && (int)REAL(VECTOR_ELT(ncov,t-1))[ri3+s3]==j)
          {
            // There are common ancestors
            ca=REAL(VECTOR_ELT(ncov,t-1))[ri3+2*s3];
            ri3++;
          }
          else
          {
            ca=0;
          }
        }
        else
        {
          ca=0;
        }
        
        ind=i-1+(j-1)*p;
        if(curr_one==1)
        {
          if(prev_one==1)
          {
            // There is an edge now, there was one before.
            ddc_log_lik[ind]=ddc_log_lik[ind]+ca*exp(REAL(c)[ind]*ca)/(1+exp(REAL(c)[ind]*ca));
            
          } else
          {
            // There is an edge now, there was none before.
            ddc_log_lik[ind]=ddc_log_lik[ind]+ca;
          }
        }
        ddc_log_lik[ind]=ddc_log_lik[ind]-ca*exp(REAL(c)[ind]*ca)/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]));
      }
    }
  }
  
  return(ddc_log_lik);
}

SEXP ddc_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j)
{
  double ddc_log_lik;
  SEXP out;
  
  // Compute derivative
  ddc_log_lik=ddc_log_likelihood_single_internal(X,ncov,c,h,nR,pR,INTEGER(i)[0],INTEGER(j)[0]);
  
  // Output
  out=PROTECT(allocVector(REALSXP,1));
  REAL(out)[0]=ddc_log_lik;
  
  UNPROTECT(1);
  return(out);
}


double ddc_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j)
{
  // Data
  int t;
  int n,p,s3;
  int* ind;
  int ca_ind,ind0;
  double ddc_log_lik;
  double ca;
  int curr_one,prev_one;
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for index list and initialize
  ind=malloc(sizeof(int)*(n+1));
  for(t=0;t<=n;t++)
    ind[t]=get_index_from_edge_list(X,i,j,t);
  ddc_log_lik=0;
  ind0=i-1+(j-1)*p;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    
    // Check if edge is presently there (currone=1 -> YES).
    if(ind[t]!=-1)
      curr_one=1;
    else
      curr_one=0;
    
    // Check if edge was previously there.
    if(ind[t-1]!=-1)
      prev_one=1;
    else
      prev_one=0;
    
    // Find the number of common successors.
    ca_ind=get_index_from_edge_list_ncov(ncov,i,j,t-1);
    
    if(ca_ind!=-1)
    {
      // There are common ancestors
      ca=REAL(VECTOR_ELT(ncov,t-1))[ca_ind+2*s3];
    }
    else
      ca=0;
    
    if(curr_one==1)
    {
      if(prev_one==1)
      {
        // There is an edge now, there was one before.
        if(REAL(c)[ind0]*ca<200)
          ddc_log_lik=ddc_log_lik+ca*exp(REAL(c)[ind0]*ca)/(1+exp(REAL(c)[ind0]*ca));
        else
          ddc_log_lik=ddc_log_lik+ca;
      } else
      {
        // There is an edge now, there was none before.
        ddc_log_lik=ddc_log_lik+ca;
      }
    }
    if(REAL(c)[ind0]*ca<200)
      ddc_log_lik=ddc_log_lik-ca*exp(REAL(c)[ind0]*ca)/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]));
    else
      ddc_log_lik=ddc_log_lik-ca;
//    Rprintf("Edge (%d,%d), t=%d current derivative value: %f, ca=%f, c=%f.\n",i,j,t,ddc_log_lik,ca,REAL(c)[ind0]);
  }
  
  return(ddc_log_lik);
}

SEXP ddh_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  double* ddh_log_lik;
  SEXP out;
  int p;
  int i;
  
  p=INTEGER(pR)[0];
  
  // Compute derivative
  ddh_log_lik=ddh_log_likelihood_internal(X,ncov,c,h,nR,pR);
  
  // Output
  out=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(out)[i]=ddh_log_lik[i];
  
  free(ddh_log_lik);
  
  UNPROTECT(1);
  return(out);
}


double* ddh_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2,ri3;
  int s1,s2,s3;
  int ind;
  double* ddh_log_lik;
  double ca;
  int curr_one,prev_one;
  
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for output and initialize
  ddh_log_lik=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
    ddh_log_lik[i]=0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri1=0;
    ri2=0;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
        // Check if edge is presently there.
        if(s1>0 && ri1<=s1-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t))[ri1]==i && INTEGER(VECTOR_ELT(X,t))[ri1+s1]==j)
          {
            // There is a one in the current time step
            curr_one=1;
            ri1++;
          }
          else
          {
            // No there is no one in the current time step
            curr_one=0;
          }
        }
        else
        {
          curr_one=0;
        }
        
        // Check if edge was previously there.
        if(s2>0 && ri2<=s2-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t-1))[ri2]==i && INTEGER(VECTOR_ELT(X,t-1))[ri2+s2]==j)
          {
            // There is a one in the previous time step
            prev_one=1;
            ri2++;
          }
          else
          {
            // No there is no one in the previous time step
            prev_one=0;
          }
        }
        else
        {
          prev_one=0;
        }
        
        // Find the number of common successors.
        if(s3>0 && ri3<=s3-1)
        {
          // There are pairs with common successors left
          if((int)REAL(VECTOR_ELT(ncov,t-1))[ri3]==i && (int)REAL(VECTOR_ELT(ncov,t-1))[ri3+s3]==j)
          {
            // There are common ancestors
            ca=REAL(VECTOR_ELT(ncov,t-1))[ri3+2*s3];
            ri3++;
          }
          else
          {
            ca=0;
          }
        }
        else
        {
          ca=0;
        }
        
        ind=i-1+(j-1)*p;
        if(curr_one==0)
        {
          if(prev_one==1)
          {
            // There is no edge now, there was one before.
            ddh_log_lik[ind]=ddh_log_lik[ind]+1;
          } else
          {
            // There is no edge now, there was none before.
            ddh_log_lik[ind]=ddh_log_lik[ind]+exp(REAL(h)[ind])/(1+exp(REAL(h)[ind]));
          }
        }
        ddh_log_lik[ind]=ddh_log_lik[ind]-exp(REAL(h)[ind])/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]));
      }
    }
  }
  
  return(ddh_log_lik);
}

SEXP ddh_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j)
{
  double ddh_log_lik;
  SEXP out;
  
  // Compute derivative
  ddh_log_lik=ddh_log_likelihood_single_internal(X,ncov,c,h,nR,pR,INTEGER(i)[0],INTEGER(j)[0]);
  
  // Output
  out=PROTECT(allocVector(REALSXP,1));
  REAL(out)[0]=ddh_log_lik;
  
  UNPROTECT(1);
  return(out);
}

double ddh_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j)
{
  // Data
  int t;
  int n,p,s3;
  int* ind;
  int ca_ind,ind0;
  double ddh_log_lik;
  double ca;
  int curr_one,prev_one;
  
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for index list and initialize
  ind=malloc(sizeof(int)*(n+1));
  for(t=0;t<=n;t++)
    ind[t]=get_index_from_edge_list(X,i,j,t);
  ddh_log_lik=0.0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    
    // Check if edge is presently there (currone=1 -> YES).
    if(ind[t]!=-1)
      curr_one=1;
    else
      curr_one=0;
    
    // Check if edge was previously there.
    if(ind[t-1]!=-1)
      prev_one=1;
    else
      prev_one=0;
    
    // Find the number of common successors.
    ca_ind=get_index_from_edge_list_ncov(ncov,i,j,t-1);
    
    if(ca_ind!=-1)
    {
      // There are common ancestors
      ca=REAL(VECTOR_ELT(ncov,t-1))[ca_ind+2*s3];
    }
    else
      ca=0;
    
    
    ind0=i-1+(j-1)*p;
    if(curr_one==0)
    {
      if(prev_one==1)
      {
        // There is no edge now, there was one before.
        ddh_log_lik=ddh_log_lik+1.0;
      } else
      {
        // There is no edge now, there was none before.
        ddh_log_lik=ddh_log_lik+exp(REAL(h)[ind0])/(1+exp(REAL(h)[ind0]));
      }
    }
    ddh_log_lik=ddh_log_lik-exp(REAL(h)[ind0])/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]));
  }
  
  return(ddh_log_lik);
}

SEXP d2dc2_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  double* d2dc2_log_lik;
  SEXP out;
  int p;
  int i;
  
  p=INTEGER(pR)[0];
  
  // Compute derivative
  d2dc2_log_lik=d2dc2_log_likelihood_internal(X,ncov,c,h,nR,pR);
  
  // Output
  out=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(out)[i]=d2dc2_log_lik[i];
  
  free(d2dc2_log_lik);
  
  UNPROTECT(1);
  return(out);
}

double* d2dc2_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2,ri3;
  int s1,s2,s3;
  int ind;
  double* d2dc2_log_lik;
  double ca;
  int curr_one,prev_one;
  
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for output and initialize
  d2dc2_log_lik=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
    d2dc2_log_lik[i]=0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri1=0;
    ri2=0;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
        // Check if edge is presently there.
        if(s1>0 && ri1<=s1-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t))[ri1]==i && INTEGER(VECTOR_ELT(X,t))[ri1+s1]==j)
          {
            // There is a one in the current time step
            curr_one=1;
            ri1++;
          }
          else
          {
            // No there is no one in the current time step
            curr_one=0;
          }
        }
        else
        {
          curr_one=0;
        }
        
        // Check if edge was previously there.
        if(s2>0 && ri2<=s2-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t-1))[ri2]==i && INTEGER(VECTOR_ELT(X,t-1))[ri2+s2]==j)
          {
            // There is a one in the previous time step
            prev_one=1;
            ri2++;
          }
          else
          {
            // No there is no one in the previous time step
            prev_one=0;
          }
        }
        else
        {
          prev_one=0;
        }
        
        // Find the number of common successors.
        if(s3>0 && ri3<=s3-1)
        {
          // There are pairs with common successors left
          if((int)REAL(VECTOR_ELT(ncov,t-1))[ri3]==i && (int)REAL(VECTOR_ELT(ncov,t-1))[ri3+s3]==j)
          {
            // There are common ancestors
            ca=REAL(VECTOR_ELT(ncov,t-1))[ri3+2*s3];
            ri3++;
          }
          else
          {
            ca=0;
          }
        }
        else
        {
          ca=0;
        }
        
        ind=i-1+(j-1)*p;
        if(curr_one==1)
        {
          if(prev_one==1)
          {
            // There is an edge now, there was one before.
            d2dc2_log_lik[ind]=d2dc2_log_lik[ind]+ca*ca*exp(REAL(c)[ind]*ca)/(1+exp(REAL(c)[ind]*ca))/(1+exp(REAL(c)[ind]*ca));
            
          }
        }
        d2dc2_log_lik[ind]=d2dc2_log_lik[ind]-(ca*ca*exp(REAL(c)[ind]*ca)*(1+exp(REAL(h)[ind])))/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]))/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]));
      }
    }
  }
  
  return(d2dc2_log_lik);
}

SEXP d2dc2_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j)
{
  double d2dc2_log_lik;
  SEXP out;
  
  // Compute derivative
  d2dc2_log_lik=d2dc2_log_likelihood_single_internal(X,ncov,c,h,nR,pR,INTEGER(i)[0],INTEGER(j)[0]);
  
  // Output
  out=PROTECT(allocVector(REALSXP,1));
  REAL(out)[0]=d2dc2_log_lik;
  
  UNPROTECT(1);
  return(out);
}


double d2dc2_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j)
{
  // Data
  int t;
  int n,p,s3;
  int* ind;
  int ca_ind,ind0;
  double d2dc2_log_lik;
  double ca;
  int curr_one,prev_one;
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for index list and initialize
  ind=malloc(sizeof(int)*(n+1));
  for(t=0;t<=n;t++)
    ind[t]=get_index_from_edge_list(X,i,j,t);
  d2dc2_log_lik=0;
  
  // Compute second derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    
    // Check if edge is presently there (currone=1 -> YES).
    if(ind[t]!=-1)
      curr_one=1;
    else
      curr_one=0;
    
    // Check if edge was previously there.
    if(ind[t-1]!=-1)
      prev_one=1;
    else
      prev_one=0;
    
    // Find the number of common successors.
    ca_ind=get_index_from_edge_list_ncov(ncov,i,j,t-1);
    
    if(ca_ind!=-1)
    {
      // There are common ancestors
      ca=REAL(VECTOR_ELT(ncov,t-1))[ca_ind+2*s3];
    }
    else
      ca=0;
    
    ind0=i-1+(j-1)*p;
    if(curr_one==1)
    {
      if(prev_one==1)
      {
        // There is an edge now, there was one before.
        if(REAL(c)[ind0]*ca<400)
          d2dc2_log_lik=d2dc2_log_lik+ca*ca*exp(REAL(c)[ind0]*ca)/(1+exp(REAL(c)[ind0]*ca))/(1+exp(REAL(c)[ind0]*ca));
        else
          d2dc2_log_lik=d2dc2_log_lik;
      }
    }
    if(REAL(c)[ind0]*ca<400)
      d2dc2_log_lik=d2dc2_log_lik-(ca*ca*exp(REAL(c)[ind0]*ca)*(1+exp(REAL(h)[ind0])))/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]))/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]));
    else
       d2dc2_log_lik=d2dc2_log_lik;
  }
  
  return(d2dc2_log_lik);
}


SEXP d2dh2_log_likelihood(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  double* d2dh2_log_lik;
  SEXP out;
  int p;
  int i;
  
  p=INTEGER(pR)[0];
  
  // Compute derivative
  d2dh2_log_lik=d2dh2_log_likelihood_internal(X,ncov,c,h,nR,pR);
  
  // Output
  out=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(out)[i]=d2dh2_log_lik[i];
  
  free(d2dh2_log_lik);
  
  UNPROTECT(1);
  return(out);
}

double* d2dh2_log_likelihood_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2,ri3;
  int s1,s2,s3;
  int ind;
  double* d2dh2_log_lik;
  double ca;
  int curr_one,prev_one;
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for output and initialize
  d2dh2_log_lik=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
    d2dh2_log_lik[i]=0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri1=0;
    ri2=0;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
        // Check if edge is presently there.
        if(s1>0 && ri1<=s1-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t))[ri1]==i && INTEGER(VECTOR_ELT(X,t))[ri1+s1]==j)
          {
            // There is a one in the current time step
            curr_one=1;
            ri1++;
          }
          else
          {
            // No there is no one in the current time step
            curr_one=0;
          }
        }
        else
        {
          curr_one=0;
        }
        
        // Check if edge was previously there.
        if(s2>0 && ri2<=s2-1)
        {
          // There are edges left to check
          if(INTEGER(VECTOR_ELT(X,t-1))[ri2]==i && INTEGER(VECTOR_ELT(X,t-1))[ri2+s2]==j)
          {
            // There is a one in the previous time step
            prev_one=1;
            ri2++;
          }
          else
          {
            // No there is no one in the previous time step
            prev_one=0;
          }
        }
        else
        {
          prev_one=0;
        }
        
        // Find the number of common successors.
        if(s3>0 && ri3<=s3-1)
        {
          // There are pairs with common successors left
          if((int)REAL(VECTOR_ELT(ncov,t-1))[ri3]==i && (int)REAL(VECTOR_ELT(ncov,t-1))[ri3+s3]==j)
          {
            // There are common ancestors
            ca=REAL(VECTOR_ELT(ncov,t-1))[ri3+2*s3];
            ri3++;
          }
          else
          {
            ca=0;
          }
        }
        else
        {
          ca=0;
        }
        
        ind=i-1+(j-1)*p;
        if(curr_one==0)
        {
          if(prev_one==0)
          {
            // There is no edge now, there was none before.
            d2dh2_log_lik[ind]=d2dh2_log_lik[ind]+exp(REAL(h)[ind])/(1+exp(REAL(h)[ind]))/(1+exp(REAL(h)[ind]));
          }
        }
        d2dh2_log_lik[ind]=d2dh2_log_lik[ind]-exp(REAL(h)[ind])*(1+exp(REAL(c)[ind]*ca))/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]))/(1+exp(REAL(c)[ind]*ca)+exp(REAL(h)[ind]));
      }
    }
  }
  
  return(d2dh2_log_lik);
}

SEXP d2dh2_log_likelihood_single(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,SEXP i,SEXP j)
{
  double d2dh2_log_lik;
  SEXP out;
  
  // Compute derivative
  d2dh2_log_lik=d2dh2_log_likelihood_single_internal(X,ncov,c,h,nR,pR,INTEGER(i)[0],INTEGER(j)[0]);
  
  // Output
  out=PROTECT(allocVector(REALSXP,1));
  REAL(out)[0]=d2dh2_log_lik;
  
  UNPROTECT(1);
  return(out);
}

double d2dh2_log_likelihood_single_internal(SEXP X,SEXP ncov,SEXP c,SEXP h,SEXP nR,SEXP pR,int i,int j)
{
  // Data
  int t;
  int n,p,s3;
  int* ind;
  int ca_ind,ind0;
  double d2dh2_log_lik;
  double ca;
  int curr_one,prev_one;
  
  
  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate memory for index list and initialize
  ind=malloc(sizeof(int)*(n+1));
  for(t=0;t<=n;t++)
    ind[t]=get_index_from_edge_list(X,i,j,t);
  d2dh2_log_lik=0;
  
  // Compute derivative of Likelihood
  for(t=1;t<=n;t++)
  {
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    
    // Check if edge is presently there (currone=1 -> YES).
    if(ind[t]!=-1)
      curr_one=1;
    else
      curr_one=0;
    
    // Check if edge was previously there.
    if(ind[t-1]!=-1)
      prev_one=1;
    else
      prev_one=0;
    
    // Find the number of common successors.
    ca_ind=get_index_from_edge_list_ncov(ncov,i,j,t-1);
    
    if(ca_ind!=-1)
    {
      // There are common ancestors
      ca=REAL(VECTOR_ELT(ncov,t-1))[ca_ind+2*s3];
    }
    else
      ca=0;
    
    ind0=i-1+(j-1)*p;
    if(curr_one==0)
      if(prev_one==0)
        d2dh2_log_lik=d2dh2_log_lik+exp(REAL(h)[ind0])/(1+exp(REAL(h)[ind0]))/(1+exp(REAL(h)[ind0]));
    d2dh2_log_lik=d2dh2_log_lik-exp(REAL(h)[ind0])*(1+exp(REAL(c)[ind0]*ca))/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]))/(1+exp(REAL(c)[ind0]*ca)+exp(REAL(h)[ind0]));
  }
  
  return(d2dh2_log_lik);
}
