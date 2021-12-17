#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>

// Functions contained in likelihood.c
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


// Functions from support.c
double max(double x1,double x2);
double min(double x1,double x2);
double myabs(double x);
int get_index_from_edge_list(SEXP X,int i,int j,int t);
int get_index_from_edge_list_ncov(SEXP ncov,int i,int j,int t);

// Functions in this file
SEXP find_lower_bound_R(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR);
SEXP find_lower_bound(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR);

SEXP find_upper_bound_R(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,SEXP limitR);
SEXP find_upper_bound(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,double limit);

SEXP find_lower_bound_hR(SEXP X,SEXP ncov,SEXP nR,SEXP pR,SEXP limitR);
SEXP find_lower_bound_h(SEXP X,SEXP ncov,SEXP nR,SEXP pR,double limit);

SEXP find_upper_bound_h(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,double limit);
SEXP find_upper_bound_hR(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,SEXP limitR);
  
SEXP optimize_c(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,SEXP limitR,SEXP tolR);
SEXP optimize_h(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,SEXP limitR,SEXP tolR);

SEXP initialh(SEXP X,SEXP nR,SEXP pR);





SEXP initialh(SEXP X,SEXP nR,SEXP pR)
{
  // Data
  int t,i,j;
  int n,p;
  int ri1,ri2;
  int s1,s2;
  int ind;
  int curr_one,prev_one;
  double* a1;
  double* a2;
  double* a3;
  double* a4;
  double q;
  SEXP out;

  // Read information from input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;

  // Allocate memory for output and initialize
  a1=malloc(sizeof(double)*p*p);
  a2=malloc(sizeof(double)*p*p);
  a3=malloc(sizeof(double)*p*p);
  a4=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
  {
    a1[i]=0.0001;
    a2[i]=0.0001;
    a3[i]=0.0001;
    a4[i]=0.0001;
  }


  // Compute four terms in the fraction a1,...,a4
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    ri1=0;
    ri2=0;

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

        // Compute the terrms
        ind=i-1+(j-1)*p;
        if(prev_one==1)
        {
          // There was an edge before
          a4[ind]=a4[ind]+1;

          if(curr_one==0)
          {
            // There is no edge now, there was one before.
            a1[ind]=a1[ind]+1;
          }
        }
        else
        {
          // There was no edge before
          a2[ind]=a2[ind]+1;

          if(curr_one==0)
          {
            // There is no edge now, there was none before.
            a3[ind]=a3[ind]+1;
          }
        }
      }
    }
  }

  // Output
  out=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
  {
    q=a1[i]*a2[i]/a3[i]/a4[i];
    REAL(out)[i]=log(q/max(0.0001,1-q));
  }
  UNPROTECT(1);

  return(out);
}



// Find lower and upper bounds for optimizing c
SEXP find_upper_bound_R(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,SEXP limitR)
{
  double limit;
  SEXP out;
  limit=REAL(limitR)[0];
  out=find_upper_bound(X,ncov,h,nR,pR,limit);
  UNPROTECT(1);
  return(out);
}
SEXP find_lower_bound_R(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR)
{
  SEXP out;
  out=find_lower_bound(X,ncov,h,nR,pR);
  UNPROTECT(1);
  return(out);
}

SEXP find_upper_bound(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,double limit)
{
  SEXP upper_bound;
  double* derivative;
  int i,j;
  int ind,p;
  int pos;
  double nv;

  p=INTEGER(pR)[0];

  // Allocate Memory
  upper_bound=PROTECT(allocVector(REALSXP,p*p));

  // Set initial values of upper bound
  for(i=0;i<=p*p-1;i++)
    REAL(upper_bound)[i]=min(0.5*log(1+exp(REAL(h)[i])),limit);

  // Check if refinement is necessary
  derivative=ddc_log_likelihood_internal(X,ncov,upper_bound,h,nR,pR);
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      ind=(i-1)+(j-1)*p;
      if((derivative[ind]<0) | (REAL(upper_bound)[ind]>=limit))
        pos=0;
      else
        pos=1;

      while(pos==1)
      {
        REAL(upper_bound)[ind]=min(1.5*(REAL(upper_bound)[ind]+0.1),limit);
        nv=ddc_log_likelihood_single_internal(X,ncov,upper_bound,h,nR,pR,i,j);
        if((nv<0) | (REAL(upper_bound)[ind]>=limit))
          pos=0;
      }
    }
  }

  free(derivative);

  return(upper_bound);
}

SEXP find_lower_bound(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR)
{
  int p,i,j,ind;
  SEXP lower_bound;
  double xm,lb,ub,nv;
  double* second_derivative;
  double tol;

  tol=0.00000000001;
  p=INTEGER(pR)[0];

  // Allocate Memory
  lower_bound=PROTECT(allocVector(REALSXP,p*p));

  // Set initial lower bound
  for(i=0;i<=p*p-1;i++)
    REAL(lower_bound)[i]=0;

  // Check if refinement is necessary
  second_derivative=d2dc2_log_likelihood_internal(X,ncov,lower_bound,h,nR,pR);
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      ind=(i-1)+(j-1)*p;

      if(second_derivative[ind]<0)
      {
        lb=0;
        ub=0;
      }
      else
      {
        lb=0;
        ub=0.5*log(1+exp(REAL(h)[ind]));
      }

      // Find the zero of the second derivative
      while(ub-lb>tol)
      {
        xm=0.5*(lb+ub);
        REAL(lower_bound)[ind]=xm;
        nv=d2dc2_log_likelihood_single_internal(X,ncov,lower_bound,h,nR,pR,i,j);
        if(nv>0)
          lb=xm;
        else
          ub=xm;
      }
    }
  }
  free(second_derivative);
  return(lower_bound);
}

// Find upper and lower bounds for optimizing h
SEXP find_lower_bound_hR(SEXP X,SEXP ncov,SEXP nR,SEXP pR,SEXP limitR)
{
  SEXP out;
  double limit;
  
  limit=REAL(limitR)[0];
  
  out=find_lower_bound_h(X,ncov,nR,pR,limit);
  UNPROTECT(1);
  return(out);
}
SEXP find_upper_bound_hR(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,SEXP limitR)
{
  SEXP out;
  double limit;
  
  limit=REAL(limitR)[0];
  
  out=find_upper_bound_h(X,ncov,c,nR,pR,limit);
  UNPROTECT(1);
  return(out);
}

SEXP find_lower_bound_h(SEXP X,SEXP ncov,SEXP nR,SEXP pR,double limit)
{
  double* xone;
  double* xtwo;
  SEXP lower_bound;
  int n,p,s1,s2,ri1,ri2,i,j,t;
  int curr_one,prev_one;
  
  // Read Information from input
  n=INTEGER(nR)[0]-1;
  p=INTEGER(pR)[0];
  lower_bound=PROTECT(allocVector(REALSXP,p*p));
  
  // Allocate memory
  xone=malloc(sizeof(double)*p*p);
  xtwo=malloc(sizeof(double)*p*p);
  for(i=0;i<=p*p-1;i++)
  {
    xone[i]=0.0;
    xtwo[i]=0.0;
  }
  
  // Compute xone and xtwo
  ri1=0;
  ri2=0;
  for(t=1;t<=n;t++)
  {
    s1=LENGTH(VECTOR_ELT(X,t))/2;
    s2=LENGTH(VECTOR_ELT(X,t-1))/2;
    ri1=0;
    ri2=0;

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
            curr_one=0;
        }
        else
          curr_one=0;
        
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
            prev_one=0;
        }
        else
          prev_one=0;
        
        // Increase xone and xtwo appropriately
        if(prev_one==1)
        {
          xone[(i-1)+(j-1)*p]=xone[(i-1)+(j-1)*p]-1.0;
          if(curr_one==0)
            xtwo[(i-1)+(j-1)*p]=xtwo[(i-1)+(j-1)*p]+1.0;
        }
        else
        {
          if(curr_one==1)
            xone[(i-1)+(j-1)*p]=xone[(i-1)+(j-1)*p]-1.0;
        }
      }
    }
  }
  
  // Compute Lower Bound
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      if(xone[(i-1)+(j-1)*p]+xtwo[(i-1)+(j-1)*p]==0)
        REAL(lower_bound)[(i-1)+(j-1)*p]=limit;
      else if(xtwo[(i-1)+(j-1)*p]==0)
        REAL(lower_bound)[(i-1)+(j-1)*p]=-limit;
      else
        REAL(lower_bound)[(i-1)+(j-1)*p]=log(-xtwo[(i-1)+(j-1)*p]/(xone[(i-1)+(j-1)*p]+xtwo[(i-1)+(j-1)*p]));
    }
  }
  
  free(xone);
  free(xtwo);
  
  return(lower_bound);
}

SEXP find_upper_bound_h(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,double limit)
{
  SEXP upper_bound;
  int n,p,i,j,t,s3,ri3,ind,pos;
  double ca,prop_ub;
  double* derivative;
  double nv;
  
  // Get Information From Input
  p=INTEGER(pR)[0];
  n=INTEGER(nR)[0]-1;
  
  // Allocate Memory and initialize
  upper_bound=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(upper_bound)[i]=0.0;
  
  
  // Set initial values of upper bound
  for(t=1;t<=n;t++)
  {
    s3=LENGTH(VECTOR_ELT(ncov,t-1))/3;
    ri3=0;
    
    for(j=1;j<=p;j++)
    {
      for(i=1;i<=p;i++)
      {
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
            ca=0.0;
        }
        else
          ca=0.0;
        
        // Set proposed upper bound
        if(REAL(c)[(i-1)+(j-1)*p]*ca<200)
          prop_ub=0.5*log(1.0+exp(REAL(c)[(i-1)+(j-1)*p]*ca));
        else
          prop_ub=0.5*REAL(c)[(i-1)+(j-1)*p]*ca;
        
        // Use the larger value
        REAL(upper_bound)[(i-1)+(j-1)*p]=max(REAL(upper_bound)[(i-1)+(j-1)*p],prop_ub);
        
//        Rprintf("Find upper bound for h: Edge (%d,%d), t=%d, current value=%f, proposed value=%f, ca=%f, c=%f.\n",i,j,t,REAL(upper_bound)[(i-1)+(j-1)*p],prop_ub,ca,REAL(c)[(i-1)+(j-1)*p]);
      }
    }
  }

  // Check if increasing is necessary
  derivative=ddh_log_likelihood_internal(X,ncov,c,upper_bound,nR,pR);
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      ind=(i-1)+(j-1)*p;
      if((derivative[ind]<0) | (REAL(upper_bound)[ind]>=limit))
        pos=0;
      else
        pos=1;
      
      while(pos==1)
      {
        REAL(upper_bound)[ind]=min(1.5*(REAL(upper_bound)[ind]+0.1),limit);
        nv=ddh_log_likelihood_single_internal(X,ncov,c,upper_bound,nR,pR,i,j);

        if((nv<0) | (REAL(upper_bound)[ind]>=limit))
          pos=0;
      }
      
      if(REAL(upper_bound)[ind]>1000.0)
        Rprintf("WARNING: Extremely Large upper bound for edge (%d,%d): %f\n",i,j,REAL(upper_bound)[ind]);
    }
  }
  
  free(derivative);
  
  return(upper_bound);
}

// Optimize c or h when the other is given
SEXP optimize_c(SEXP X,SEXP ncov,SEXP h,SEXP nR,SEXP pR,SEXP limitR,SEXP tolR)
{
  SEXP lb;
  SEXP ub;
  double limit;
  double tol;
  SEXP x0;
  SEXP xnull;
  double* likelihood;
  double* likelihood_null;
  double* derivatives;
  double derivative;
  double second_derivative;
  int p;
  double newton;
  int i,j,ind;

  limit=REAL(limitR)[0];
  tol=REAL(tolR)[0];
  p=INTEGER(pR)[0];

  // Compute lower and upper bounds on the maxima (the concave interval)
  lb=find_lower_bound(X,ncov,h,nR,pR);
  ub=find_upper_bound(X,ncov,h,nR,pR,limit);

  // Check if derivative is only decreasing on concave interval
  derivatives=ddc_log_likelihood_internal(X,ncov,lb,h,nR,pR);
  for(i=0;i<=p*p-1;i++)
    if(derivatives[i]<0)
      REAL(ub)[i]=REAL(lb)[i];
  free(derivatives);

  // Check if derivative is only increasing on concave interval
  derivatives=ddc_log_likelihood_internal(X,ncov,ub,h,nR,pR);
  for(i=0;i<=p*p-1;i++)
    if(derivatives[i]>0)
      REAL(lb)[i]=REAL(ub)[i];
  free(derivatives);

  // Compute Next Candidate point
  x0=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(x0)[i]=0.5*(REAL(lb)[i]+REAL(ub)[i]);

  // Find the maximisers in the concave interval
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      ind=(i-1)+(j-1)*p;
//      Rprintf("c: Optimize edge (%d,%d).\n",i,j);
      // Check if convergence is reached
      while(REAL(ub)[ind]-REAL(lb)[ind]>=tol)
      {
        // Tolerance target has not been reached yet.
        // Evaluate the derivatives at their current spot
        derivative=ddc_log_likelihood_single_internal(X,ncov,x0,h,nR,pR,i,j);
        second_derivative=d2dc2_log_likelihood_single_internal(X,ncov,x0,h,nR,pR,i,j);
        
//        Rprintf("c: Edge (%d,%d), lb=%f, ub=%f, x0=%f, h=%f, deriv=%f.\n",i,j, REAL(lb)[ind],REAL(ub)[ind],REAL(x0)[ind],REAL(h)[ind],derivative);
        if(!(derivative==derivative)) {
          Rprintf("Optimize c: SEVERE ERROR!!!! SEVERE ERROR!!!!\n");
          return(x0);
        }
        if(!(second_derivative==second_derivative)) {
          Rprintf("Optimize c: SEVERE ERROR!!!! SEVERE ERROR!!!! Second Derivative equals NaN for edge (%d,%d)\n",i,j);
          return(x0);
        }
        
        // Update boundaries
        if(derivative>0)
          REAL(lb)[ind]=REAL(x0)[ind];
        if(derivative<0)
          REAL(ub)[ind]=REAL(x0)[ind];
        if(derivative==0)
        {
          REAL(lb)[ind]=REAL(x0)[ind];
          REAL(ub)[ind]=REAL(x0)[ind];
        }

         // Compute Newton Step (unless zero second derivative)
        if(second_derivative!=0)
          newton=REAL(x0)[ind]-derivative/second_derivative;
        else
          newton=0.5*(REAL(lb)[ind]+REAL(ub)[ind]);

        // Check if Newton step falls out of the interval or gives too little improvement, if so do bisection
        if((newton<=REAL(lb)[ind]) | (newton>=REAL(ub)[ind])  | (myabs(newton-REAL(x0)[ind])<=0.001*(REAL(ub)[ind]-REAL(lb)[ind])))
          REAL(x0)[ind]=0.5*(REAL(lb)[ind]+REAL(ub)[ind]);
        else
          REAL(x0)[ind]=newton;
      }
    }
  }

  // Check if the boundary yields a better smaller likelihood
  xnull=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<p*p-1;i++)
    REAL(xnull)[i]=0.0;
  likelihood=log_likelihood_internal(X,ncov,x0,h,nR,pR,0);
  likelihood_null=log_likelihood_internal(X,ncov,xnull,h,nR,pR,0);
  for(i=0;i<=p*p-1;i++)
    if(likelihood_null[i]>likelihood[i])
      REAL(x0)[i]=0.0;
  free(likelihood);
  free(likelihood_null);


  UNPROTECT(4);
  return(x0);
}

SEXP optimize_h(SEXP X,SEXP ncov,SEXP c,SEXP nR,SEXP pR,SEXP limitR,SEXP tolR)
{
  SEXP lb;
  SEXP ub;
  double limit;
  double tol;
  SEXP x0;
  double* derivatives_lb;
  double* derivatives_ub;
  double derivative;
  double second_derivative;
  int p;
  double newton;
  int i,j,ind;
  
//  Rprintf("Hallo1\n");
  
  limit=REAL(limitR)[0];
  tol=REAL(tolR)[0];
  p=INTEGER(pR)[0];
  
//  Rprintf("Hallo2\n");
  
  // Compute lower and upper bounds on the maxima
  lb=find_lower_bound_h(X,ncov,nR,pR,limit);
  ub=find_upper_bound_h(X,ncov,c,nR,pR,limit);
  
//  Rprintf("Hallo3\n");
  
  // Check if derivative has the same sign and act accordingly
  derivatives_lb=ddh_log_likelihood_internal(X,ncov,c,lb,nR,pR);
  derivatives_ub=ddh_log_likelihood_internal(X,ncov,c,ub,nR,pR);
  
//  Rprintf("Hallo4\n"); 
  
  for(i=0;i<=p*p-1;i++)
  {
    if((derivatives_lb[i]<0) & (derivatives_ub[i]<0))
      REAL(ub)[i]=REAL(lb)[i];
    
    if((derivatives_lb[i]>0) & (derivatives_ub[i]>0))
      REAL(lb)[i]=REAL(ub)[i];
  }
  
//  Rprintf("Hallo5\n");
  
  free(derivatives_lb);
  free(derivatives_ub);
  
//  Rprintf("Hallo6\n");
    
  // Compute Next Candidate points
  x0=PROTECT(allocVector(REALSXP,p*p));
  for(i=0;i<=p*p-1;i++)
    REAL(x0)[i]=0.5*(REAL(lb)[i]+REAL(ub)[i]);
  
//  Rprintf("Hallo7\n");
  
  // Find the maximisers
  for(i=1;i<=p;i++)
  {
    for(j=1;j<=p;j++)
    {
      ind=(i-1)+(j-1)*p;
      // Check if convergence is reached
      while(REAL(ub)[ind]-REAL(lb)[ind]>=tol)
      {
        // Tolerance target has not been reached yet.
        // Evaluate the derivatives at their current spot
        derivative=ddh_log_likelihood_single_internal(X,ncov,c,x0,nR,pR,i,j);
        second_derivative=d2dh2_log_likelihood_single_internal(X,ncov,c,x0,nR,pR,i,j);
        
//        Rprintf("h: Edge (%d,%d), lb=%f, ub=%f, x0=%f, deriv=%f.\n",i,j, REAL(lb)[ind],REAL(ub)[ind],REAL(x0)[ind],derivative);
        if(!(derivative==derivative)) {
          Rprintf("Optimize h: SEVERE ERROR!!!! SEVERE ERROR!!!!\n");
          return(x0);
        }
        
        // Update boundaries
        if(derivative>0)
          REAL(lb)[ind]=REAL(x0)[ind];
        if(derivative<0)
          REAL(ub)[ind]=REAL(x0)[ind];
        if(derivative==0)
        {
          REAL(lb)[ind]=REAL(x0)[ind];
          REAL(ub)[ind]=REAL(x0)[ind];
        }
            
        // Compute Newton Step if in concave area
        if(second_derivative<0)
          newton=REAL(x0)[ind]-derivative/second_derivative;
        else
          newton=0.5*(REAL(lb)[ind]+REAL(ub)[ind]);
            
        // Check if Newton step falls out of the interval or gives too little improvement, if so do bisection
        if((newton<=REAL(lb)[ind]) | (newton>=REAL(ub)[ind]) | (myabs(newton-REAL(x0)[ind])<=0.001*(REAL(ub)[ind]-REAL(lb)[ind])))
          REAL(x0)[ind]=0.5*(REAL(lb)[ind]+REAL(ub)[ind]);
        else
          REAL(x0)[ind]=newton;
      }
    }
  }
  
//  Rprintf("Hallo8\n");
      
  UNPROTECT(3);
  return(x0);
}