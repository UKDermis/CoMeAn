
/*
   The following code is taken wholesale from the WGCNA package.
   The original paper can be found here:
   Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.”
   BMC Bioinformatics, 559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559.

   The reason why we aren't including the WGCNA package as a dependency is because the package relies on a few libraries
   which aren't available for R >= 4.0.0 as of writing this code.
*/


/*
Copyright (C) 2008 Peter Langfelder; parts based on R by R Development team
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


// #include "corFunctions.h"


#ifndef __corFunctions_h__

#define __corFunctions_h__

#define minSizeForThreading	100

#include <R.h>
#include <Rinternals.h>



void cor1Fast(double * x, int * nrow, int * ncol, double * weights, double * quick,
          int * cosine,
          double * result, int *nNA, int * err,
          int * nThreads,
          int * verbose, int * indent);

void bicor1Fast(double * x, int * nrow, int * ncol,
            double * maxPOutliers, double * quick,
            int * fallback, int * cosine,
            double * result, int *nNA, int * err,
            int * warn,
            int * nThreads,
            int * verbose, int * indent);

void bicorFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           int * robustX, int * robustY,
           double * maxPOutliers, double * quick,
           int * fallback,
           int * cosineX, int * cosineY,
           double * result, int *nNA, int * err,
           int * warnX, int * warnY,
           int * nThreads,
           int * verbose, int * indent);

void corFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           double * weights_x, double * weights_y,
           double * quick,
           int * cosineX, int * cosineY,
           double * result, int *nNA, int * err,
           int * nThreads,
           int * verbose, int * indent);

SEXP cor1Fast_call(SEXP x_s, SEXP weights, SEXP quick_s, SEXP cosine_s,
                   SEXP nNA_s, SEXP err_s,
                   SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP corFast_call(SEXP x_s, SEXP y_s,
                 SEXP weights_x_s, SEXP weights_y_s,
                 SEXP quick_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP bicor1_call(SEXP x_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s, SEXP cosine_s,
                 SEXP nNA_s, SEXP err_s, SEXP warn_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP bicor2_call(SEXP x_s, SEXP y_s,
                 SEXP robustX_s, SEXP robustY_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP warnX_s, SEXP warnY_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);


#endif

// #include "conditionalThreading.h"



#ifndef __conditionalThreading_h__
#define __conditionalThreading_h__

#define MxThreads      128

#ifdef WITH_THREADS

  // #warning Including pthread headers.

  #include <unistd.h>
  #include <pthread.h>

#else

  // define fake pthread functions so we don't have to put a #ifdef everywhere
  //
  // This prevents competing definitions of pthread types to be included
  #define _BITS_PTHREADTYPES_H

  typedef int pthread_mutex_t;
  typedef int pthread_t;
  // in the original code this was called "pthread_attr_t",
  // which causes issues with arm64 macs, as the OS uses this variable internally already, hence the change
  typedef int pthread_attr;

  #define PTHREAD_MUTEX_INITIALIZER 0

  static inline void pthread_mutex_lock ( pthread_mutex_t * lock ) { }
  static inline void pthread_mutex_unlock ( pthread_mutex_t * lock ) { }

  static inline int pthread_join ( pthread_t t, void ** p) { return 0; }

#endif


// Conditional pthread routines

static inline void pthread_mutex_lock_c( pthread_mutex_t * lock, int threaded)
{
  if (threaded) pthread_mutex_lock(lock);
}

static inline void pthread_mutex_unlock_c(pthread_mutex_t * lock, int threaded)
{
  if (threaded) pthread_mutex_unlock(lock);
}

static inline int pthread_create_c(pthread_t *thread, const pthread_attr *attr,
    void *(*start_routine)(void*), void *arg, int threaded)
{
  #ifdef WITH_THREADS
  if (threaded)
    return pthread_create(thread, attr, start_routine, arg);
  else
  #endif
    (*start_routine)(arg);
  return 0;
}

static inline int pthread_join_c(pthread_t thread, void * * value_ptr, int threaded)
{
  if (threaded) return pthread_join(thread, (void * *) value_ptr);
  return 0;
}


#endif

// #include "pivot.h"

#ifndef __pivot_h__

#define __pivot_h__

#include <R.h>
#include <Rinternals.h>

void RprintV(double * v, size_t l);
double vMax(double * v, size_t len);
double vMin(double * v, size_t len);
double pivot(double * v, size_t len, double target);

typedef struct
{
  double val;
  size_t index;
} orderStructure;

int compareOrderStructure(const orderStructure * os1, const orderStructure * os2);
void qorder_internal(double * x, size_t n, orderStructure * os);

SEXP qorder(SEXP data);

#endif

// #include "pivot.c"



#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

#include "pivot.h"

void RprintV(double * v, size_t l)
{
  for (size_t i=0; i<l; i++) Rprintf("%5.3f ", v[i]);
  Rprintf("\n");
}

double vMax(double * v, size_t len)
{
  double mx = v[0];
  for (size_t i=1; i<len; i++)
    if (v[i] > mx) mx = v[i];
  return mx;
}

double vMin(double * v, size_t len)
{
  double mn = v[0];
  for (size_t i=1; i<len; i++)
    if (v[i] < mn) mn = v[i];
  return mn;
}


double pivot(double * v, size_t len, double target)
{
  // Rprintf("Entering pivot with len=%d and target=%f\n   ", len, target);
  // RprintV(v, len);

  if (len > 2)
  {
    // pick the pivot, say as the median of the first, middle and last
    size_t i1 = 0, i2 = len-1, i3 = (len-1)/2, ip;
    if (v[i1] <= v[i2])
    {
      if (v[i2] <= v[i3])
        ip = i2;
      else if (v[i3] >= v[i1])
         ip = i3;
      else
         ip = i1;
    } else {
      if (v[i1] <= v[i3])
        ip = i1;
      else if (v[i2] <= v[i3])
        ip = i3;
      else
        ip = i2;
    }

    // put ip at the end
    double vp = v[ip];
    v[ip] = v[len-1];
    v[len-1] = vp;

    // Rprintf("   pivot value: %5.3f, index: %d\n", vp, ip);

    // pivot everything else
    size_t bound = 0;
    for (size_t i=0; i<len; i++) if (v[i] < vp)
    {
      double x = v[bound];
      v[bound] = v[i];
      v[i] = x;
      bound++;
    }

    v[len-1] = v[bound];
    v[bound] = vp;

    // Rprintf("   After pivoting: bound:%d and vector: ", bound); // RprintV(v, len);

    // Did we find the target?

    double crit = target - bound;
    // Rprintf("   crit: %5.3f\n", crit);
    if (fabs(crit) > 1.0)
    {
      if (crit < 0)
        return pivot(v, bound, target);
      else
        return pivot(v+bound+1, len-bound-1, target-bound-1);
    }
    // Rprintf("vMax(v, bound): %5.3f, vMin(v+bound+1, len-bound-1): %5.3f, vp: %5.3f\n", vMax(v, bound),
                // vMin(v+bound+1, len-bound-1), vp);
    if (crit < 0)
    {
       double v1 = vMax(v, bound);
       return (v1 *(-crit) + vp * (1+crit));
    } // else
    double v2 = vMin(v+bound+1, len-bound-1);
    return (vp * (1-crit) + v2 * crit);
  }
  else if (len==2)
  {
      // Rprintf("  Short v, returning a direct value.\n");
      double v1 = vMin(v, 2);
      double v2 = vMax(v, 2);
      if (target < 0) return v1;
      else if (target > 1) return v2;
      else return (target * v2 + (1-target) * v1);
  }
  else
  {
     // Rprintf("  length 1 v, returning a direct value.\n");
     return v[0];
  }
}


/*====================================================================================================
 *
 * Weighted pivot.
 *
 * arguments:
 *    v is the vector of values;
 *    from, to: indices between which to pivot. v[from],..., v[to-1] will be worked on.
 *    target the quantile which is to be calculated;
 *    w are the weights, assumed to be of length len;
 *    csw is the cumulative sum of weights, csw[i] = sum(w, from=0, to=i)
 *
 *    NOT FINISHED
 *
 *===================================================================================================*/

#define swap(a, b, temp)      { temp = a; a = b; b = temp; }

double pivot_weighted(double * v, size_t from, size_t to, double target,
                      double * w, double * csw)
{
  // Rprintf("Entering pivot with len=%d and target=%f\n   ", len, target);
  // RprintV(v, len);

  size_t len = to-from;
  if (len > 2)
  {
    // pick the pivot, say as the median of the first, middle and last
    size_t i1 = from, i2 = to-1, i3 = (from + to)/2, ip;
    if (v[i1] <= v[i2])
    {
      if (v[i2] <= v[i3])
        ip = i2;
      else if (v[i3] >= v[i1])
         ip = i3;
      else
         ip = i1;
    } else {
      if (v[i1] <= v[i3])
        ip = i1;
      else if (v[i2] <= v[i3])
        ip = i3;
      else
        ip = i2;
    }

    // put v[ip] at the end
    double vp, wp;
    swap(v[ip], v[to-1], vp);
    swap(w[ip], w[to-1], wp);

    // Rprintf("   pivot value: %5.3f, index: %d\n", vp, ip);

    // pivot everything else
    size_t bound = from;
    for (size_t i=from; i<to; i++) if (v[i] < vp)
    {
      double temp;
      swap(v[bound], v[i], temp);
      swap(w[bound], w[i], temp)
      bound++;
    }

    v[len-1] = v[bound]; v[bound] = vp;
    w[len-1] = w[bound]; w[bound] = wp;

    // Update the cumulative sums

    double sum = from > 0 ? csw[from-1] : 0;
    for (size_t i=from; i<to; i++)
    {
       sum = sum + w[i];
       csw[i] = sum;
    }

    // Rprintf("   After pivoting: bound:%d and vector: ", bound); // RprintV(v, len);

    // Did we find the target?

    double crit = target - bound;
    // Rprintf("   crit: %5.3f\n", crit);
    if (fabs(crit) > 1.0)
    {
      if (crit < 0)
        return pivot(v, bound, target);
      else
        return pivot(v+bound+1, len-bound-1, target-bound-1);
    }
    // Rprintf("vMax(v, bound): %5.3f, vMin(v+bound+1, len-bound-1): %5.3f, vp: %5.3f\n", vMax(v, bound),
                // vMin(v+bound+1, len-bound-1), vp);
    if (crit < 0)
    {
       double v1 = vMax(v, bound);
       return (v1 *(-crit) + vp * (1+crit));
    } // else
    double v2 = vMin(v+bound+1, len-bound-1);
    return (vp * (1-crit) + v2 * crit);
  }
  else if (len==2)
  {
      // Rprintf("  Short v, returning a direct value.\n");
      double v1 = vMin(v, 2);
      double v2 = vMax(v, 2);
      if (target < 0) return v1;
      else if (target > 1) return v2;
      else return (target * v2 + (1-target) * v1);
  }
  else
  {
     // Rprintf("  length 1 v, returning a direct value.\n");
     return v[0];
  }
}


/*
 *
 * This isn't needed for now.
 *
 *
void testPivot(double * v, size_t * len, double * target, double * result)
{
   * result = pivot(v, *len, *target);
}
*/

/*****************************************************************************************************
 *
 * Implement order via qsort.
 *
 *****************************************************************************************************/

int compareOrderStructure(const orderStructure * os1, const orderStructure * os2)
{
  if (ISNAN(os1->val)) return 1;
  if (ISNAN(os2->val)) return -1;
  if (os1->val < os2->val) return -1;
  if (os1->val > os2->val) return 1;
  return 0;
}

void qorder_internal(double * x, size_t n, orderStructure * os)
{
  for (R_xlen_t i = 0; i<n; i++)
  {
    (os+i)->val = *(x+i);
    (os+i)->index = i;
  }

  // Rprintf("qorder: calling qsort..");
  qsort(os, (size_t) n, sizeof(orderStructure),
           ((int (*) (const void *, const void *)) compareOrderStructure));
}


SEXP qorder(SEXP data)
{
  R_xlen_t n = Rf_xlength(data);

  // Rprintf("qorder: length of input data is %ld.\n", n);

  double * x = REAL(data);

  orderStructure * os = Calloc((size_t) n, orderStructure);

  qorder_internal(x, (size_t) n, os);

  SEXP ans;
  if (n<(size_t) 0x80000000)
  {
    // Rprintf("..returning integer order.\n");
    PROTECT (ans = allocVector(INTSXP, n));
    int * ansp = INTEGER(ans);
    for (R_xlen_t i = 0; i<n; i++) ansp[i] = (int) ( (os+i)->index+1);
  } else {
    // Rprintf("..returning floating point (double) order.\n");
    PROTECT (ans = allocVector(REALSXP, n));
    double * ansp = REAL(ans);
    for (R_xlen_t i = 0; i<n; i++) ansp[i] = (double) ((os+i)->index+1);
  }
  Free(os);
  UNPROTECT(1);
  return ans;
}

// #include "corFunctions-typeDefs.h"

#ifndef __corFunctions_internal_h__

#define __corFunctions_internal_h__

#define LDOUBLE 	long double

typedef struct
{
   volatile size_t i, n;
}  progressCounter;

/* For each parallel operation will presumably need a separate structure to hold its
 * information, but can define a common structure holding the general information that is needed to
 * calculate correlation. Can keep two versions, one for calculating cor(x), one for cor(x,y).
 * Each specific thread-task specific struct can contain a pointer to the general structure.
 */

// General information for a [bi]cor(x) calculation

typedef struct
{
   double * x, * weights;
   size_t nr, nc;
   double * multMat, * result;
   double * aux;
   size_t *nNAentries;
   int *NAme;
   int zeroMAD;
   int * warn;
   double maxPOutliers;
   double quick;
   int robust, fallback;
   int cosine;
   int id;
   int threaded; 	// This flag will be used to indicate whether the calculation really is threaded.
			// For small problems it doesn't make sense to use threading.
}  cor1ThreadData;

// General information for a [bi]cor(x,y) calculation

typedef struct
{
   cor1ThreadData * x, * y;
}  cor2ThreadData;

// Information for column preparation

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pc;
   pthread_mutex_t * lock;
}  colPrepThreadData;

// Information for symmetrization

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pc;
}  symmThreadData;

// Information for threaded slow calculations for cor1

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pci, * pcj;
   size_t * nSlow, * nNA;
   pthread_mutex_t * lock;
}  slowCalcThreadData;

/*==============================================================================================
 *
 * Threaded 2-variable versions of the correlation functions
 *
 *==============================================================================================
*/


typedef struct
{
   cor2ThreadData * x;
   progressCounter * pci, *pcj;
   size_t * nSlow, * nNA;
   pthread_mutex_t * lock;
   double quick;
}  slowCalc2ThreadData;

// Data for NAing out appropriate rows and columns

typedef struct
{
   cor2ThreadData * x;
   progressCounter * pci, *pcj;
}  NA2ThreadData;

#endif

// #include "corFunctions-utils.h"

#ifndef __corFunctions_utils_h__

#define __corFunctions_utils_h__

enum { noWarning, warnZeroMAD };

double median(double * x, size_t n, int copy, int * err);
double quantile(double * x, size_t n, double q, int copy, int * err);
double quantile_noCopy(double * x, size_t n, double q);
void testMedian(double *x, int * n, double * res);
void testQuantile(double *x, int *n, double *q, double *res);

void prepareColBicor(double * col, size_t nr, double maxPOutliers, int fallback,
                     int cosine,
                     double * res, size_t * nNAentries,
                     int * NAmed, volatile int * zeroMAD,
                     double * aux, double * aux2);
void prepareColCor(double * x, size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean);

void prepareColCor_weighted(double * x, double * weights,
      size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean);

int basic2variableCorrelation(
   double *xx, double *yy,
   size_t nr,
   double *res,
   int cosineX, int cosineY);

int basic2variableCorrelation_weighted(
   double *xx, double *yy,
   double *wx, double *wy,
   size_t nr,
   double *res,
   int cosineX, int cosineY);

void * threadPrepColBicor(void * par);
void * threadPrepColCor(void * par);
void * threadSymmetrize(void * par);
void * threadSlowCalcBicor(void * par);
void * threadSlowCalcCor(void * par);
void * threadNAing(void * par);
void * threadSlowCalcBicor2(void * par);
void * threadSlowCalcCor2(void * par);
void * threadPrepColCor_weighted(void * par);
void * threadSlowCalcCor_weighted(void * par);
void * threadSlowCalcCor2_weighted(void * par);

#endif


#define RefUX	0.5

/*===================================================================================
 *
 * median
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the median of the remaining
// (finite) entries.

double median(double * x, size_t n, int copy, int * err)
{
  double * xx, res;
  if (copy)
  {
    if ( (xx=(double *) malloc(n * sizeof(double)))==NULL )
    {
      Rprintf("Memory allocation error in median(). Could not allocate %d kB.\n",
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;


  *err = 0;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; )
  {
    i--;
    if (ISNAN(xx[i]))
    {
       bound--;
       xx[i] = xx[bound];
       xx[bound] = NA_REAL;
    }
  }

  // Rprintf("Median: n: %d, bound: %d\n", n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot.
    res = pivot(xx, bound, ( 1.0 * (bound-1))/2);

  if (copy) free(xx);

  return res;

}


/*===================================================================================
 *
 * quantile
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the appropriate
// quantile of the remaining (finite) entries.

// q is the quantile: 1/2 will give exactly the median above.

double quantile(double * x, size_t n, double q, int copy, int * err)
{
  double * xx;
  double res;

  if (copy)
  {
    if ( (xx=(double *) malloc(n * sizeof(double)))==NULL )
    {
      Rprintf("Memory allocation error in quantile(). Could not allocate %d kB.\n",
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;


  *err = 0;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; )
  {
    i--;
    if (ISNAN(xx[i]))
    {
       bound--;
       xx[i] = xx[bound];
       xx[bound] = NA_REAL;
    }
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot.
    res = pivot(xx, bound, ( 1.0 * (bound-1))*q);

  if (copy) free(xx);

  return res;

}

double quantile_noCopy(double * x, size_t n, double q)
{
  double res;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; )
  {
    i--;
    if (ISNAN(x[i]))
    {
       bound--;
       x[i] = x[bound];
       x[bound] = NA_REAL;
    }
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot.
    res = pivot(x, bound, ( 1.0 * (bound-1))*q);

  return res;

}

/*==========================================================================================
 *
 * testMedian
 *
 * =========================================================================================*/


void testMedian(double *x, int * n, double * res)
{
  int err;
  *res = median(x, (size_t) *n, 0, &err);
}

/*==========================================================================================
 *
 * testQuantile
 *
 * =========================================================================================*/


void testQuantile(double *x, int *n, double *q, double *res)
{
  int err;
  *res = quantile(x, (size_t) *n, *q, 0, &err);
}


/*==========================================================================================
 *
 * prepareColBicor
 *
 * =========================================================================================*/

// prepareColBicor: calculate median and mad of x and put
// (1-u^2)^2 * (x - median(x))/(9.0 * qnorm75 * mad(x))/ appropriate normalization
// into res.
// res must have enough space allocated to hold the result;
// aux and aux2 each must also have enough space to hold a copy of x.

// maxQoutliers is the maximum allowed proportion of outliers on either side of the median.

// fallback: 1: none, 2: individual, 3: all, 4: force Pearson calculation. 4 is necessary for remedial
// calculations.

// In this case: Pearson pre-calculation entails normalizing columns by mean and variance.

void prepareColBicor(double * col, size_t nr, double maxPOutliers, int fallback,
                     int cosine,
                     double * res, size_t * nNAentries,
                     int * NAmed, volatile int * zeroMAD,
                     double * aux, double * aux2)
{
  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;
  // Note to self: asymptCorr * qnorm75 is very close to 1 and should equal 1 theoretically. Should
  // probably leave them out completely.

  if (fallback==4)
  {
    prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
    return;
  }

  int err = 0;

  // Calculate the median of col

  memcpy((void *)res, (void *)col, nr * sizeof(double));
  double med = median(res, nr, 0, &err);

  // Create a conditional copy of the median
  double medX;
  if (cosine) medX = 0; else medX = med;

  *zeroMAD = 0;
  // calculate absolute deviations from the median

  if (ISNAN(med))
  {
    *NAmed = 1;
    for (size_t k=0; k<nr; k++) *(res + k) = 0;
  } else {
    *NAmed = 0;
    *nNAentries = 0;
    for (size_t k=0; k<nr; k++)
      if (ISNAN(col[k]))
      {
        (*nNAentries)++;
        res[k] = NA_REAL;
        aux[k] = NA_REAL;
      } else {
        res[k] = col[k] - medX;
        aux[k] = fabs(col[k] - med);
      }

    // calculate mad, i.e. median absolute deviation
    double mad = median(aux, nr, 0, &err);

    // If mad is zero, value of fallback decides what is it we will do.
    if (mad==0)
    {
       *zeroMAD = 1;
       switch (fallback)
       {
          case 1:
          {
             // Return after zeoring out results and setting the NAmed flag
             for (size_t k=0; k<nr; k++) res[k] = 0;
             *NAmed = 1;
             return;
          }
          case 2:
             // Switch to Pearson correlation and return
             // Rprintf("mad is zero in a column. Switching to Pearson for this column.\n");
             prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
          case 3:
             // Do nothing: the setting of *zeroMAD above is enough.
             return;
       }
    }

    // We now re-use aux to store a copy of the weights ux. To calculate them, first get (x-med)/(9*mad).

    // Rprintf("median: %6.4f, mad: %6.4f, cosine: %d\n", med, mad, cosine);

    double denom = 9.0 * mad;
    for (size_t k=0; k<nr; k++)
      if (!ISNAN(col[k]))
        aux[k] = (col[k] - med) / denom;
      else
        aux[k] = NA_REAL;

    // Get the low and high quantiles  of ux
    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double lowQ = quantile(aux2, nr, maxPOutliers, 0, &err);

    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double hiQ = quantile(aux2, nr, 1-maxPOutliers, 0, &err);

    // Rprintf("prepareColBicor: lowQ=%f, hiQ = %f\n", lowQ, hiQ);
    // If the low quantile is below -1, rescale the aux (that serve as ux below)
    // such that the low quantile will fall at -1; similarly for the high quantile

    if (lowQ > -RefUX) lowQ = -RefUX;
    if (hiQ < RefUX) hiQ = RefUX;
    lowQ = fabs(lowQ);

    for (size_t k=0; k<nr; k++) if (!ISNAN(aux[k]))
    {
      if (aux[k] < 0)
        aux[k] = aux[k] * RefUX / lowQ;
      else
        aux[k] = aux[k] * RefUX / hiQ;
    }

    // Calculate the (1-ux^2)^2 * (x-median(x))

    LDOUBLE sum = 0;
    for (size_t k=0; k<nr; k++)
      if (!ISNAN(res[k]))
      {
        double ux = aux[k];
        if (fabs(ux) > 1) ux = 1;  // sign of ux doesn't matter.
        ux = 1-ux*ux;
        res[k] *= ux*ux ;
        sum += res[k]*res[k];
      } else
        res[k] = 0;
    sum = sqrtl(sum);
    if (sum==0)
    {
       for (size_t k=0; k<nr; k++)
          res[k] = 0;
       *NAmed = 1;
    } else {
       for (size_t k=0; k<nr; k++)
          res[k] = res[k] / sum;
    }
  }
}

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for Pearson correlation fallback
// and when bicor is called with robustsX or robustY = 0
// if cosine is not zero, the cosine correlation will be calculated.

void prepareColCor(double * x, size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean)
{
  *nNAentries = 0;
  size_t count = 0;
  LDOUBLE mean = 0, sum = 0;
  for (size_t k = 0; k < nr; k++)
    if (!ISNAN(x[k]))
    {
      mean += x[k];
      sum += ((LDOUBLE) x[k])*( (LDOUBLE) x[k]);
      count ++;
    }
  if (count > 0)
  {
    *NAmean = 0;
    *nNAentries = nr-count;
    if (cosine) mean = 0; else mean = mean/count;
    sum = sqrtl(sum - count * mean*mean);
    if (sum > 0)
    {
      // Rprintf("sum: %Le\n", sum);
       for (size_t k=0; k<nr; k++)
         if (!ISNAN(x[k]))
            res[k] = (x[k] - mean)/sum;
         else
            res[k] = 0;
    } else {
       // Rprintf("prepareColCor: have zero variance.\n");
       *NAmean = 1;
       for (size_t k=0; k<nr; k++) res[k] = 0;
    }
  } else {
    *NAmean = 1;
    *nNAentries = nr;
    for (size_t k=0; k<nr; k++)
       res[k] = 0;
  }
}

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for Pearson correlation fallback
// and when bicor is called with robustsX or robustY = 0
// if cosine is not zero, the cosine correlation will be calculated.

void prepareColCor_weighted(double * x, double * weights,
      size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean)
{
  *nNAentries = 0;
  size_t count = 0;
  LDOUBLE mean = 0, wsum = 0, wsumSq = 0, sumSq = 0, sumxwSq = 0;
  for (size_t k = 0; k < nr; k++)
    if (!ISNAN(x[k]) && !ISNAN(weights[k]))
    {
      wsum += weights[k];
      mean += x[k] * weights[k];
      sumSq += ((LDOUBLE) x[k])* x[k] * weights[k] * weights[k];
      sumxwSq += ( (LDOUBLE) x[k]) * weights[k] * weights[k];
      wsumSq += ( (LDOUBLE) weights[k]) * weights[k];
      count ++;
    }
  if (count > 0)
  {
    *NAmean = 0;
    *nNAentries = nr-count;
    if (cosine) mean = 0; else mean = mean/wsum;
    sumSq = sqrtl(sumSq - 2*mean * sumxwSq  + mean*mean * wsumSq);
    //Rprintf("\nprepareColCor_weighted: \n");
    //Rprintf("  mean: %5.3Lf, sumSq: %5.3Lf\n", mean, sumSq);
    //Rprintf("  x: "); RprintV(x, nr);
    //Rprintf("  weights: "); RprintV(weights, nr);
    if ((wsum > 0) && (sumSq > 0))
    {
      // Rprintf("sum: %Le\n", sum);
       for (size_t k=0; k<nr; k++)
         if (!ISNAN(x[k]))
            res[k] = weights[k] * (x[k] - mean)/sumSq;
         else
            res[k] = 0;
    } else {
       // Rprintf("prepareColCor: have zero variance.\n");
       *NAmean = 1;
       for (size_t k=0; k<nr; k++) res[k] = 0;
    }
    //Rprintf("res: "); RprintV(res, nr);
  } else {
    *NAmean = 1;
    *nNAentries = nr;
    for (size_t k=0; k<nr; k++) res[k] = 0;
  }
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for weighted pearson correlation
 *
 *===================================================================================================
*/

// basic function that calculates the (unweighted) correlation of two columns.

// The input pointers must point to the start of the rows in x, y
// The res pointer must point to the appropriate component of the output.
// The return value is 1 if the result is undefined (NA) and 0 if it is valid.

int basic2variableCorrelation(
   double *xx, double *yy,
   size_t nr,
   double *res,
   int cosineX, int cosineY)
{
  LDOUBLE sumxy = 0, sumx = 0, sumy = 0, sumxs = 0, sumys = 0;
  size_t count = 0;
  double vx, vy;
  for (size_t k=0; k<nr; k++)
  {
     vx = *xx; vy = *yy;
     if (!ISNAN(vx) && !ISNAN(vy))
     {
       count ++;
       sumxy += vx * vy;
       sumx += vx;
       sumy += vy;
       sumxs += vx*vx;
       sumys += vy*vy;
     }
     xx++; yy++;
  }
  if (count==0)
  {
      *res = NA_REAL;
      return 1;
  } else {
      if (cosineX) sumx = 0;
      if (cosineY) sumy = 0;
      LDOUBLE varx = sumxs - sumx*sumx/count,
              vary = sumys - sumy*sumy/count;
      if (varx==0 || vary==0)
      {
         *res = NA_REAL;
         return 1;
      } else
         *res = (double) ( (sumxy - sumx * sumy/count)/ sqrtl( varx * vary));
  }
  return 0;
}

// basic function that calculates the weighted correlation of two columns with two different sets of weights.
// The input pointers must point to the start of the rows in x, y, weights for x and weights for y.
//
// The res pointer must point to the appropriate component of the output.
//
// The return value is 1 if the result is undefined (NA) and 0 if it is valid.

int basic2variableCorrelation_weighted(
   double *xx, double *yy,
   double *wx, double *wy,
   size_t nr,
   // output
   double *res,
   // options
   int cosineX, int cosineY)
{
  LDOUBLE
     sum_wx_x = 0, sum_w_x = 0, sum_wSq_x = 0, sum_xwSq_x = 0, sum_Sq_x = 0,
     sum_wx_y = 0, sum_w_y = 0, sum_wSq_y = 0, sum_xwSq_y = 0, sum_Sq_y = 0,
     sum_wwxy = 0, sum_wwy = 0, sum_wwx = 0, sum_ww = 0;

  size_t count = 0;

  for (size_t k=0; k<nr; k++)
  {
     double vx = *xx, vy = *yy, vwx = *wx, vwy = *wy;
     if (!ISNAN(vx) && !ISNAN(vy) && !ISNAN(vwx) && !ISNAN(vwy))
     {
       double ww = vwx * vwy;
       count ++;
       sum_wx_x += vx * vwx;
       sum_Sq_x += ((LDOUBLE) vx) * vx * vwx * vwx;
       sum_xwSq_x += ( (LDOUBLE) vx) * vwx * vwx;
       sum_w_x += vwx;
       sum_wSq_x += ( (LDOUBLE) vwx) * vwx;

       sum_wx_y += vy * vwy;
       sum_Sq_y += ((LDOUBLE) vy) * vy * vwy * vwy;
       sum_xwSq_y += ( (LDOUBLE) vy) * vwy * vwy;
       sum_w_y += vwy;
       sum_wSq_y += ( (LDOUBLE) vwy) * vwy;

       sum_wwxy += ( (LDOUBLE) vx) * vy * ww;
       sum_wwx += ( (LDOUBLE) vx) * ww;
       sum_wwy += ( (LDOUBLE) vy) * ww;
       sum_ww += ww;
     }
     xx++; yy++;
     wx++; wy++;
  }
  // Rprintf("Count of included observations: %d\n", count);
  if (count==0)
  {
    *res = NA_REAL;
    return 1;
  } else {
    double
      mean_x = cosineX ? 0 : sum_wx_x/sum_w_x,
      mean_y = cosineY ? 0 : sum_wx_y/sum_w_y,
      varx = sum_Sq_x - 2*mean_x * sum_xwSq_x + mean_x * mean_x * sum_wSq_x,
      vary = sum_Sq_y - 2*mean_y * sum_xwSq_y + mean_y * mean_y * sum_wSq_y;

    if (varx==0 || vary==0)
    {
      *res = NA_REAL;
      return 1;
    } else
       *res = (double) ( (sum_wwxy - mean_x*sum_wwy - mean_y * sum_wwx + sum_ww *  mean_x * mean_y)/
                     sqrt( varx * vary));

  }
  return 0;
}

/*=========================================================================================================
 *
 *
 * Threaded functions
 *
 *
 *=========================================================================================================*/


/*======================================================================================
 *
 * prepareColBicor
 *
 * =====================================================================================*/

void * threadPrepColBicor(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;

  // Rprintf("Preparing columns: nr = %d, nc = %d\n", x->nr, x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded);
      if (td->pc->i < td->pc->n)
      {
         size_t col = td->pc->i;
         // Rprintf("...working on column %d in thread %d\n", col, td->x->id);
         td->pc->i++;
         pthread_mutex_unlock_c( td->lock, x->threaded );

         prepareColBicor(x->x + col * x->nr,
                         x->nr,
                         x->maxPOutliers,
                         x->fallback,
                         x->cosine,
                         x->multMat + col * x->nr,
                         x->nNAentries + col,
                         x->NAme + col,
                         &(x->zeroMAD),
                         x->aux,
                         x->aux + x->nr);
         // if (x->zeroMAD > 0) { Rprintf("threadPrepColBicor: mad was zero in column %d.\n", col); }
         if (x->zeroMAD > 0) *(x->warn) = warnZeroMAD;
         if ( (x->zeroMAD > 0) && (x->fallback==3))
         {
           pthread_mutex_lock_c( td->lock, x->threaded );
           // Rprintf("threadPrepColBicor: Moving counter from %d %d to end at %d in thread %d.\n",
                   // col, td->pc->i, td->pc->n, x->id);
           x->zeroMAD = col+1; td->pc->i = td->pc->n;
           pthread_mutex_unlock_c( td->lock, x->threaded );
         }
      } else
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
}


/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for the fast calculation of Pearson correlation
// and when bicor is called with robustsX or robustY = 0


void * threadPrepColCor(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;
  //Rprintf("threadPrepColCor: starting in thread %d: counter.i = %d, counter.n = %d, nc = %d.\n",
  //         td->x->id, td->pc->i, td->pc->n, td->x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded );
      int col = td->pc->i;
      if (col < td->x->nc)
      {
         td->pc->i++;
    //     Rprintf("threadPrepColCor: preparing column %d in thread %d.\n", col, td->x->id);
         pthread_mutex_unlock_c( td->lock, x->threaded );

         prepareColCor(x->x + col * x->nr,
                       x->nr,
                       x->cosine,
                       x->multMat + col * x->nr,
                       x->nNAentries + col,
                       x->NAme + col);
      } else
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded symmetrization and NA'ing out of rows and columns with NA means
 *
 *===================================================================================================
*/

void * threadSymmetrize(void * par)
{
  symmThreadData * td = (symmThreadData *) par;
  cor1ThreadData * x = td->x;

  size_t nc = x->nc;
  double * result = x->result;
  int * NAmean = x->NAme;
  size_t col = 0;
  while ( (col = td->pc->i) < nc)
  {
      // Symmetrize the column
      // point counter to the next column
      td->pc->i = col+1;
      // and update the matrix. Start at j=col to check for values greater than 1.
      if (NAmean[col] == 0)
      {
        double * resx = result + col*nc + col;
        // Rprintf("Symmetrizing row %d to the same column.\n", col);
        for (size_t j=col; j<nc; j++)
        {
          // Rprintf("..symmetrizing element %d...\n", j);
          if (NAmean[j] == 0)
          {
            if (!ISNAN(*resx))
            {
               if (*resx > 1.0) *resx = 1.0;
               if (*resx < -1.0) *resx = -1.0;
            }
            result[j*nc + col] = *resx;
          }
          resx ++;
        }
      } else {
        // Rprintf("NA-ing out column and row %d\n", col);
        for (size_t j=0; j<nc; j++)
        {
           result[col*nc + j] = NA_REAL;
           result[j*nc + col] = NA_REAL;
        }
      }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for bicor
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcBicor(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * multMat = td->x->multMat;
  double * result = td->x->result;
  int fbx = td->x->fallback;
  int cosine = td->x->cosine;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  double maxPOutliers = td->x->maxPOutliers;

  double * xx = td->x->aux, * yy = xx + nr;
  double * xxx = xx + 2*nr, * yyy = xx + 3*nr;
  double * xx2 = xx + 4*nr, * yy2 = xx + 5*nr;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  if (fbx==3) fbx = 2; // For these calculations can't go back and redo everything


  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc)
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) &&
              ((NAmean[i] > 0) || (NAmean[j] > 0) ||
               ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );

     if ((i < nc1) && (j < nc) )
     {
        // Rprintf("Recalculating row %d and column %d, column size %d\n", i, j, nr);
        memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
        memcpy((void *)yy, (void *)(x + j*nr), nr * sizeof(double));

        size_t nNAx = 0, nNAy = 0;
        for (size_t k=0; k<nr; k++)
        {
           if (ISNAN(xx[k])) yy[k] = NA_REAL;
           if (ISNAN(yy[k])) xx[k] = NA_REAL;
           if (ISNAN(xx[k])) nNAx++;
           if (ISNAN(yy[k])) nNAy++;
        }
        int NAx = 0, NAy = 0;

        if ((nNAx - nNAentries[i] > maxDiffNA) || (nNAy-nNAentries[j] > maxDiffNA))
        {
            // must recalculate the auxiliary variables for both columns
            size_t temp = 0;
            int zeroMAD = 0;
            if (nNAx - nNAentries[i] > maxDiffNA)
               {
                  prepareColBicor(xx, nr, maxPOutliers, fbx, cosine, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->warn) = warnZeroMAD;
               }
               else
                  memcpy((void *) xxx, (void *) (multMat + i * nr),  nr * sizeof(double));
            if (nNAy-nNAentries[j] > maxDiffNA)
               {
                  prepareColBicor(yy, nr, maxPOutliers, fbx, cosine, yyy, &temp, &NAy, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->warn) = warnZeroMAD;
               }
               else
                  memcpy((void *) yyy, (void *) (multMat + j * nr),  nr * sizeof(double));
            if (NAx + NAy==0)
            {
               LDOUBLE sumxy = 0;
               size_t count = 0;
               for (size_t k=0; k<nr; k++)
               {
                 double vx = *(xxx + k), vy = *(yyy +  k);
                 // Rprintf("i: %d, j: %d, k: %d: vx: %e, vy: %e\n", i,j,k,vx,vy);
                 if (!ISNAN(vx) && !ISNAN(vy))
                 {
                   sumxy += vx * vy;
                   count++;
                 }
               }
               if (count==0)
               {
                  result[i*nc + j] = NA_REAL;
                  (*nNA)++;
               } else {
                  result[i*nc + j] = (double) sumxy;
               }
             } else {
                result[i*nc + j] = NA_REAL;
                (*nNA)++;
             }
             // result[j*nc + i] = result[i*nc + j];
             (*nSlow)++;
        }
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for pearson correlation
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * result = td->x->result;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int cosine = td->x->cosine;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  double *xx, *yy;
  double vx, vy;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  // Rprintf("quick:%f\n", td->x->quick);


  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc)
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) &&
               ((NAmean[i] > 0) || (NAmean[j] > 0) ||
                ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));

     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );

     if ((i < nc1) && (j < nc))
     {
        // Rprintf("Recalculating column %d and row %d, column size %d\n", i, j, nr);
        *nNA += basic2variableCorrelation(
                       x + i * nr, x + j * nr,
                       nr, result + i*nc + j,
                       cosine, cosine);
        (*nSlow)++;
     }
  }
  return NULL;
}


/*===================================================================================================
 *
 * Threaded NA-ing
 *
 *===================================================================================================
*/

void * threadNAing(void * par)
{
  NA2ThreadData * td = (NA2ThreadData *) par;

  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc;
  int * NAmedX = td->x->x->NAme;

  size_t ncy = td->x->y->nc;
  int * NAmedY = td->x->y->NAme;

  progressCounter * pci = td->pci;
  progressCounter * pcj = td->pcj;

  // Go row by row

  size_t row = 0, col = 0;

  while  ((row = pci->i) < ncx)
  {
      pci->i = row + 1;
      if (NAmedX[row])
      {
         // Rprintf("NA-ing out column and row %d\n", col);
         for (size_t j=0; j<ncy; j++)
              result[row + j * ncx] = NA_REAL;
      }
  }

  // ... and column by column

  while ( (col = pcj->i) < ncy)
  {
     pcj->i = col + 1;
     if (NAmedY[col])
     {
        // Rprintf("NA-ing out column and row %d\n", col);
        for (size_t i=0; i<ncx; i++)
             result[i + col * ncx] = NA_REAL;
     } else {
        double *resx = result + col*ncx;
        for (size_t i=0; i<ncx; i++)
        {
           if (!ISNAN(*resx))
           {
             if (*resx > 1.0) *resx = 1.0;
             if (*resx < -1.0) *resx = -1.0;
           }
           resx++;
        }
     }
  }

  return NULL;
}


/*===================================================================================================
 *
 * Threaded "slow" calculations for bicor(x,y)
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcBicor2(void * par)
{
  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int robustX = td->x->x->robust;
  int fbx = td->x->x->fallback;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int robustY = td->x->y->robust;
  int fby = td->x->y->fallback;
  int cosineY = td->x->y->cosine;

  double maxPOutliers = td->x->x->maxPOutliers;

  progressCounter * pci = td->pci, * pcj = td->pcj;

  double * xx = td->x->x->aux;
  double * xxx = xx + nr;
  double * xx2 = xx + 2*nr;

  double * yy = td->x->y->aux;
  double * yyy = yy + nr;
  double * yy2 = yy + 2*nr;

  double * xx3, *yy3;

  int maxDiffNA = (int) (td->x->x->quick * nr);

  if (fbx==3) fbx = 2;
  if (fby==3) fby = 2;

  if (!robustX) fbx = 4;
  if (!robustY) fby = 4;

  // Rprintf("Remedial calculation thread #%d: starting at %d and %d\n", td->x->x->id,
  //            pci->i, pcj->i);
  //

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy)
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) &&
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) ||
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );

     if ((i < ncx) && (j < ncy))
     {
        memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
        memcpy((void *)yy, (void *)(y + j*nr), nr * sizeof(double));

        size_t nNAx = 0, nNAy = 0;
        for (size_t k=0; k<nr; k++)
        {
           if (ISNAN(xx[k])) yy[k] = NA_REAL;
           if (ISNAN(yy[k])) xx[k] = NA_REAL;
           if (ISNAN(xx[k])) nNAx++;
           if (ISNAN(yy[k])) nNAy++;
        }
        int NAx = 0, NAy = 0;

        if ((nNAx - nNAentriesX[i] > maxDiffNA) || (nNAy-nNAentriesY[j] > maxDiffNA))
        {
            // Rprintf("Recalculating row %d and column %d, column size %d in thread %d\n", i, j, nr,
            //         td->x->x->id);
            // must recalculate the auxiliary variables for both columns
            size_t temp = 0;
            int zeroMAD = 0;
            if (nNAx - nNAentriesX[i] > maxDiffNA)
            {
               // Rprintf("...Recalculating row... \n");
               //if (robustX && (fbx!=4))
                  prepareColBicor(xx, nr, maxPOutliers, fbx, cosineX, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->x->warn) = warnZeroMAD;
               //else
               //   prepareColCor(xx, nr, xxx, &temp, &NAx);
               xx3 = xxx;
            } else
               xx3 = multMatX + i * nr;
            if (nNAy-nNAentriesY[j] > maxDiffNA)
            {
               // Rprintf("...Recalculating column... \n");
               //if (robustY && (fby!=4))
                  prepareColBicor(yy, nr, maxPOutliers, fby, cosineY, yyy, &temp, &NAy, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->y->warn) = warnZeroMAD;
               //else
               //   prepareColCor(yy, nr, yyy, &temp, &NAy);
               yy3 = yyy;
            } else
               yy3 = multMatY + j * nr;
            if (NAx + NAy==0)
            {
               // LDOUBLE sumxy = 0;
               double sumxy = 0;
               size_t count = 0;
               for (size_t k=0; k<nr; k++)
               {
                 double vx = *(xx3 + k), vy = *(yy3 +  k);
                 // Rprintf("i: %d, j: %d, k: %d: vx: %e, vy: %e\n", i,j,k,vx,vy);
                 if (!ISNAN(vx) && !ISNAN(vy))
                 {
                   sumxy += vx * vy;
                   count++;
                 }
               }
               if (count==0)
               {
                  result[i + j*ncx] = NA_REAL;
                  (*nNA)++;
               } else {
                  result[i + j*ncx] = (double) sumxy;
                  // Rprintf("Recalculated row %d and column %d, column size %d in thread %d: result = %12.6f\n",
                  //         i, j, nr, td->x->x->id,  result[i + j*ncx]);
               }
            } else {
                result[i + j*ncx] = NA_REAL;
                (*nNA)++;
            }
            (*nSlow)++;
        }
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for pearson correlation of 2 variables.
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor2(void * par)
{
  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
//  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
//  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int cosineY = td->x->y->cosine;

  size_t maxDiffNA = (size_t) (td->x->x->quick * nr);

  progressCounter * pci = td->pci, * pcj = td->pcj;
  double * xx, * yy;

  // Rprintf("Will tolerate %d additional NAs\n",  maxDiffNA);
  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  //

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy)
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) &&
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) ||
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );

     if ((i < ncx) && (j < ncy))
     {
        // Rprintf("Recalculating row %d and column %d, column size %d; cosineX: %d, cosineY: %d\n",
        //         i, j, nr, cosineX, cosineY);
        *nNA += basic2variableCorrelation(
                       x + i * nr, y + j * nr,
                       nr, result + i + j*ncx,
                       cosineX, cosineY);
        (*nSlow)++;
     }
  }
  return NULL;
}

/*======================================================================================
 *
 * threaded prepareColCor_weighted
 *
 * =====================================================================================*/

// Used for the fast calculation of Pearson correlation
// and when bicor is called with robustsX or robustY = 0


void * threadPrepColCor_weighted(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;
  //Rprintf("threadPrepColCor: starting in thread %d: counter.i = %d, counter.n = %d, nc = %d.\n",
  //         td->x->id, td->pc->i, td->pc->n, td->x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded );
      int col = td->pc->i;
      if (col < td->x->nc)
      {
         td->pc->i++;
    //     Rprintf("threadPrepColCor: preparing column %d in thread %d.\n", col, td->x->id);
         pthread_mutex_unlock_c( td->lock, x->threaded );

         prepareColCor_weighted(x->x + col * x->nr,
                       x->weights + col * x->nr,
                       x->nr,
                       x->cosine,
                       x->multMat + col * x->nr,
                       x->nNAentries + col,
                       x->NAme + col);
      } else
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
}

void * threadSlowCalcCor_weighted(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * weights = td->x->weights;
  double * result = td->x->result;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int cosine = td->x->cosine;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  // Rprintf("quick:%f\n", td->x->quick);


  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc)
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) &&
               ((NAmean[i] > 0) || (NAmean[j] > 0) ||
                ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));

     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );

     if ((i < nc1) && (j < nc))
     {
        // Rprintf("Recalculating column %d and row %d, column size %d\n", i, j, nr);
        *nNA += basic2variableCorrelation_weighted(x + i * nr, x + j * nr,
                        weights + i * nr, weights + j * nr,
                        nr, result + i*nc + j,
                        cosine, cosine);
        (*nSlow)++;
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for weighted pearson correlation of 2 variables.
 *
 *===================================================================================================
*/
// The search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor2_weighted(void * par)
{

  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
  double * weights_x = td->x->x->weights;
//  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
  double * weights_y = td->x->y->weights;
//  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int cosineY = td->x->y->cosine;

  size_t maxDiffNA = (size_t) (td->x->x->quick * nr);

  progressCounter * pci = td->pci, * pcj = td->pcj;

  double * xx, * yy;
  double vx = 0, vy = 0;

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy)
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) &&
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) ||
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );

     if ((i < ncx) && (j < ncy))
     {
        // Rprintf("Recalculating row %d and column %d, column size %d; cosineX: %d, cosineY: %d\n",
        //         i, j, nr, cosineX, cosineY);
        *nNA += basic2variableCorrelation_weighted(x + i * nr, y + j * nr,
                        weights_x + i * nr, weights_y + j * nr,
                        nr, result + i + j * ncx,
                        cosineX, cosineY);
        (*nSlow)++;
     }
  }
  return NULL;
}

#include <stdio.h>
//#include <stdlib.h>

#include <sys/time.h>

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/libextern.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>


/*========================================================================
 *
 * Short test code to see whether parallel code can be incorporated into R
 *
 * =======================================================================
 */

int nProcessors()
{
#ifdef WITH_THREADS
#ifdef _SC_NPROCESSORS_CONF
  long nProcessorsOnline = sysconf(_SC_NPROCESSORS_ONLN);
#else
  long nProcessorsOnline = 2;
#endif
#else
  long nProcessorsOnline = 1;
#endif
  return (int) nProcessorsOnline;
}

// Function to calculate suitable number of threads to use.

int useNThreads(size_t n, int nThreadsRequested)
{
#ifdef WITH_THREADS
  int nt = nThreadsRequested;
  if ((nt < 1) || (nt > MxThreads))
  {
    nt = nProcessors();
    if (nt >MxThreads) nt = MxThreads;
  }
  if (n < nt * minSizeForThreading) nt = (n/minSizeForThreading) + 1;
  return nt;
#else
  // Silence "unused argument" warning
  n = n+1;
  return 1;
#endif
}

//===================================================================================================

// Pearson correlation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's
// and uses threading to speed up the rest of the calculation.

//===================================================================================================


// C-level correlation calculation

void cor1Fast(double * x, int * nrow, int * ncol,
          double * weights, double * quick,
          int * cosine,
          double * result, int *nNA, int * err,
          int * nThreads,
          int * verbose, int * indent)
{
  size_t nr = (size_t) *nrow, nc = (size_t) *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *err = 0;

  size_t nNA_ext = 0;

  // Allocate space for various variables

  double * multMat;
  size_t * nNAentries;
  int *NAmean;

  // This matrix will hold preprocessed entries that can be simply multiplied together to get the
  // numerator

  if ( (multMat = (double *) malloc(nc*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor1: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = (size_t *) malloc(nc * sizeof(size_t)))==NULL )
  {
    free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the mean of each column is NA

  if ( (NAmean = (int *) malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( nc*nc, *nThreads);

  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  // double * aux[MxThreads];

  // for (int t=0; t < nt; t++)
  // {
     // if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     // {
       // *err = 1;
       // Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       // for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       // free(NAmean); free(nNAentries); free(multMat);
       // return;
     // }
  // }

  // Put the general data of the correlation calculation into a structure that can be passed on to
  // threads.

  cor1ThreadData thrdInfo[MxThreads];
  for (int t = 0; t < nt; t++)
  {
     thrdInfo[t].x = x;
     thrdInfo[t].weights = weights;
     thrdInfo[t].nr = nr;
     thrdInfo[t].nc = nc;
     thrdInfo[t].multMat = multMat;
     thrdInfo[t].result = result;
     thrdInfo[t].nNAentries = nNAentries;
     thrdInfo[t].NAme = NAmean;
     thrdInfo[t].quick = *quick;
     thrdInfo[t].cosine = *cosine;
     thrdInfo[t].id = t;
     thrdInfo[t].threaded = (nt > 1);
  }

  // Column preparation (calculation of the matrix to be multiplied) in a threaded form.

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
  progressCounter pc;

  pc.i = 0;
  pc.n = nc;

  // Rprintf("Preparing columns...\n");
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfo[t];
    cptd[t].pc = &pc;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t],
                    NULL,
                    weights==NULL ? threadPrepColCor : threadPrepColCor_weighted,
                    (void *) &cptd[t],
                    thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
      *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);

  // Rprintf("done...\n");
  // Rprintf("NAmean:");
  // for (int i=0; i<nc; i++) Rprintf(" %d,", NAmean[i]);
  // Rprintf("\n");

  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  F77_NAME(dsyrk)("L", "T", ncol, nrow, & alpha, multMat, nrow, & beta, result, ncol FCONE FCONE);

  size_t nSlow = 0;

  // Rprintf("nNAentries values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", nNAentries[i]);
  // Rprintf("\n");
  //
  //

  if (*quick < 1.0)
  {
      // Parallelized slow calculations
      slowCalcThreadData  sctd[MxThreads];
      progressCounter pci, pcj;
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pci.i = 0;
      pci.n = nc;
      pcj.i = 1;
      pcj.n = nc;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pci;
        sctd[t].pcj = &pcj;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        status[t] = pthread_create_c(&thr3[t], NULL,
                weights==NULL? threadSlowCalcCor : threadSlowCalcCor_weighted,
                (void *) &sctd[t], thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
         if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfo[t].threaded);

      // Rprintf("done...\n");

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow*2 ) / (nc*(nc-1)) );
  }

  // Symmetrize the result and set all rows and columns with NA means to zero

  symmThreadData  std[MxThreads];
  // reset the progress counter
  pc.i = 0;
  pc.n = nc;

  pthread_t  thr2[MxThreads];

  // Rprintf("symmetrizing... nt=%d\n", nt);
  for (int t=0; t<nt; t++)
  {
    std[t].x = &thrdInfo[t];
    std[t].pc = &pc;
    status[t] = pthread_create_c(&thr2[t], NULL, threadSymmetrize, (void *) &std[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
     if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfo[t].threaded);

  // Rprintf("done... nt=%d\n", nt);
  // Here I need to recalculate results that have NA's in them.

  // for (int t=nt-1; t >= 0; t--) free(aux[t]);

  //Rprintf("End of cor1Fast (1): err = %d\n", *err);
  //*nNA = 1234;
  //Rprintf("End of cor1Fast (2): err = %d\n", *err);
  *nNA = (int) nNA_ext;
  //Rprintf("End of cor1Fast (3): err = %d\n", *err);
  free(NAmean);
  free(nNAentries);
  free(multMat);
}


//===================================================================================================

// bicorrelation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's
// and is threaded to speed up the rest of the calculation.

//===================================================================================================


void bicor1Fast(double * x, int * nrow, int * ncol, double * maxPOutliers,
            double * quick, int * fallback, int * cosine,
            double * result, int *nNA, int * err,
            int * warn,
            int * nThreads,
            int * verbose, int * indent)
{
  size_t nr = *nrow, nc = *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *nNA = 0;
  *warn = noWarning;
  *err = 0;

  size_t nNA_ext = 0;

  // Allocate space for various variables

  double * multMat;
  size_t * nNAentries;
  int  *NAmed;

  if ( (multMat = (double *) malloc(nc*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor1: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = (size_t *) malloc(nc * sizeof(size_t)))==NULL )
  {
    free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the mean of each column is NA

  if ( (NAmed = (int *) malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( nc*nc, *nThreads);

  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  double * aux[MxThreads];

  for (int t=0; t < nt; t++)
  {
     if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     {
       *err = 1;
       Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       free(NAmed); free(nNAentries); free(multMat);
       return;
     }
  }

  // Put the general data of the correlation calculation into a structure that can be passed on to
  // threads.

  cor1ThreadData thrdInfo[MxThreads];
  for (int t = 0; t < nt; t++)
  {
     thrdInfo[t].x = x;
     thrdInfo[t].weights = NULL;
     thrdInfo[t].nr = nr;
     thrdInfo[t].nc = nc;
     thrdInfo[t].multMat = multMat;
     thrdInfo[t].result = result;
     thrdInfo[t].nNAentries = nNAentries;
     thrdInfo[t].NAme = NAmed;
     thrdInfo[t].zeroMAD = 0;
     thrdInfo[t].warn = warn;   // point the pointer
     thrdInfo[t].aux = aux[t];
     thrdInfo[t].robust = 0;
     thrdInfo[t].fallback = *fallback;
     thrdInfo[t].quick = *quick;
     thrdInfo[t].cosine = *cosine;
     thrdInfo[t].maxPOutliers = *maxPOutliers;
     thrdInfo[t].id = t;
     thrdInfo[t].threaded = (nt > 1);
  }

  // Column preparation (calculation of the matrix to be multiplied) in a threaded form.

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
  progressCounter pc;

  pc.i = 0;
  pc.n = nc;

  // Rprintf("Preparing columns...\n");
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfo[t];
    cptd[t].pc = &pc;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "WARNING: RETURNED RESULTS WILL BE INCORRECT.");
          *err = 2;
    }

  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);

  int pearson = 0;

  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfo[t].zeroMAD > 0)
    {
      pearson = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x): Thread %d (of %d) reported zero MAD in column %d. %s",
                t, nt, thrdInfo[t].zeroMAD, "Switching to Pearson correlation.\n");
    }
    if (pearson==1) // Re-do all column preparations using Pearson preparation.
    {
      // Set fallback to 4 for slow calculations below.
      for (int t = 0; t < nt; t++) thrdInfo[t].fallback = 4;

      pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
      pc.i = 0;
      pc.n = nc;

      for (int t=0; t<nt; t++)
      {
        cptd[t].lock = &mutex2;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t], thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++)
         if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);
    }
  }

  // Rprintf("done...\n");
  // Rprintf("NAmed:");
  // for (int i=0; i<nc; i++) Rprintf(" %d,", NAmed[i]);
  // Rprintf("\n");

  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  // Rprintf("alpha: %f\n", alpha);
  F77_NAME(dsyrk)("L", "T", ncol, nrow, & alpha, multMat, nrow, & beta, result, ncol FCONE FCONE);

  // Here I need to recalculate results that have NA's in them.

  size_t nSlow = 0;

  // Rprintf("nNAentries values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", nNAentries[i]);
  // Rprintf("\n");
  //
  //

  if (*quick < 1)
  {
      // Parallelized slow calculations
      slowCalcThreadData  sctd[MxThreads];
      progressCounter pci, pcj;
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pci.i = 0;
      pci.n = nc;
      pcj.i = 1;
      pcj.n = nc;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pci;
        sctd[t].pcj = &pcj;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        status[t] = pthread_create_c(&thr3[t], NULL, threadSlowCalcBicor, (void *) &sctd[t],
                    thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
         if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfo[t].threaded);

      // Rprintf("done...\n");

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow*2 ) / (nc*(nc-1)) );
  }

  // Symmetrize the result and set all rows and columns with NA means to zero

  symmThreadData  std[MxThreads];
  // reset the progress counter
  pc.i = 0;
  pc.n = nc;

  pthread_t  thr2[MxThreads];

  // Rprintf("symmetrizing... nt=%d\n", nt);
  for (int t=0; t<nt; t++)
  {
    std[t].x = &thrdInfo[t];
    std[t].pc = &pc;
    status[t] = pthread_create_c(&thr2[t], NULL, threadSymmetrize, (void *) &std[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfo[t].threaded);

  for (int t=nt-1; t >= 0; t--) free(aux[t]);

  *nNA = (int) nNA_ext;

  free(NAmed);
  free(nNAentries);
  free(multMat);
}

//===================================================================================================
//
// Two-variable bicorrelation. Basically the same as bicor1, just must calculate the whole matrix.
// If robustX,Y is zero, the corresponding variable will be treated as in pearson correlation.
//
//===================================================================================================

void bicorFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           int * robustX, int * robustY, double *maxPOutliers,
           double * quick, int * fallback,
           int * cosineX, int * cosineY,
           double * result, int *nNA, int * err,
           int * warnX, int * warnY,
           int * nThreads,
           int * verbose, int * indent)
{
  size_t nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *warnX = noWarning;
  *warnY = noWarning;
  *err = 0;

  size_t nNA_ext = 0;

  double * multMatX, * multMatY;
  size_t * nNAentriesX, * nNAentriesY;
  int *NAmedX, *NAmedY;


  // Rprintf("nr: %d, ncx: %d, ncy: %d\n", nr, ncx, ncy);
  // Rprintf("robustX: %d, robustY: %d, cosineX: %d, cosineY: %d\n", *robustX, *robustY, *cosineX, *cosineY);
  // Rprintf("quick: %12.6f, maxPOutliers: %12.6f\n", *quick, *maxPOutliers);

  // Rprintf("Last few entries of x:\n");
  // for (int i = nr-2; i<nr; i++)
  // {
  //   for (int j = ncx-3; j<ncx; j++)
  //     Rprintf("%12.6f ", x[j*nr + i]);
  //   Rprintf("\n");
  // }

  // Rprintf("Last few entries of y:\n");
  // for (int i = nr-2; i<nr; i++)
  // {
  //   for (int j = ncy-3; j<ncy; j++)
  //     Rprintf("%12.6f ", y[j*nr + i]);
  //   Rprintf("\n");
  // }

  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("bicor: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = (double *) malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = (size_t *) malloc(ncx * sizeof(size_t)))==NULL )
  {
    free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = (size_t *) malloc(ncy * sizeof(size_t)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmedX = (int *) malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmedY = (int *) malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmedX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( ncx* ncy, *nThreads);

  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  double * aux[MxThreads];

  for (int t=0; t < nt; t++)
  {
     if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     {
       *err = 1;
       Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       free(NAmedY); free(NAmedX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
       return;
     }
  }

  cor1ThreadData thrdInfoX[MxThreads];
  cor1ThreadData thrdInfoY[MxThreads];
  cor2ThreadData thrdInfo[MxThreads];

  for (int t = 0; t < nt; t++)
  {
     thrdInfoX[t].x = x;
     thrdInfoX[t].weights = NULL;
     thrdInfoX[t].nr = nr;
     thrdInfoX[t].nc = ncx;
     thrdInfoX[t].multMat = multMatX;
     thrdInfoX[t].result = result;
     thrdInfoX[t].nNAentries = nNAentriesX;
     thrdInfoX[t].NAme = NAmedX;
     thrdInfoX[t].zeroMAD = 0;
     thrdInfoX[t].aux = aux[t];
     thrdInfoX[t].robust = *robustX;
     thrdInfoX[t].fallback = *fallback;
     thrdInfoX[t].maxPOutliers = *maxPOutliers;
     thrdInfoX[t].quick = *quick;
     thrdInfoX[t].cosine = *cosineX;
     thrdInfoX[t].warn = warnX;
     thrdInfoX[t].id = t;
     thrdInfoX[t].threaded = (nt > 1);


     thrdInfoY[t].x = y;
     thrdInfoY[t].weights = NULL;
     thrdInfoY[t].nr = nr;
     thrdInfoY[t].nc = ncy;
     thrdInfoY[t].multMat = multMatY;
     thrdInfoY[t].result = result;
     thrdInfoY[t].nNAentries = nNAentriesY;
     thrdInfoY[t].NAme = NAmedY;
     thrdInfoY[t].zeroMAD = 0;
     thrdInfoY[t].aux = aux[t] + 3 * nr;
     thrdInfoY[t].robust = *robustY;
     thrdInfoY[t].fallback = *fallback;
     thrdInfoY[t].maxPOutliers = *maxPOutliers;
     thrdInfoY[t].quick = *quick;
     thrdInfoY[t].cosine = *cosineY;
     thrdInfoY[t].warn = warnY;
     thrdInfoY[t].id = t;
     thrdInfoY[t].threaded = (nt > 1);

     thrdInfo[t].x = thrdInfoX + t;
     thrdInfo[t].y = thrdInfoY + t;
  }

  // Prepare the multMat columns in X and Y

  // Rprintf(" ..preparing columns in x\n");

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  progressCounter pcX, pcY;
  int pearsonX = 0, pearsonY = 0;

  // Prepare columns in X

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

  pcX.i = 0;
  pcX.n = ncx;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoX[t];
    cptd[t].pc = &pcX;
    cptd[t].lock = &mutex1;
    if (* robustX)
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t],
                                      thrdInfoX[t].threaded);
       else
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t],
                                      thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }
  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // If the fallback method is to re-do everything in Pearson, check whether any columns had zero MAD.
  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfoX[t].zeroMAD > 0)
    {
      pearsonX = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x, y): thread %d of %d reported zero MAD in column %d of x. %s",
                t, nt, thrdInfoX[t].zeroMAD, "Switching to Pearson calculation for x.\n");
    }
    if (pearsonX==1) // Re-do all column preparations
    {
      for (int t = 0; t < nt; t++) thrdInfoX[t].fallback = 4;

      pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
      pcX.i = 0;
      pcX.n = ncx;

      for (int t=0; t<nt; t++)
      {
        cptd[t].lock = &mutex2;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t],
                                     thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++) if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);
    }
  }


  // Prepare columns in Y

  // Rprintf(" ..preparing columns in y\n");
  pthread_mutex_t mutex1Y = PTHREAD_MUTEX_INITIALIZER;

  pcY.i = 0;
  pcY.n = ncy;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoY[t];
    cptd[t].pc = &pcY;
    cptd[t].lock = &mutex1Y;
    if (* robustY)
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t],
                                      thrdInfoX[t].threaded);
       else
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t],
                                     thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0)  pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // If the fallback method is to re-do everything in Pearson, check whether any columns had zero MAD.
  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfoY[t].zeroMAD > 0)
    {
      pearsonY = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x, y): thread %d of %d reported zero MAD in column %d of y. %s",
                t, nt, thrdInfoY[t].zeroMAD, "Switching to Pearson calculation for y.\n");
    }
    if (pearsonY==1) // Re-do all column preparations
    {
      for (int t = 0; t < nt; t++) thrdInfoY[t].fallback = 4;

      pthread_mutex_t mutex2Y = PTHREAD_MUTEX_INITIALIZER;
      pcY.i = 0;
      pcY.n = ncy;

      for (int t=0; t<nt; t++)
      {
      //  Rprintf("Starting pearson re-calculation in thread %d of %d.\n", t, nt);
        cptd[t].lock = &mutex2Y;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t], thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++) if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);
    }
  }

   // Rprintf("nNAentriesX:");
   // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesX[i]);
   // Rprintf("\n");
   // Rprintf("nNAentriesY:");
   // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesY[i]);
   // Rprintf("\n");


  // The main calculation: matrix multiplication

  double alpha = 1.0, beta = 0.0;
  F77_NAME(dgemm)("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx FCONE FCONE);

  // Rprintf("matrix multiplication result:\n");
  // for (int i=0; i<ncx; i++)
  // {
  //   for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", result[i + ncx*j]);
  //   Rprintf("\n");
 //  }
/*
  Rprintf("Last few entries of result just after multiplication:\n");
  for (int i = 0; i<ncx; i++)
  {
    for (int j = 0; j<ncy; j++)
      Rprintf("%12.6f ", result[j*ncx + i]);
    Rprintf("\n");
  }
  Rprintf("multMatX:\n");
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<ncx; j++) Rprintf(" %12.6f ", multMatX[i + nr*j]);
    Rprintf("\n");
  }

  Rprintf("multMatY:\n");
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", multMatY[i + nr*j]);
    Rprintf("\n");
  }
*/

  // Remedial calculations

  size_t nSlow = 0;
  if (*quick < 1.0)
  {
      slowCalc2ThreadData  sctd[MxThreads];
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pcX.i = 0; pcY.i = 0;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pcX;
        sctd[t].pcj = &pcY;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        sctd[t].quick = *quick;
        status[t] = pthread_create_c(&thr3[t], NULL, threadSlowCalcBicor2, (void *) &sctd[t],
                                     thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
        if (status[t]==0)  pthread_join_c(thr3[t], NULL, thrdInfoX[t].threaded);

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow) / (ncx*ncy) );
  }

  // NA out all rows and columns that need it and check for values outside of [-1, 1]

  NA2ThreadData  natd[MxThreads];
  // reset the progress counter
  pcX.i = 0;
  pcY.i = 0;

  pthread_t  thr2[MxThreads];

  for (int t=0; t<nt; t++)
  {
    natd[t].x = &thrdInfo[t];
    natd[t].pci = &pcX;
    natd[t].pcj = &pcY;
    status[t] = pthread_create_c(&thr2[t], NULL, threadNAing, (void *) &natd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
     if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfoX[t].threaded);

  *nNA = (int) nNA_ext;


  // Rprintf("Last few entries of result:\n");
  // for (int i = ncx-4; i<ncx; i++)
  // {
  //  for (int j = ncy-4; j<ncy; j++)
  //    Rprintf("%12.6f ", result[j*ncx + i]);
  //  Rprintf("\n");
  //}

  // Clean up

  for (int t=nt-1; t >= 0; t--) free(aux[t]);
  free(NAmedY);
  free(NAmedX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
}

/*======================================================================================================
 *
 * corFast: fast correlation of 2 matrices
 *
 *======================================================================================================
One important note: if weights_x is not NULL, weights_y is also assumed to be valid.
*/

void corFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           double * weights_x, double * weights_y,
           double * quick,
           int * cosineX, int * cosineY,
           double * result, int *nNA, int * err,
           int * nThreads,
           int * verbose, int * indent)
{
  size_t nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  size_t nNA_ext = 0;
  *err = 0;

  double * multMatX, * multMatY;
  size_t * nNAentriesX, * nNAentriesY;
  int *NAmeanX, *NAmeanY;

  if ( (weights_x == NULL) != (weights_y == NULL))
  {
    *err = 2;
    error("corFast: weights_x and weights_y must both be either NULL or non-NULL.\n");
  }

  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = (double *) malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = (size_t *) malloc(ncx * sizeof(size_t)))==NULL )
  {
    free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = (size_t *) malloc(ncy * sizeof(size_t)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanX = (int *) malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanY = (int *) malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmeanX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( ncx* ncy, *nThreads);

  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  cor1ThreadData thrdInfoX[MxThreads];
  cor1ThreadData thrdInfoY[MxThreads];
  cor2ThreadData thrdInfo[MxThreads];

  for (int t = 0; t < nt; t++)
  {
     thrdInfoX[t].x = x;
     thrdInfoX[t].weights = weights_x;
     thrdInfoX[t].nr = nr;
     thrdInfoX[t].nc = ncx;
     thrdInfoX[t].multMat = multMatX;
     thrdInfoX[t].result = result;
     thrdInfoX[t].nNAentries = nNAentriesX;
     thrdInfoX[t].NAme = NAmeanX;
     thrdInfoX[t].quick = *quick;
     thrdInfoX[t].cosine = *cosineX;
     thrdInfoX[t].maxPOutliers = 1;
     thrdInfoX[t].id = t;
     thrdInfoX[t].threaded = (nt > 1);

     thrdInfoY[t].x = y;
     thrdInfoY[t].weights = weights_y;
     thrdInfoY[t].nr = nr;
     thrdInfoY[t].nc = ncy;
     thrdInfoY[t].multMat = multMatY;
     thrdInfoY[t].result = result;
     thrdInfoY[t].nNAentries = nNAentriesY;
     thrdInfoY[t].NAme = NAmeanY;
     thrdInfoY[t].quick = *quick;
     thrdInfoY[t].cosine = *cosineY;
     thrdInfoY[t].maxPOutliers = 1;
     thrdInfoY[t].id = t;
     thrdInfoY[t].threaded = (nt > 1);

     thrdInfo[t].x = thrdInfoX + t;
     thrdInfo[t].y = thrdInfoY + t;
  }

  // Prepare the multMat columns in X and Y

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  progressCounter pcX, pcY;

  // Prepare columns in X

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

  pcX.i = 0;
  pcX.n = ncx;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoX[t];
    cptd[t].pc = &pcX;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t], NULL,
                   weights_x==NULL ? threadPrepColCor : threadPrepColCor_weighted,
                   (void *) &cptd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }
  for (int t=0; t<nt; t++)
    if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // Prepare columns in Y

  pthread_mutex_t mutex1Y = PTHREAD_MUTEX_INITIALIZER;

  pcY.i = 0;
  pcY.n = ncy;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoY[t];
    cptd[t].pc = &pcY;
    cptd[t].lock = &mutex1Y;
    status[t] = pthread_create_c(&thr[t], NULL,
           weights_y==NULL ? threadPrepColCor: threadPrepColCor_weighted,
           (void *) &cptd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  //Rprintf("multMatX:\n");
  //for (int i=0; i<nr; i++)
  //{
    //for (int j=0; j<ncx; j++) Rprintf(" %12.6f ", multMatX[i + nr*j]);
    //Rprintf("\n");
  //}

  //Rprintf("multMatY:\n");
  //for (int i=0; i<nr; i++)
  //{
    //for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", multMatY[i + nr*j]);
    //Rprintf("\n");
  //}

  // Rprintf("nNAentriesX:");
  // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesX[i]);
  // Rprintf("\n");
  // Rprintf("nNAentriesY:");
  // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesY[i]);
  // Rprintf("\n");


  // The main calculation: matrix multiplication

  double alpha = 1.0, beta = 0.0;
  F77_NAME(dgemm)("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx FCONE FCONE);

  // Remedial calculations

  size_t nSlow = 0;
  if (*quick < 1.0)
  {
      slowCalc2ThreadData  sctd[MxThreads];
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pcX.i = 0; pcY.i = 0;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pcX;
        sctd[t].pcj = &pcY;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        sctd[t].quick = *quick;
        status[t] = pthread_create_c(&thr3[t], NULL,
                 weights_x==NULL ? threadSlowCalcCor2 : threadSlowCalcCor2_weighted,
                 (void *) &sctd[t], thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
          if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfoX[t].threaded);

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow) / (ncx*ncy) );
  }

  // NA out all rows and columns that need it and check for values outside of [-1, 1]
  //
  NA2ThreadData  natd[MxThreads];
  // reset the progress counters
  pcX.i = 0;
  pcY.i = 0;

  pthread_t  thr2[MxThreads];

  for (int t=0; t<nt; t++)
  {
    natd[t].x = &thrdInfo[t];
    natd[t].pci = &pcX;
    natd[t].pcj = &pcY;
    status[t] = pthread_create_c(&thr2[t], NULL, threadNAing, (void *) &natd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0)  pthread_join_c(thr2[t], NULL, thrdInfoX[t].threaded);

  *nNA = (int) nNA_ext;

  // clean up and return

  free(NAmeanY);
  free(NAmeanX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
}

//===================================================================================================
//
// Two-variable Pearson correlation.
//
//===================================================================================================

SEXP bicor2_call(SEXP x_s, SEXP y_s,
                 SEXP robustX_s, SEXP robustY_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP warnX_s, SEXP warnY_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dimX, dimY, cor_s;

  int nr, ncx, ncy;
  int *cosineX, *cosineY;
  int *err, *nThreads, *verbose, *indent, *fallback;
  int *warnX, *warnY, *robustX, *robustY;
  int *nNA;

  double *x, *y, *corMat, *quick, *maxPOutliers;

  /* Get dimensions of 'x'. */
  PROTECT(dimX = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dimX)[0];
  ncx = INTEGER(dimX)[1];
  // Rprintf("Matrix x dimensions: %d %d\n", nr, ncx);
  /* Get dimensions of 'y'. */
  PROTECT(dimY = getAttrib(y_s, R_DimSymbol));
  ncy = INTEGER(dimY)[1];
  // Rprintf("Matrix y dimensions: %d %d\n", INTEGER(dimY)[0], ncy);

  x = REAL(x_s);
  y = REAL(y_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  maxPOutliers = REAL(maxPOutliers_s);
  cosineX = INTEGER(cosineX_s);
  cosineY = INTEGER(cosineY_s);
  robustX = INTEGER(robustX_s);
  robustY = INTEGER(robustY_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);
  fallback = INTEGER(fallback_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, ncx, ncy));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);
  warnX = INTEGER(warnX_s);
  warnY = INTEGER(warnY_s);

  bicorFast(x, &nr, &ncx, y, &ncy,
           robustX, robustY,
           maxPOutliers, quick, fallback,
           cosineX, cosineY,
           corMat, nNA, err,
           warnX, warnY,
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(3);
  return cor_s;
}

SEXP corFast_call(SEXP x_s, SEXP y_s,
                 SEXP weights_x_s, SEXP weights_y_s,
                 SEXP quick_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dimX, dimY, cor_s;

  int nr, ncx, ncy;
  int *cosineX, *cosineY;
  int *err, *nThreads, *verbose, *indent;
  int *nNA;

  double *x, *y, *weights_x, *weights_y, *corMat, *quick;

  /* Get dimensions of 'x'. */
  PROTECT(dimX = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dimX)[0];
  ncx = INTEGER(dimX)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);
  /* Get dimensions of 'y'. */
  PROTECT(dimY = getAttrib(y_s, R_DimSymbol));
  ncy = INTEGER(dimY)[1];

  x = REAL(x_s);
  y = REAL(y_s);

  if (isNull(weights_x_s)) weights_x = NULL; else weights_x = REAL(weights_x_s);
  if (isNull(weights_y_s)) weights_y = NULL; else weights_y = REAL(weights_y_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  cosineX = INTEGER(cosineX_s);
  cosineY = INTEGER(cosineY_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, ncx, ncy));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);

  // Rprintf("Calling cor1Fast...\n");
  corFast(x, &nr, &ncx, y, &ncy,
          weights_x, weights_y,
          quick,
          cosineX, cosineY,
          corMat, nNA, err,
          nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(3);
  return cor_s;
}

// Re-write cor1Fast as a function that can be called using .Call
// Since I don't know how to create and fill lists in C code, I will for now return the nNA and err results
// via supplied arguments. Not ideal but will do.

SEXP cor1Fast_call(SEXP x_s, SEXP weights_s, SEXP quick_s, SEXP cosine_s,
                   SEXP nNA_s, SEXP err_s,
                   SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dim, cor_s;
  // SEXP out, nNA_s, err_s;

  int nr, nc;
  int *cosine, *err, *nThreads, *verbose, *indent;
  int *nNA;

  double *x, *weights, *corMat, *quick;

  /* Get dimensions of 'x'. */
  PROTECT(dim = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dim)[0];
  nc = INTEGER(dim)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);

  x = REAL(x_s);

  if (isNull(weights_s)) weights = NULL; else weights = REAL(weights_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  cosine = INTEGER(cosine_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, nc, nc));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);

  // Rprintf("Difference of nNA and err pointers: %d\n", (int) (err - nNA));

  // Rprintf("Calling cor1Fast...\n");
  cor1Fast(x, &nr, &nc, weights, quick, cosine,
           corMat, nNA, err,
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(2);
  return cor_s;
}

SEXP bicor1_call(SEXP x_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s, SEXP cosine_s,
                 SEXP nNA_s, SEXP err_s, SEXP warn_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dim, cor_s;
  // SEXP out, nNA_s, err_s;

  int nr, nc;
  int *cosine, *err, *nThreads, *verbose, *indent, *fallback;
  int *nNA, *warn;

  double *x, *corMat, *quick, *maxPOutliers;

  /* Get dimensions of 'x'. */
  PROTECT(dim = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dim)[0];
  nc = INTEGER(dim)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);

  x = REAL(x_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  maxPOutliers = REAL(maxPOutliers_s);
  cosine = INTEGER(cosine_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);
  fallback = INTEGER(fallback_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, nc, nc));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);
  warn = INTEGER(warn_s);

  // Rprintf("Calling cor1Fast...\n");
  bicor1Fast(x, &nr, &nc,
           maxPOutliers, quick,
           fallback, cosine,
           corMat, nNA, err,
           warn,
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(2);
  return cor_s;
}