#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "NumericsSparseMatrix.h"
#include "numerics_verbose.h"
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"

#ifdef WITH_FCLIB
static int fccounter =-1;
#endif
static double * alloc_memory_double(unsigned int size, double *p)
{
  double * r = (double *) malloc (size * sizeof(double));
  memcpy(r, p, size * sizeof(double));
  return r;
}

static csi * alloc_memory_csi(unsigned int size, unsigned int *p)
{
  csi * r = (csi *) malloc (size * sizeof(csi));
  for(unsigned int i=0; i<size; ++i)
  {
    r[i] = (csi) p[i];
  }
  return r;
}

int globalFrictionContact_fclib_write(
  GlobalFrictionContactProblem* problem,
  char * title,
  char * description,
  char * mathInfo,
  const char *path);

int gfc3d_LmgcDriver(double *reaction,
                     double *velocity,
                     double *globalVelocity,
                     double *q,
                     double *b,
                     double *mu,
                     double *Mdata,
                     unsigned int nzM,
                     long long int *rowM,
                     long long int *colM,
                     double* Hdata,
                     unsigned int nzH,
                     long long int *rowH,
                     long long int *colH,
                     unsigned int n,
                     unsigned int nc,
                     int solver_id,
                     int isize,
                     int *iparam,
                     int dsize,
                     double *dparam,
                     int verbose_in,
                     int outputFile,
                     int freq_output)
{


  verbose = verbose_in;
 
  /* NumericsMatrix M, H; */
  NumericsMatrix * M =NM_new();
  M->storageType = 2; /* sparse */
  M->size0 = n;
  M->size1 = n;


  NumericsSparseMatrix * SM =newNumericsSparseMatrix();
  M->matrix2 = SM;
  SM->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _M = SM->triplet;
  SM->origin = NS_TRIPLET;

  /* csi * _colM = alloc_memory_csi(nzM, colM); */
  /* csi * _rowM = alloc_memory_csi(nzM, rowM); */
  
  csi * _colM = (csi *) colM;
  csi * _rowM = (csi *) rowM;
  
  
  _M->nzmax = nzM;
  _M->nz = nzM;
  _M->m = M->size0;
  _M->n = M->size1;
  _M->p = (csi *) _colM;
  _M->i = (csi *) _rowM;

  
  
  /* double * _Mdata = alloc_memory_double(nzM, Mdata); */
  double * _Mdata = Mdata;
  _M->x = _Mdata;

  DEBUG_PRINTF("_M->n=%lli\t",_M->n);
  DEBUG_PRINTF("_M->m=%lli\n",_M->m);





  NumericsMatrix * H =NM_new();
  H->storageType = 2;
  H->size0 = M->size0;
  H->size1 = 3 * nc;

  
  NumericsSparseMatrix * SH =newNumericsSparseMatrix();
  H->matrix2 = SH;
  SH->triplet =   (CSparseMatrix * )malloc(sizeof(CSparseMatrix));
  CSparseMatrix * _H = SH->triplet;
  SH->origin = NS_TRIPLET;

  
  DEBUG_EXPR(
    for (unsigned int i =0 ; i <nzH ; i++ )
    {
      printf("%i , rowH = %i, colH = %i, Hdata = %e\n ", i, rowH[i], colH[i], Hdata[i]);
    }
    );
  
  
  /* csi * _colH = alloc_memory_csi(nzH, colH); */
  /* csi * _rowH = alloc_memory_csi(nzH, rowH); */
  /* csi * _colH = (csi *) colH; */
  /* csi * _rowH = (csi *) rowH; */

  _H->nzmax = nzH;
  _H->nz = nzH;
  _H->m = H->size0;
  _H->n = H->size1;
  DEBUG_PRINTF("nzH=%lli\n",nzH);

  _H->p = (csi *) colH;
  _H->i = (csi *) rowH;
  /* double * _Hdata = alloc_memory_double(nzH, Hdata); */
  double * _Hdata = Hdata;
  _H->x = _Hdata;
  DEBUG_EXPR(NM_display(H););
  DEBUG_EXPR(
    for (unsigned int i =0; i <nzH ; i++ )
    {
      printf("%i , _H = %i, _colH = %i, _Hdata = %e\n ", i,  _H->i[i], _H->p[i], _H->x[i]);
    }
    );

  for (int i=0; i< _M->nz; ++i)
  {
    /* _M->p[i] --; */
    /* _M->i[i] --; */
    DEBUG_PRINTF("%d -> %d,%d val : %e\n", i, _M->p[i], _M->i[i], _M->x[i]);

  }

  for (int i=0; i< _H->nz; ++i)
  {
    /* _H->p[i] --; */
    /* _H->i[i] --; */
    DEBUG_PRINTF("%d -> %d,%d val : %e \n", i, _H->p[i], _H->i[i], _H->x[i]);
  }

  GlobalFrictionContactProblem * problem =(GlobalFrictionContactProblem*)malloc(sizeof(GlobalFrictionContactProblem));

  problem->dimension = 3;
  problem->numberOfContacts = nc;
  problem->env = NULL;
  problem->workspace = NULL;

  problem->M = M;
  problem->H = H;
  problem->q = q;
  problem->b = b;
  problem->mu = mu;

  SolverOptions numerics_solver_options;
  
  int infi = gfc3d_setDefaultSolverOptions(&numerics_solver_options, solver_id);
  assert(!infi);
  int iSize_min = isize < numerics_solver_options.iSize ? isize : numerics_solver_options.iSize;
  DEBUG_PRINTF("iSize_min = %i", iSize_min);
  for (int i = 0; i < iSize_min; ++i) 
    numerics_solver_options.iparam[i] = iparam[i];

  int dSize_min = dsize <  numerics_solver_options.dSize ? dsize : numerics_solver_options.dSize;
  for (int i=0; i < dSize_min; ++i)
    numerics_solver_options.dparam[i] = dparam[i];

  int rinfo =  gfc3d_driver(problem,
                            reaction,
                            velocity,
                            globalVelocity,
                            &numerics_solver_options);


  if(outputFile == 1)
  {
    /* dump in C format */
  }
  else if (outputFile == 2)
  {
    /* dump in Numerics .dat format */
  }
  else if (outputFile == 3)
  {
#ifdef WITH_FCLIB
    fccounter++;
    if (fccounter % freq_output == 0)
    {
      char fname[256];
      snprintf(fname, sizeof(fname), "LMGC_GFC3D-i%.5d-%i-%.5d.hdf5", numerics_solver_options.iparam[7], nc, fccounter);
      printf("Dump LMGC_GFC3D-i%.5d-%i-%.5d.hdf5.\n", numerics_solver_options.iparam[7], nc, fccounter);
      /* printf("ndof = %i.\n", ndof); */

      FILE * foutput  =  fopen(fname, "w");
      int n = 100;
      char * title = (char *)malloc(n * sizeof(char *));
      strncpy(title, "LMGC dump in hdf5", n);
      char * description = (char *)malloc(n * sizeof(char *));

      snprintf(description, n, "Rewriting in hdf5 through siconos of %s in FCLIB format", fname);
      char * mathInfo = (char *)malloc(n * sizeof(char *));
      strncpy(mathInfo, "unknown", n);

      globalFrictionContact_fclib_write(problem,
                                        title,
                                        description,
                                        mathInfo,
                                        fname);


      fclose(foutput);
    }
#else
    printf("Fclib is not available ...\n");
#endif

  }


  /* NM_free(M); */
  /* NM_free(H); */
  /* free(M); */
  /* free(H); */
  free(problem);

  /* free(_colM); */
  /* free(_colH); */

  /* free(_rowM); */
  /* free(_rowH); */

  return rinfo;
}


/* int gfc3d_LmgcDriver_SBM(double *reaction, */
/*                                        double *velocity, */
/*                                        double *globalVelocity, */
/*                                        double *q, */
/*                                        double *b, */
/*                                        double *mu, */
/*                                        double *Mdata, */
/*                                        unsigned int nzM, */
/*                                        unsigned int *rowM, */
/*                                        unsigned int *colM, */
/*                                        double* Hdata, */
/*                                        unsigned int nzH, */
/*                                        unsigned int *rowH, */
/*                                        unsigned int *colH, */
/*                                        unsigned int n, */
/*                                        unsigned int nc, */
/*                                        int solver_id, */
/*                                        int isize, */
/*                                        int *iparam, */
/*                                        int dsize, */
/*                                        double *dparam, */
/*                                        int verbose, */
/*                                        int outputFile) */
/* { */

/* } */
