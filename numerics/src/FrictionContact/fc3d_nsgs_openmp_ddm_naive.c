/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_projection.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <alloca.h>
#include "op3x3.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"


#include "fc3d_nsgs_openmp.h"
void fc3d_nsgs_openmp_ddm_naive_build_interface(unsigned int max_threads, unsigned int nominal_interface_size, unsigned int domain_size, unsigned int  nc,
                                                unsigned int * interface, unsigned int * interface_size,
                                                unsigned int * interface_out,   unsigned int * interface_out_size);

void fc3d_nsgs_openmp_ddm_naive_build_interface(unsigned int max_threads , unsigned int nominal_interface_size, unsigned int domain_size,  unsigned int  nc,
                                                unsigned int * interface, unsigned int * interface_size,
                                                unsigned int * interface_out,   unsigned int * interface_out_size)
{
  unsigned int n_inter=max_threads+1;
  unsigned int cmp=0;

  for(unsigned int s =0; s< nominal_interface_size/2; s++)
  {
    interface[cmp] = s;
    cmp++;
  }

  for (unsigned int ii =1 ; ii < n_inter-1 ; ii++)
  {
    for(int s =-(int)nominal_interface_size/2 ; s<(int)nominal_interface_size/2; s++)
    {
      interface[cmp] = (ii)*domain_size+s;
      cmp++;
      if (cmp > *interface_size)
      {
        printf("problem in building interface, allocated memory is too small");
        exit(1);
      }
    }
  }

  for(int s =- (int)nominal_interface_size/2; s <0; s++)
  {
    interface[cmp] = nc+s;
    cmp++;
    if (cmp > *interface_size)
    {
      printf("problem in building interface, allocated memory is too small");
      exit(1);
    }
  }
  printf("interface_size = %i\n", cmp);

  for (unsigned int ii =0 ; ii < cmp ; ii++)
  {
    if (ii%6 ==0)  printf("\n");
    printf("interface[%i] = %i\t", ii, interface[ii]);
  }
  printf("\n");

  *interface_size=cmp;


  cmp=0;
  unsigned int cmp2 =0;
  for(unsigned int i=0; i<nc; i++)
  {
    if (i == interface[cmp])
    {
      cmp++;
    }
    else
    {
     interface_out[cmp2]=i;
     cmp2++;
     if (cmp2 > *interface_out_size)
     {
       printf("problem in building interface, allocated memory is too small");
       exit(1);
     }
    }
  }
  *interface_out_size=cmp2;
   printf("interface_out_size = %i\n", *interface_out_size);

  for (unsigned int ii =0 ; ii < *interface_out_size ; ii++)
  {
    if (ii%6 ==0)  printf("\n");
    printf("interface_out[%i] = %i\t", ii, interface_out[ii]);
  }
  printf("\n");
}
void fc3d_nsgs_openmp_ddm_naive_build_domain(unsigned int max_threads , unsigned int nominal_interface_size, unsigned int domain_size,  unsigned int  nc,
                                             unsigned int ** domains, unsigned int * domains_size,
                                             unsigned int ** domains_out,   unsigned int * domains_out_size);
void fc3d_nsgs_openmp_ddm_naive_build_domain(unsigned int max_threads , unsigned int nominal_interface_size, unsigned int domain_size,  unsigned int  nc,
                                             unsigned int ** domains, unsigned int * domains_size,
                                             unsigned int ** domains_out,   unsigned int * domains_out_size)
{
    int  istart, istop;

    for (unsigned int i=0; i < max_threads; i++)
    {

      domains[i] = (unsigned int *)malloc((domain_size+1)*sizeof(unsigned int));

      istart = i*domain_size;
      istop = i*domain_size + domain_size;
      if (i ==  max_threads-1)  istop =nc;
      printf("id = %i \t, istart = %i\t istop = %i \n", i, istart, istop);
      domains[i] = (unsigned int *)malloc((istop-istart)*sizeof(unsigned int));
      printf("domains building for i = %i\n", i);

      /* contruct index_local */
      int kk=0;
      for (int jj = istart; jj < istop; jj++ )
      {
        domains[i][kk]=jj;
        if (jj%6 ==0)  printf("\n");
        printf("domains[%i][%i] = %i\t",i,kk,domains[i][kk]);
        kk++;
      }

      domains_size[i] = kk;
      printf("\n domains_size[%i] = %i\n",i,domains_size[i]);

      kk = 0;
      printf("nc-(istop-istart) = %i \n", nc-(istop-istart));
      domains_out[i] = (unsigned int *)malloc( (nc-(istop-istart)) *sizeof(unsigned int));

      for (int jj = 0; jj < istart; jj++ )
      {
        domains_out[i][kk]=jj;
        kk++;
      }
      for (unsigned int jj = istop; jj < nc; jj++ )
      {
        domains_out[i][kk]=jj;
        if (jj%6 ==0)  printf("\n");
        printf("threadout__index[%i][%i] = %i\t",i,kk,domains_out[i][kk]);
        kk++;
      }
      domains_out_size[i] = kk;
      printf("\n domains_out_size[%i] = %i\n",i,domains_out_size[i]);

    }
    //getchar();
}

void fc3d_nsgs_openmp_ddm_naive(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options)
{


#if defined(USE_OPENMP) && defined(_OPENMP)

#else
  printf("fc3d_nsgs_openmp_ddm_stupid cannot be used without openmp");
  exit(1);
#endif

  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);
  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs_redblack_openmp", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;

  SolverPtr local_solver = NULL;
  Update_indexPtr update_domain_problem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Allocate space for local solver and local problem */
  unsigned int max_threads = 1;
#if defined(USE_OPENMP) && defined(_OPENMP)
  if (iparam[10] > 0)
  {
    max_threads = iparam[10];
    omp_set_num_threads(max_threads);
  }
  else
    max_threads = omp_get_max_threads();
  FrictionContactProblem **domain_problems = alloca(max_threads*sizeof(void*));
  SolverOptions **domain_solver_options = alloca(max_threads*sizeof(void*));
  FrictionContactProblem **interface_local_problems = alloca(max_threads*sizeof(void*));
  SolverOptions **interface_local_solver_options = alloca(max_threads*sizeof(void*));
#else
  FrictionContactProblem *domain_problems[1];
  SolverOptions *domain_solver_options[1];
  FrictionContactProblem *interface_local_problems[1];
  SolverOptions *interface_local_solver_options[1];
#endif


  if (verbose > 0) printf("----------------------------------- number of threads %i\n", omp_get_max_threads()  );
  if (verbose > 0) printf("----------------------------------- number of contacts %i\n", nc );

  double * q_k = (double *) malloc(nc*3*sizeof(double));

  int domain_itermax=options->iparam[12], interface_itermax=options->iparam[13],  domain_iter_total=0;


  for (unsigned int i=0; i < max_threads; i++)
  {
    printf(" initilialization of interface_local_problem and local solver options\n");
    FrictionContactProblem *interface_local_problem = malloc(sizeof(FrictionContactProblem));
    interface_local_problems[i] = interface_local_problem;
    interface_local_problem->numberOfContacts = 1;
    interface_local_problem->dimension = 3;
    interface_local_problem->q = (double*)malloc(3 * sizeof(double));
    interface_local_problem->mu = (double*)malloc(sizeof(double));

    if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
    {
      interface_local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                                      malloc(9 * sizeof(double)));
    }
    else /* NM_SPARSE_BLOCK */
    {
      interface_local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
    }

    interface_local_solver_options[i] = malloc(sizeof(SolverOptions));
    solver_options_nullify(interface_local_solver_options[i]);
    interface_local_solver_options[i]->dparam = NULL;
    interface_local_solver_options[i]->iparam = NULL;
    solver_options_copy(localsolver_options,interface_local_solver_options[i]);

    fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_domain_problem,
                                            (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                            problem, interface_local_problems[i],
                                            options, interface_local_solver_options[i]);
    printf(" initilialization of domain_problem and thread solver options\n");

    FrictionContactProblem *domain_problem = malloc(sizeof(FrictionContactProblem));
    domain_problems[i] = domain_problem;
    domain_problem->M = problem->M;
    domain_problem->numberOfContacts = nc;
    domain_problem->dimension = 3;
    domain_problem->q = q_k;
    domain_problem->mu = problem->mu;;


    /* printSolverOptions(options); */
    /* getchar(); */

    domain_solver_options[i] = malloc(sizeof(SolverOptions));
    solver_options_nullify(domain_solver_options[i]);
    solver_options_copy(options, domain_solver_options[i]);
    domain_solver_options[i]->dparam[0] /= 10.0;

    domain_solver_options[i]->iparam[0]=domain_itermax;
    domain_solver_options[i]->iparam[1]=1; // light error

    /* fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_domain_problem, */
    /*                                   (FreeSolverNSGSPtr *)&freeSolver, &computeError, */
    /*                                   problem, domain_problems[i], */
    /*                                   options, domain_solver_options[i]); */

  }


    FrictionContactProblem *interface_problem = malloc(sizeof(FrictionContactProblem));
    interface_problem->M = problem->M;
    interface_problem->numberOfContacts = nc;
    interface_problem->dimension = 3;
    interface_problem->q = q_k;
    interface_problem->mu = problem->mu;;


    /***********************************************
     * build contact domain sets and interface set *
     ***********************************************/

    unsigned int nominal_domain_size = nc/(max_threads); // domain size
    // estimation of interface size
    unsigned int nominal_nominal_interface_size = nominal_domain_size/max_threads;// nominal_domain_size/4; // nominal_domain_size/(p) ;

    if (max_threads==1)
      nominal_nominal_interface_size = nominal_domain_size/(max_threads+1);

    if (nominal_nominal_interface_size%2 == 1)
    {
      nominal_nominal_interface_size ++;
    }

    printf("nominal_domain_size = %i\n", nominal_domain_size);
    printf("nc = %i\n", nc);
    printf("nominal_nominal_interface_size = %i\n", nominal_nominal_interface_size);

    unsigned int interface_size =max_threads*nominal_nominal_interface_size;
    printf("precomputed interface_size = %i\n", interface_size);
    unsigned int * interface = (unsigned int*) calloc(interface_size,sizeof(unsigned int));
    unsigned int interface_out_size = nc-interface_size;
    printf("precomputed interface_out_size = %i\n", interface_out_size);
    unsigned int * interface_out = (unsigned int*) calloc((interface_out_size),sizeof(unsigned int));

    fc3d_nsgs_openmp_ddm_naive_build_interface(max_threads , nominal_nominal_interface_size, nominal_domain_size,  nc,
                                               interface, &interface_size,
                                               interface_out,   &interface_out_size);

    unsigned int ** domains =     (unsigned int **)malloc(max_threads*sizeof(unsigned int *));
    unsigned int ** domains_out = (unsigned int **)malloc(max_threads*sizeof(unsigned int *));

    unsigned int * domains_size     = (unsigned int *)malloc(max_threads*sizeof(unsigned int));
    unsigned int * domains_out_size = (unsigned int *)malloc(max_threads*sizeof(unsigned int));

    fc3d_nsgs_openmp_ddm_naive_build_domain(max_threads , nominal_nominal_interface_size, nominal_domain_size, nc,
                                            domains, domains_size,
                                            domains_out, domains_out_size);

    /**********************************************
     *****  NSGS Iterations ****
     **********************************************/
    int iter = 0; /* Current iteration number */
    double error = 1.; /* Current error */
    int hasNotConverged = 1;
    double error_delta_reaction=0.0;
    double error_nat=0.0;


  double * reaction_k = (double*)malloc(nc*3*sizeof(double));
  double * velocity_k = (double*)malloc(nc*3*sizeof(double));
  for (unsigned int i =0 ; i < 3*nc; i++)    reaction_k[i] = reaction[i];


  double normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1);
  printf("normreaction_k = %e \n", normreaction_k);

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    if (verbose > 0) printf(" \n ################### START GLOBAL LOOP ################################ \n");

    ++iter;

    { /* thread block */
      error_delta_reaction=0.0;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      /* for (unsigned int kk=0; kk < 3*nc; kk++ ) reaction_k[kk]=reaction[kk]; */

      #pragma omp  parallel                           \
        shared(reaction, reaction_k, velocity_k, q_k, \
               domain_problems, domain_solver_options,\
               domains, domains_size,       \
               error_delta_reaction, domain_iter_total)
      {
        int id = omp_get_thread_num();
        int domain_iter=0;
        int domain_info =1;
        for (unsigned int i = 0; i < domains_size[id]; i++ )
        {
          int contact = domains[id][i];
          fc3d_nsgs_index_computeqLocal(problem, reaction, contact,
                                        domains_out[id], domains_out_size[id],
                                        &(q_k[3*contact]) );
        }
        /* for (unsigned int i =0 ; i < 3*nc; i++) q_k[i] = problem->q[i]; */
        /* double normq_k = cblas_dnrm2(nc*3 , q_k , 1); */
        /* printf("############ normq = %e\t normq_k = %e\n", normq, normq_k); */

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        fc3d_nsgs_index(domain_problems[id],
                        reaction_k, velocity_k,
                        &domain_info, domain_solver_options[id],
                        domains[id], domains_size[id]);

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        /* fc3d_nsgs(domain_problems[id], */
        /*           reaction_k, velocity_k, */
        /*           &domain_info, domain_solver_options[id]); */
        domain_iter = domain_solver_options[id]->iparam[7];
        #pragma omp critical
        domain_iter_total += domain_iter;
        error_delta_reaction +=  domain_solver_options[id]->dparam[1];;
      }

      /* normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1); */
      /* printf("################### normreaction_k = %e \n", normreaction_k); */

      if (verbose > 0) printf("----------------------------------- FC3D - NSGS DDM NAIVE - End of thread problems after %i iterations with error_delta_reaction =%e \n", domain_iter_total, error_delta_reaction);
    } /* end of thread block */


    /* ------------------------------------------------------- */
    /* ----------------- interface loop ---------------------- */
    /* ------------------------------------------------------- */
    {
      double error_delta_reaction_interface=0.0;
      for (unsigned int i = 0; i < interface_size; i++ )
      {
        int  contact = interface[i];
        q_k[3*contact]=0.0;
        q_k[3*contact+1]=0.0;
        q_k[3*contact+2]=0.0;
        fc3d_nsgs_index_computeqLocal(problem, reaction_k, contact,
                                      interface_out, interface_out_size,
                                      &(q_k[3*contact]));
      }
      /* double normq_k = cblas_dnrm2(nc*3 , q_k , 1); */
      /* printf("normq = %e\t normq_k = %e\n", normq, normq_k); */

      int interface_hasNotConverged =1;
      int interface_iter=0;

      while ((interface_iter < interface_itermax) && ( interface_hasNotConverged > 0 ))
      {
        ++interface_iter;
        error_delta_reaction_interface=0.0;
        #pragma omp parallel for reduction(+:error_delta_reaction_interface)
        for ( unsigned int i = 0 ; i < interface_size ; i++)
        {
          unsigned int tid = omp_get_thread_num();
          int contact = interface[i];

          if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);

          (*update_domain_problem)(contact, interface_problem, interface_local_problems[tid],
                                   reaction_k,
                                   interface_local_solver_options[tid],
                                   interface, interface_size);

          interface_local_solver_options[tid]->iparam[4] = contact;

          /* version without localreaction */
          double localreaction[3];
          {
            localreaction[0] = reaction_k[3 * contact+0];
            localreaction[1] = reaction_k[3 * contact+1];
            localreaction[2] = reaction_k[3 * contact+2];
          };

          (*local_solver)(interface_local_problems[tid], localreaction,
                          interface_local_solver_options[tid]);

          /* for( int kkk=0; kkk <3 ; kkk++) printf("localreaction[%i] = %4.2e\n", kkk,localreaction[kkk] ); */
          {
            error_delta_reaction_interface += pow(reaction_k[3 * contact] - localreaction[0], 2) +
              pow(reaction_k[3 * contact + 1] - localreaction[1], 2) +
              pow(reaction_k[3 * contact + 2] - localreaction[2], 2);

            reaction_k[3 * contact+0] = localreaction[0];
            reaction_k[3 * contact+1] = localreaction[1];
            reaction_k[3 * contact+2] = localreaction[2];

          }
        }

        /* if (error_delta_reaction_interface < tolerance) interface_hasNotConverged = 0;*/
        error_delta_reaction_interface = sqrt(error_delta_reaction_interface);
        double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
        if (fabs(norm_r) > DBL_EPSILON)
          error_delta_reaction_interface /= norm_r;

        if (error_delta_reaction_interface < tolerance)
        {
          interface_hasNotConverged = 0;
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS DDM NAIVE interface - Iteration %i Residual = %14.7e < %7.3e\n", interface_iter, error_delta_reaction_interface, options->dparam[0]);
        }
        else
        {
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS DDM NAIVE interface - Iteration %i Residual = %14.7e > %7.3e\n", interface_iter, error_delta_reaction_interface, options->dparam[0]);
        }

      }

      error_delta_reaction+= error_delta_reaction_interface;

      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM NAIVE iter = %i,  error_delta_reaction = %14.7e \n", iter, error_delta_reaction);
      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM NAIVE iter = %i,  error_delta_reaction_interface = %14.7e \n", iter, error_delta_reaction_interface);
      //double normq = cblas_dnrm2(nc*3 , problem->q , 1);

      error_delta_reaction = sqrt(error_delta_reaction);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);

      if (fabs(norm_r) > DBL_EPSILON)
        error_delta_reaction /= norm_r;


    }
    /* ----------------- end of interface loop --------------- */


    for (unsigned int i =0 ; i < 3*nc; i++)      reaction[i] = reaction_k[i];


    error = 0.0;
    double normq = cblas_dnrm2(nc*3 , problem->q , 1);
    (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS DDM NAIVE - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
        double normq = cblas_dnrm2(nc*3 , problem->q , 1);
        (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error_nat);
        printf("----------------------------------- FC3D - NSGS DDM NAIVE - Iteration %i Full Residual = %14.7e \n", iter, error_nat);
        /* test of consistency */
        double c  = 10.0;
        if   ( (error_nat/c >=  error_delta_reaction) || (error_delta_reaction >= c *error_nat))
        {
          printf("%e %e %e   \n",error_nat/c, error_delta_reaction, c *error_nat     );
          printf(" WARNING: rel error_delta_reaction is not consistent with natural map error  \n");
        }

      }
    }
    else
    {
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS DDM NAIVE- Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
      }
    }

    *info = hasNotConverged;

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, 3 * nc,
                                               reaction, velocity,
                                               error, NULL);
    }
  }

  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  for (unsigned int i=0; i < max_threads; i++)
  {
    (*freeSolver)(problem,interface_local_problems[i],interface_local_solver_options[i]);
    if (problem->M->storageType == NM_DENSE && interface_local_problems[i]->M->matrix0)
    {
      free(interface_local_problems[i]->M->matrix0);
    }
    interface_local_problems[i]->M->matrix0 = NULL;
    freeFrictionContactProblem(interface_local_problems[i]);
    solver_options_delete(interface_local_solver_options[i]);
    free(interface_local_solver_options[i]);
  }


}
