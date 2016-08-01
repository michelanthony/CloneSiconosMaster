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
#include "fc3d_Solvers.h"
#include "fc3d_nsgs_openmp.h"
#include "fc3d_nsgs_openmp_ddm_naive_domain_interface.h"
#include "SiconosBlas.h"



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <alloca.h>

#pragma GCC diagnostic ignored "-Wmissing-prototypes"



//#define USE_OPENMP_HERE 1


void fc3d_nsgs_openmp_ddm(FrictionContactProblem* problem, double *reaction,
                          double *velocity, int* info, SolverOptions* options,
                          unsigned int max_threads,
                          unsigned int ** domains, unsigned int * domains_size,
                          unsigned int ** domains_out,   unsigned int * domains_out_size,
                          unsigned int * interface, unsigned int * interface_size,
                          unsigned int * interface_out,   unsigned int * interface_out_size)
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
#if defined(USE_OPENMP) && defined(_OPENMP)
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
    if (verbose > 0) printf(" initilialization of interface_local_problem and local solver options\n");
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

    fc3d_nsgs_domain_initialize_local_solver(&local_solver, &update_domain_problem,
                                            (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                            problem, interface_local_problems[i],
                                            options, interface_local_solver_options[i]);
    if (verbose > 0) printf(" initilialization of domain_problem and thread solver options\n");

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

    /* the choice of the accuracy for the domain solver has to be studied in more details */
    domain_solver_options[i]->dparam[0] /= 10.0;
    /* domain_solver_options[i]->dparam[0] /= max_threads; */

    domain_solver_options[i]->iparam[0]=domain_itermax;
    domain_solver_options[i]->iparam[1]=
      SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION_WITH_FULL_FINAL; // light with full final  error
    domain_solver_options[i]->iparam[1]=
      SICONOS_FRICTION_3D_NSGS_LIGHT_ERROR_EVALUATION; // light error

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


  /* double normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1); */
  /* printf("normreaction_k = %e \n", normreaction_k); */

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    if (verbose > 0) printf("\n----------------------------------- FC3D - NSGS DDM - Start global loop \n");
    ++iter;

    { /* thread block */
      error_delta_reaction=0.0;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      /* for (unsigned int kk=0; kk < 3*nc; kk++ ) reaction_k[kk]=reaction[kk]; */
#ifdef USE_OPENMP_HERE
      #pragma omp  parallel                             \
        shared(reaction, reaction_k, velocity_k, q_k,   \
               domain_problems, domain_solver_options,  \
               domains, domains_size,                   \
               error_delta_reaction, domain_iter_total)
      {
        int id = omp_get_thread_num();
#else
        for (int id =0; id < (int)max_threads; id++)
      {
#endif
        int domain_iter =0;
        int domain_info =1;
        for (unsigned int i = 0; i < domains_size[id]; i++ )
        {
          int contact = domains[id][i];
          fc3d_nsgs_domain_computeqLocal(problem, reaction, contact,
                                        domains_out[id], domains_out_size[id],
                                        &(q_k[3*contact]));
        }
        /* for (unsigned int i =0 ; i < 3*nc; i++) q_k[i] = problem->q[i]; */
        /* double normq_k = cblas_dnrm2(nc*3 , q_k , 1); */
        /* printf("############ normq = %e\t normq_k = %e\n", normq, normq_k); */

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        fc3d_nsgs_domain(domain_problems[id],
                        reaction_k, velocity_k,
                        &domain_info, domain_solver_options[id],
                        domains[id], domains_size[id]);

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        /* fc3d_nsgs(domain_problems[id], */
        /*           reaction_k, velocity_k, */
        /*           &domain_info, domain_solver_options[id]); */
        domain_iter = domain_solver_options[id]->iparam[7];
#ifdef USE_OPENMP_HERE
        #pragma omp critical
#endif
        {
        domain_iter_total += domain_iter;
        error_delta_reaction +=  domain_solver_options[id]->dparam[2];
        }
      }

      /* normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1); */
      /* printf("################### normreaction_k = %e \n", normreaction_k); */
      /* double normq = cblas_dnrm2(nc*3 , problem->q , 1); */
      /* printf("normq = %e\n",normq); */

     if (verbose > 0)
       printf("----------------------------------- FC3D - NSGS DDM  - End of domain problems after %i iterations\n",
              domain_iter_total);
    } /* end of thread block */


    /* ------------------------------------------------------- */
    /* ----------------- interface loop ---------------------- */
    /* ------------------------------------------------------- */
    {
      double error_delta_reaction_interface=0.0, error_interface=0.0;
      for (unsigned int i = 0; i < *interface_size; i++ )
      {
        int  contact = interface[i];
        fc3d_nsgs_domain_computeqLocal(problem, reaction_k, contact,
                                      interface_out, *interface_out_size,
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
#ifdef USE_OPENMP_HERE
         #pragma omp parallel for reduction(+:error_delta_reaction_interface)
#endif
        for ( unsigned int i = 0 ; i < *interface_size ; i++)
        {
          unsigned int tid = omp_get_thread_num();
          int contact = interface[i];

          if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);

          (*update_domain_problem)(contact, interface_problem, interface_local_problems[tid],
                                   reaction_k,
                                   interface_local_solver_options[tid],
                                   interface, *interface_size);

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
        error_interface = sqrt(error_delta_reaction_interface);
        double norm_r = cblas_dnrm2(nc*3 , reaction , 1);

        if (fabs(norm_r) > DBL_EPSILON)
          error_interface /= norm_r;

        if (error_interface < tolerance)
        {
          interface_hasNotConverged = 0;
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS INTERFACE - Iteration %i Residual = %14.7e < %7.3e\n",
                   interface_iter, error_interface, options->dparam[0]);
        }
        else
        {
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS INTERFACE - Iteration %i Residual = %14.7e > %7.3e\n",
                   interface_iter, error_interface, options->dparam[0]);
        }
      } /* ----------------- end of interface loop --------------- */


      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM  error_delta_reaction = %14.7e \n",
                               error_delta_reaction);
      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM  error_delta_reaction_interface = %14.7e \n",
                               error_delta_reaction_interface);




      error_delta_reaction+= error_delta_reaction_interface;
     //double normq = cblas_dnrm2(nc*3 , problem->q , 1);

      error_delta_reaction = sqrt(error_delta_reaction);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);

      if (fabs(norm_r) > DBL_EPSILON)
        error_delta_reaction /= norm_r;

      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM iter = %i,  error_delta_reaction = %14.7e \n",
                               iter, error_delta_reaction);


    } /* ----------------- end of interface block --------------- */


    /* normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1); */
    /* printf("normreaction_k = %e \n", normreaction_k); */

    for (unsigned int i =0 ; i < 3*nc; i++)     reaction[i] = reaction_k[i];

    /* error = 0.0; */
    /* double normq = cblas_dnrm2(nc*3 , problem->q , 1); */
    /* (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error); */
    error =  error_delta_reaction;


    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS DDM - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
        double normq = cblas_dnrm2(nc*3 , problem->q , 1);
        //printf("normq = %e\n",normq);
        (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error_nat);
        printf("----------------------------------- FC3D - NSGS DDM - Iteration %i Full Residual = %14.7e \n", iter, error_nat);
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
        printf("----------------------------------- FC3D - NSGS DDM - Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
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
  /* ----------------- end of the global loop --------------- */

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
