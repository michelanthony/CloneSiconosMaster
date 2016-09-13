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
#include <stdio.h>
#include <stdlib.h>
#include "fc3d_Solvers.h"
#include "fc3d_nsgs_openmp.h"
#include "fc3d_nsgs_openmp_ddm.h"
#include "fc3d_nsgs_openmp_ddm_naive_domain_interface.h"


void fc3d_nsgs_openmp(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;

  unsigned int max_threads = 1;

  #if defined(USE_OPENMP) && defined(_OPENMP)
  if (options->iparam[10] > 0)
  {
    omp_set_num_threads(options->iparam[10]);
  }
  max_threads = omp_get_max_threads();
  #endif

  if (iparam[11] == 0)
  {
    fc3d_nsgs_openmp_for(problem, reaction, velocity, info, options) ;
  }
  else if (iparam[11] == 1)
  {
    fc3d_nsgs_openmp_redblack(problem, reaction, velocity, info, options) ;
  }
  else if (iparam[11] == 2)
  {
    /* Number of contacts */
    unsigned int nc = problem->numberOfContacts;

    /* build contact domain sets and interface set in a naive way */
    Ddm_domain_interface domain_interface;
    fc3d_nsgs_openmp_ddm_naive_build_domain_interface(nc, max_threads, &domain_interface);

    ddm_domain_interface_display(&domain_interface);
    if (ddm_domain_interface_check_sorted(&domain_interface))
      numericsError("fc3d_nsgs_openmp_ddm", "The method needs sorted index for domain and interface definition");
      
    /* call for ddm solver */
    fc3d_nsgs_openmp_ddm(problem, reaction, velocity, info, options,
                         &domain_interface);

    ddm_domain_interface_free(&domain_interface);
  }
  else if (iparam[11] == 3)
  {
    fc3d_nsgs_openmp_iterfor(problem, reaction, velocity, info, options) ;
  }
  else if (iparam[11] == 10)
  {
    fc3d_nsgs_error_comparison(problem, reaction, velocity, info, options) ;
  }
  else
  {
    numericsError("fc3d_nsgs_openmp", "The method defined by iparam[11] is not recognized");
  }

}



int fc3d_nsgs_openmp_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS OPENMP Solver\n");
  }

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_FRICTION_3D_NSGS_OPENMP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 15;
  options->dSize = 15;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  for (i = 0; i < 15; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
