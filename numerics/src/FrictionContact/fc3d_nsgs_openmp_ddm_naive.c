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




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#pragma GCC diagnostic ignored "-Wmissing-prototypes"


void fc3d_nsgs_openmp_ddm_naive(FrictionContactProblem* problem, double *reaction,
                                double *velocity, int* info, SolverOptions* options)
{

  int* iparam = options->iparam;

  unsigned int max_threads = 1;
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;



#if defined(USE_OPENMP) && defined(_OPENMP)
  if (iparam[10] > 0)
  {
    max_threads = iparam[10];
    omp_set_num_threads(max_threads);
  }
  else
    max_threads = omp_get_max_threads();
#endif

  /**************************************************************
   * build contact domain sets and interface set in a naive way *
   **************************************************************/

  unsigned int nominal_domain_size = nc/(max_threads); // domain size
  // estimation of interface size
  unsigned int nominal_interface_size = nominal_domain_size/max_threads;// nominal_domain_size/4; // nominal_domain_size/(p) ;

  if (max_threads==1)
    nominal_interface_size = nominal_domain_size/(max_threads+1);

  if (nominal_interface_size%2 == 1)
  {
    nominal_interface_size ++;
  }
  if (nominal_interface_size  == 0)
  {
    nominal_interface_size =2;
  }



  printf("nominal_domain_size = %i\n", nominal_domain_size);
  printf("nc = %i\n", nc);
  printf("nominal_interface_size = %i\n", nominal_interface_size);

  unsigned int interface_size =max_threads*nominal_interface_size;
  printf("precomputed interface_size = %i\n", interface_size);

  unsigned int * interface = (unsigned int*) calloc(interface_size,sizeof(unsigned int));
  unsigned int interface_out_size = nc-interface_size;
  printf("precomputed interface_out_size = %i\n", interface_out_size);
  unsigned int * interface_out = (unsigned int*) calloc((interface_out_size),sizeof(unsigned int));

  fc3d_nsgs_openmp_ddm_naive_build_interface(max_threads , nominal_interface_size, nominal_domain_size,  nc,
                                             interface, &interface_size,
                                             interface_out,   &interface_out_size);

  unsigned int ** domains =     (unsigned int **)malloc(max_threads*sizeof(unsigned int *));
  unsigned int ** domains_out = (unsigned int **)malloc(max_threads*sizeof(unsigned int *));

  unsigned int * domains_size     = (unsigned int *)malloc(max_threads*sizeof(unsigned int));
  unsigned int * domains_out_size = (unsigned int *)malloc(max_threads*sizeof(unsigned int));

  fc3d_nsgs_openmp_ddm_naive_build_domain(max_threads , nominal_interface_size, nominal_domain_size, nc,
                                          domains, domains_size,
                                          domains_out, domains_out_size);


  /**************************************************************
   * call for ddm solver *
   **************************************************************/
  fc3d_nsgs_openmp_ddm(problem, reaction,
                       velocity, info, options,
                       max_threads,
                       domains, domains_size,
                       domains_out, domains_out_size,
                       interface, &interface_size,
                       interface_out, &interface_out_size);

}
