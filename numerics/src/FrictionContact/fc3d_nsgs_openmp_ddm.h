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
#include <SiconosConfig.h>
#if defined(WITH_OPENMP) && defined(_OPENMP)
#define USE_OPENMP 1
#include <omp.h>
#endif

#include "FrictionContactProblem.h"
#include "SolverOptions.h"

#ifndef FC3D_NSGS_OPENMP_DDM_H
#define FC3D_NSGS_OPENMP_DDM_H


typedef struct
{
  unsigned int number_of_domains;
  unsigned int ** domains;
  unsigned int * domains_size;
  unsigned int ** domains_out;
  unsigned int * domains_out_size;
  unsigned int * interface;
  unsigned int interface_size;
  unsigned int * interface_out;
  unsigned int interface_out_size;
} Ddm_domain_interface;


void ddm_domain_interface_free(Ddm_domain_interface * domain_interface);
void ddm_domain_interface_display(Ddm_domain_interface * domain_interface);
int ddm_domain_interface_check_sorted(Ddm_domain_interface * domain_interface);
void fc3d_nsgs_openmp_ddm(FrictionContactProblem* problem, double *reaction,
                          double *velocity, int* info, SolverOptions* options,
                          Ddm_domain_interface * domain_interface);
#endif
