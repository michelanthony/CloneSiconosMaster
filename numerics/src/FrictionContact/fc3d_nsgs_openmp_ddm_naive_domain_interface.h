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

#include "fc3d_nsgs_openmp_ddm.h"


void fc3d_nsgs_openmp_ddm_naive_build_interface(unsigned int number_of_domains, unsigned int nominal_interface_size, unsigned int domain_size, unsigned int  nc,
                                                unsigned int * interface, unsigned int * interface_size,
                                                unsigned int * interface_out,   unsigned int * interface_out_size);

void fc3d_nsgs_openmp_ddm_naive_build_domain(unsigned int number_of_domains , unsigned int nominal_interface_size, unsigned int domain_size,  unsigned int  nc,
                                             unsigned int ** domains, unsigned int * domains_size,
                                             unsigned int ** domains_out,   unsigned int * domains_out_size);

void fc3d_nsgs_openmp_ddm_naive_build_domain_interface(unsigned int  nc,
                                                       unsigned number_of_domains,
                                                       Ddm_domain_interface * domain_interface);
