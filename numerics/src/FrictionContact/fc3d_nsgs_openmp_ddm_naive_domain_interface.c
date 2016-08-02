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

#include "fc3d_nsgs_openmp_ddm_naive_domain_interface.h"
#include "intersection_union.h"
#include "stdio.h"
#include "stdlib.h"

void fc3d_nsgs_openmp_ddm_naive_build_interface(unsigned int number_of_domains , unsigned int nominal_interface_size, unsigned int domain_size,  unsigned int  nc,
                                                unsigned int * interface, unsigned int * interface_size,
                                                unsigned int * interface_out,   unsigned int * interface_out_size)
{
  unsigned int n_inter=number_of_domains+1;
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

  *interface_size=cmp;
  cmp=0;
  unsigned int cmp2 =0;
  /* Warning: we assume for the moment the index of the interface component that
   * the indices if the interface are ordered.
   */

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


}

void fc3d_nsgs_openmp_ddm_naive_build_domain(unsigned int number_of_domains , unsigned int nominal_interface_size,
                                             unsigned int domain_size,  unsigned int  nc,
                                             unsigned int ** domains, unsigned int * domains_size,
                                             unsigned int ** domains_out,   unsigned int * domains_out_size)
{
    int  istart, istop;

    for (unsigned int i=0; i < number_of_domains; i++)
    {

      domains[i] = (unsigned int *)malloc((domain_size+1)*sizeof(unsigned int));

      istart = i*domain_size;
      istop = i*domain_size + domain_size;
      if (i ==  number_of_domains-1)  istop =nc;
      /* printf("id = %i \t, istart = %i\t istop = %i \n", i, istart, istop); */
      domains[i] = (unsigned int *)malloc((istop-istart)*sizeof(unsigned int));
      /* contruct index_local */
      int kk=0;
      for (int jj = istart; jj < istop; jj++ )
      {
        domains[i][kk]=jj;
        kk++;
      }
      domains_size[i] = kk;

      kk = 0;
      //printf("nc-(istop-istart) = %i \n", nc-(istop-istart));

      domains_out[i] = (unsigned int *)malloc( (nc-(istop-istart)) *sizeof(unsigned int));
      for (int jj = 0; jj < istart; jj++ )
      {
        domains_out[i][kk]=jj;
        kk++;
      }
      for (unsigned int jj = istop; jj < nc; jj++ )
      {
        domains_out[i][kk]=jj;
        kk++;
      }
      domains_out_size[i] = kk;
    }
}


void fc3d_nsgs_openmp_ddm_naive_build_domain_interface(unsigned int  nc, unsigned int number_of_domains, Ddm_domain_interface * domain_interface)
{

  unsigned int nominal_domain_size = nc/(number_of_domains); // domain size
  // estimation of interface size
  unsigned int nominal_interface_size = nominal_domain_size/number_of_domains;// nominal_domain_size/4; // nominal_domain_size/(p) ;

  if (number_of_domains==1)
    nominal_interface_size = nominal_domain_size/(number_of_domains+1);

  if (nominal_interface_size%2 == 1)
  {
    nominal_interface_size ++;
  }
  if (nominal_interface_size  == 0)
  {
    nominal_interface_size =2;
  }

  domain_interface->number_of_domains =  number_of_domains;

  domain_interface->interface_size = number_of_domains*nominal_interface_size;
  domain_interface->interface = (unsigned int*) calloc(domain_interface->interface_size,sizeof(unsigned int));
  domain_interface->interface_out_size = nc-domain_interface->interface_size;

  domain_interface->interface_out = (unsigned int*) calloc((domain_interface->interface_out_size),sizeof(unsigned int));

  fc3d_nsgs_openmp_ddm_naive_build_interface(number_of_domains ,
                                             nominal_interface_size, nominal_domain_size,  nc,
                                             domain_interface->interface, &(domain_interface->interface_size),
                                             domain_interface->interface_out, &(domain_interface->interface_out_size));

  domain_interface->domains =     (unsigned int **)malloc(number_of_domains*sizeof(unsigned int *));
  domain_interface->domains_out = (unsigned int **)malloc(number_of_domains*sizeof(unsigned int *));

  domain_interface->domains_size     = (unsigned int *)malloc(number_of_domains*sizeof(unsigned int));
  domain_interface->domains_out_size = (unsigned int *)malloc(number_of_domains*sizeof(unsigned int));

  fc3d_nsgs_openmp_ddm_naive_build_domain(number_of_domains , nominal_interface_size, nominal_domain_size, nc,
                                          domain_interface->domains, domain_interface->domains_size,
                                          domain_interface->domains_out, domain_interface->domains_out_size);
}
