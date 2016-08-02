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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SiconosBlas.h"
#include "intersection_union.h"
int main(void)
{

  printf("========= Starts SBM tests 6 for SBM ========= \n");
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM3x3.dat", "r");
  newFromFileSBM(&M, file);
  //printSBM(&M);
  fclose(file);
  /*alloc enough memory */

  int res = 1 ;

  //int res = test_SBMRowToDense(&M);

  unsigned int i;
  unsigned int m = M.blocksize1[M.blocknumber1 - 1];
  double * q = (double *)malloc(m * sizeof(double));
  double * y = (double *)malloc(m * sizeof(double));
  double * yref = (double *)malloc(m * sizeof(double));
  double * ytmp = (double *)malloc(m * sizeof(double));
  double * y3 = (double *)malloc(3 * sizeof(double));
  for(unsigned int j =0; j<m; j++) yref[j]=0.0;

  for (i = 0; i < m; i++)
  {
    q[i] = i + 1;
  }
  printf("test rowProdNoDiagSBM\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM(m, 3, blockRow, &M, q, &yref[3*blockRow],0);
  }
  for(int j =0; j<6; j++) printf("yref[%i] = %e\n", j, yref[j]);
  double normref = cblas_dnrm2(m , yref , 1);
  printf("normref  = %e \n", normref);

  printf("\ntest rowProdNoDiagSBM3x3 on row Block  0\n");
  for(unsigned int j =0; j<3; j++) y3[j]=0.0;
  rowProdNoDiagSBM3x3(m, 3, 0, &M, q, y3);
  for(unsigned int j =0; j<3; j++) printf("y3[%i] = %e\n", j, y3[j]);


  printf("\ntest rowProdNoDiagSBM\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM(m, 3, blockRow, &M, q, &y[3*blockRow],0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);


  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  double normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);




  printf("\ntest rowProdNoDiagSBM3x3 on all row Blocks\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM3x3(m, 3, blockRow, &M, q, &y[3*blockRow]);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);
  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  for (i = 0; i < m; i++)
  {
    yref[i] = y[i];
  }




  printf("\ntest rowProdNoDiagSBM3x3_index_block for full index\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  unsigned int * index = (unsigned int *) malloc(M.blocknumber0*sizeof(unsigned int));
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    index[blockRow] = blockRow;
  }
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM3x3_index_block(m, 3, blockRow, &M, q, &y[3*blockRow],index,M.blocknumber0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);
  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);


  printf("\ntest rowProdNoDiagSBM3x3_index_block for empty index\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM3x3_index_block(m, 3, blockRow, &M, q, &y[3*blockRow],index,0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);


  printf("\ntest rowProdNoDiagSBM3x3_index_block for a given index\n");

  for(unsigned int blockRow =0; blockRow < 10; blockRow++)
  {
    index[blockRow] = blockRow+20;
  }
  printf("index= ");
  uint_array_print(index, 10);

  unsigned int * index_out = (unsigned int *) malloc(M.blocknumber0*sizeof(unsigned int));
  unsigned index_out_size = 0;
  unsigned int cmp =0;
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    if (blockRow == index[cmp])
    {
      cmp++;
    }
    else
    {
      index_out[index_out_size] = blockRow ;
      index_out_size++;

    }
  }

  printf("index_out= ");
  uint_array_print(index_out, index_out_size );
  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM3x3_index_block(m, 3, blockRow, &M, q, &y[3*blockRow],index,10);
  }
  for(int j =3*20; j<3*20+6; j++) printf("y[%i] = %e\n", j, y[j]);

  for(unsigned int blockRow =0; blockRow < M.blocknumber0; blockRow++)
  {
    rowProdNoDiagSBM3x3_index_block(m, 3, blockRow, &M, q, &y[3*blockRow],index_out,index_out_size);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);

  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);

  if (res)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    return 1;
  }

  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);

}
