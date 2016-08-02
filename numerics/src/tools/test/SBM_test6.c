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
  double tolerance= 1e-9;
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
  
  SparseBlockStructuredMatrix M_copy;
  copySBM(&M,&M_copy, 1);
  //printSBM(&M_copy);
  size_t colNumber = 0;

  for(unsigned int currentRowNumber =0; currentRowNumber < M_copy.blocknumber0; currentRowNumber++)
  {
    /* nullify diagional block */
    for (size_t blockNum = M_copy.index1_data[currentRowNumber];
         blockNum < M_copy.index1_data[currentRowNumber + 1];
         ++blockNum)
    {
      /* Get row/column position of the current block */
      colNumber = M_copy.index2_data[blockNum];
      if (colNumber == currentRowNumber)
      {
        for (i = 0; i < 9; i++)      M_copy.block[blockNum][i]=0.0;
      }
    }
  }
  
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag(m, 3, currentRowNumber, &M, q, &yref[3*currentRowNumber],0);
  }
  for(int j =0; j<6; j++) printf("yref[%i] = %e\n", j, yref[j]);
  double normref = cblas_dnrm2(m , yref , 1);
  printf("normref  = %e \n", normref);

  
  printf("test SBM_row_prod_no_diag w.r.t SBM_row_prod\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M_copy.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod(m, 3, currentRowNumber, &M_copy, q, &y[3*currentRowNumber],0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);
  for (i = 0; i < 3; i++)    ytmp[i] = y[i] - yref[i];
  double normtmp = cblas_dnrm2(3 , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);

  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  


  
  printf("\ntest SBM_row_prod_no_diag_3x3 on row Block  0\n");
  for(unsigned int j =0; j<3; j++) y3[j]=0.0;
  SBM_row_prod_no_diag_3x3(m, 3, 0, &M, q, y3);
  for(unsigned int j =0; j<3; j++) printf("y3[%i] = %e\n", j, y3[j]);
  for (i = 0; i < 3; i++)
  {
    ytmp[i] = y3[i] - yref[i];
  }
  normtmp = cblas_dnrm2(3 , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);

  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  
  printf("\ntest SBM_row_prod_no_diag\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber],0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);


  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }



  printf("\ntest SBM_row_prod_no_diag_3x3 on all row Blocks\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber]);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);
  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  


  
  for (i = 0; i < m; i++)
  {
    yref[i] = y[i];
  }




  printf("\ntest SBM_row_prod_no_diag_3x3_index_block for full index\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  unsigned int * index = (unsigned int *) malloc(M.blocknumber0*sizeof(unsigned int));
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    index[currentRowNumber] = currentRowNumber;
  }
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber],index,M.blocknumber0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);
  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  

  printf("\ntest SBM_row_prod_no_diag_3x3_index_block for empty index\n");
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber],index,0);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);


  printf("\ntest SBM_row_prod_no_diag_3x3_index_block for a given index\n");

  for(unsigned int currentRowNumber =0; currentRowNumber < 10; currentRowNumber++)
  {
    index[currentRowNumber] = currentRowNumber+20;
  }
  printf("index= ");
  uint_array_print(index, 10);

  unsigned int * index_out = (unsigned int *) malloc(M.blocknumber0*sizeof(unsigned int));
  unsigned index_out_size = 0;
  unsigned int cmp =0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    if (currentRowNumber == index[cmp])
    {
      cmp++;
    }
    else
    {
      index_out[index_out_size] = currentRowNumber ;
      index_out_size++;

    }
  }

  printf("index_out= ");
  uint_array_print(index_out, index_out_size );
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber],index,10);
  }
  for(int j =3*20; j<3*20+6; j++) printf("y[%i] = %e\n", j, y[j]);

  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &y[3*currentRowNumber],index_out,index_out_size);
  }
  for(int j =0; j<6; j++) printf("y[%i] = %e\n", j, y[j]);

  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  
  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }

  
  printf("\ntest SBM_row_prod_no_diag_3x3_index_block for a given index w.r.t. SBM_row_prod_no_diag_3x3\n");

  
  SparseBlockStructuredMatrix M_index;
  copySBM(&M,&M_index, 1);
  //printSBM(&M_copy);

  for(unsigned int currentRowNumber =0; currentRowNumber < M_index.blocknumber0; currentRowNumber++)
  {
    /* printf("--- currentRowNumber = %i\n", currentRowNumber); */
    /* nullify block in column that corresponds to index */
    for (size_t blockNum = M_index.index1_data[currentRowNumber];
         blockNum < M_index.index1_data[currentRowNumber + 1];
         ++blockNum)
    {
      /* Get row/column position of the current block */
      colNumber = M_index.index2_data[blockNum];
      /* printf("colnumber = %zu\n", colNumber ); */
      int inside_index =0;
      for (int i = 0 ; i <10 ; i++)
      {
        if (colNumber == index[i])
        {
          inside_index=1;
        }
      }
      if (inside_index)
      {
        /* printf("keep block number = %i, i = %i \n", blockNum); */
      }
      else
      {
        /* printf("nullify block number = %i \n", blockNum); */
        for (i = 0; i < 9; i++)      M_index.block[blockNum][i]=0.0;
      }
    }
  }
  double * yref2 = (double *)malloc(m * sizeof(double));
  for(unsigned int j =0; j<m; j++) yref2[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    /* printf("--- currentRowNumber = %i\n", currentRowNumber); */
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &yref2[3*currentRowNumber],index,10);
  }
  for(int j =3*20; j<3*20+6; j++) printf("yref2[%i] = %e\n", j, yref2[j]);


  
  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M_index.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3(m, 3, currentRowNumber, &M_index, q, &y[3*currentRowNumber]);
  }
  for(int j =3*20; j<3*20+6; j++) printf("y[%i] = %e\n", j, y[j]);

  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref2[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);

  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  printf("\ntest SBM_row_prod_no_diag_3x3_index_block for a given index w.r.t. SBM_row_prod_no_diag_3x3 index_out\n");

  
  SparseBlockStructuredMatrix M_index_out;
  copySBM(&M,&M_index_out, 1);
  //printSBM(&M_copy);

  for(unsigned int currentRowNumber =0; currentRowNumber < M_index_out.blocknumber0; currentRowNumber++)
  {
    /* printf("--- currentRowNumber = %i\n", currentRowNumber); */
    /* nullify block in column that corresponds to index */
    for (size_t blockNum = M_index_out.index1_data[currentRowNumber];
         blockNum < M_index_out.index1_data[currentRowNumber + 1];
         ++blockNum)
    {
      /* Get row/column position of the current block */
      colNumber = M_index_out.index2_data[blockNum];
      /* printf("colnumber = %zu\n", colNumber ); */
      int inside_index_out =0;
      for (int i = 0 ; i <index_out_size ; i++)
      {
        if (colNumber == index_out[i])
        {
          inside_index_out=1;
        }
      }
      if (inside_index_out)
      {
        /* printf("keep block number = %i, i = %i \n", blockNum); */
      }
      else
      {
        /* printf("nullify block number = %i \n", blockNum); */
        for (i = 0; i < 9; i++)      M_index_out.block[blockNum][i]=0.0;
      }
    }
  }
  
  double * yref3 = (double *)malloc(m * sizeof(double));
  for(unsigned int j =0; j<m; j++) yref3[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M.blocknumber0; currentRowNumber++)
  {
    /* printf("--- currentRowNumber = %i\n", currentRowNumber); */
    SBM_row_prod_no_diag_3x3_index_block(m, 3, currentRowNumber, &M, q, &yref3[3*currentRowNumber],index_out,index_out_size);
  }
  for(int j =3*20; j<3*20+6; j++) printf("yref3[%i] = %e\n", j, yref3[j]);

  for(unsigned int j =0; j<m; j++) y[j]=0.0;
  for(unsigned int currentRowNumber =0; currentRowNumber < M_index.blocknumber0; currentRowNumber++)
  {
    SBM_row_prod_no_diag_3x3(m, 3, currentRowNumber, &M_index_out, q, &y[3*currentRowNumber]);
  }
  for(int j =3*20; j<3*20+6; j++) printf("y[%i] = %e\n", j, y[j]);

  for (i = 0; i < m; i++)
  {
    ytmp[i] = y[i] - yref3[i];
  }

  normtmp = cblas_dnrm2(m , ytmp , 1);
  printf("normtmp  = %e \n", normtmp);
  if (normtmp > tolerance)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    exit(1);
  }
  res=0 ;
  if (res)
  {
    printf("========= Failed SBM tests 6 for SBM  ========= \n");
    return 1;
  }

  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);

}
