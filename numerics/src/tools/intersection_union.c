#include "intersection_union.h"
#include<stdio.h>
#include<stdlib.h>

void int_array_print(int * idx, int idx_size)
{
  printf("[");
  for (int i = 0; i <  idx_size; i++  ) printf("%i ", idx[i]);
  printf("]\n");
}

void int_array_of_array_print(int ** idx, int idx_size, int array_size)
{
  printf("[");
  for (int i = 0; i <  idx_size; i++  )
  {
    printf("[");
    for (int j = 0; j < array_size  ; j++  ) printf("%i ", idx[i][j]);
    printf("] ");
  } 
  printf("]\n");
}

void uint_array_print(unsigned int * idx, unsigned int idx_size)
{
  printf("[");
  for (unsigned int i = 0; i <  idx_size; i++  ) printf("%i ", idx[i]);
  printf("]\n");
}
void uint_array_of_array_print(unsigned int ** idx, unsigned int idx_size, unsigned int array_size)
{
  printf("[");
  for (unsigned int i = 0; i <  idx_size; i++  )
  {
    printf("[");
    for (unsigned int j = 0; j < array_size  ; j++  ) printf("%i ", idx[i][j]);
    printf("] ");
  }
  printf("]\n");
}
int int_array_check_sorted(int * idx, int idx_size)
{
  for (int i = 1; i <  idx_size; i++  )
    if (idx[i] < idx[i-1])
      return 1;
  return 0;
}

int int_array_of_array_check_sorted(int ** idx, int idx_size, int pos)
{
  for (int i = 1; i <  idx_size; i++  )
  {
    if (idx[i][pos] < idx[i-1][pos])
      return 1;
  }
  return 0;
}
int uint_array_check_sorted(unsigned int * idx, unsigned int idx_size)
{
  for (unsigned int i = 1; i <  idx_size; i++  )
    if (idx[i] < idx[i-1])
      return 1;
  return 0;
}

int uint_array_of_array_check_sorted(unsigned int ** idx, unsigned int idx_size, unsigned int pos)
{
  for (unsigned int i = 1; i <  idx_size; i++  )
  {
    if (idx[i][pos] < idx[i-1][pos])
      return 1;
  }
  return 0;
}
/* Function prints Intersection of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
void sorted_array_intersection_print(int arr1[], int arr2[], int m, int n)
{
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      printf(" %d ", arr2[j++]);
      i++;
    }
  }
}

void sorted_array_of_array_intersection_with_array(int * arr1, int ** arr2, int m, int n, int array_size, int pos,  int ** intersection_set, int * intersection_set_size )
{
  int i = 0, j = 0;
  int size= 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j][pos])
      i++;
    else if (arr2[j][pos] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      intersection_set[size][0] = arr2[j++][0];
      for (int k=1;  k< array_size ; k++)
        intersection_set[size][k] = arr2[j-1][k];
      size++;
      i++;
    }
  }
  *intersection_set_size = size;
}

void sorted_array_intersection(int * arr1, int * arr2, int m, int n, int * intersection_set, int * intersection_set_size )
{
  
  int i = 0, j = 0;
  int size= 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      /* printf(" %d ", arr2[j++]); */
      intersection_set[size] = arr2[j++];
      size++;
      i++;
    }
  }
  *intersection_set_size = size;
}



/* Function prints union of arr1[] and arr2[]
   m is the number of elements in arr1[]
   n is the number of elements in arr2[] */
void printUnion(int arr1[], int arr2[], int m, int n)
{
  int i = 0, j = 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      printf(" %d ", arr1[i++]);
    else if (arr2[j] < arr1[i])
      printf(" %d ", arr2[j++]);
    else
    {
      printf(" %d ", arr2[j++]);
      i++;
    }
  }
 
  /* Print remaining elements of the larger array */
  while(i < m)
   printf(" %d ", arr1[i++]);
  while(j < n)
   printf(" %d ", arr2[j++]);
}

void array_union(int arr1[], int arr2[], int m, int n, int * union_set, int * union_set_size)
{
  int i = 0, j = 0;
  int size =0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
    {
      //printf(" %d ", arr1[i++]);
      union_set[size]=arr1[i++];
      size++;
    }
    else if (arr2[j] < arr1[i])
    {
      //printf(" %d ", arr2[j++]);
      union_set[size]=arr2[j++];
      size++;
    }
    else
    {
      //printf(" %d ", arr2[j++]);
      union_set[size]=arr2[j++];
      size++;
      i++;
    }
  }
 
  /* Print remaining elements of the larger array */
  while(i < m)
    //printf(" %d ", arr1[i++]);
  {
    union_set[size] = arr1[i++];
    size++;
  }
  while(j < n)
  {
    union_set[size] = arr2[j++];
    size++;
  }
  * union_set_size = size;
  //printf(" %d ", arr2[j++]);
  
}



void uint_sorted_array_of_array_intersection_with_array(unsigned int * arr1, unsigned int ** arr2, unsigned int m, unsigned int n,
                                                 unsigned int array_size, unsigned int pos,  unsigned int ** intersection_set, unsigned int * intersection_set_size )
{
  unsigned int i = 0, j = 0;
  unsigned int size= 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j][pos])
      i++;
    else if (arr2[j][pos] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      intersection_set[size][0] = arr2[j++][0];
      for (unsigned int k=1;  k< array_size ; k++)
        intersection_set[size][k] = arr2[j-1][k];
      size++;
      i++;
    }
  }
  *intersection_set_size = size;
}


void uint_sorted_array_intersection(unsigned int * arr1, unsigned int * arr2, unsigned int m, unsigned int n,
                             unsigned int * intersection_set, unsigned int * intersection_set_size )
{

  unsigned int i = 0, j = 0;
  unsigned int size= 0;
  while (i < m && j < n)
  {
    if (arr1[i] < arr2[j])
      i++;
    else if (arr2[j] < arr1[i])
      j++;
    else /* if arr1[i] == arr2[j] */
    {
      /* printf(" %d ", arr2[j++]); */
      intersection_set[size] = arr2[j++];
      size++;
      i++;
    }
  }
  *intersection_set_size = size;
}






/* /\* Driver program to test above function *\/ */
/* int main() */
/* { */
/*   int arr1[] = {1, 2, 4, 5, 6}; */
/*   int arr2[] = {2, 3, 5, 7}; */
/*   int m = sizeof(arr1)/sizeof(arr1[0]); */
/*   int n = sizeof(arr2)/sizeof(arr2[0]); */
/*   printUnion(arr1, arr2, m, n); */


/*   int max_size=m+n ; */
/*   printf("\n"); */
/*   int * union_set = (int*) malloc(max_size*sizeof(int)); */
/*   int union_set_size; */
/*   compute_union(arr1, arr2, m,  n,  union_set, &union_set_size ); */
/*   printf("union_set_size = %i\n", union_set_size); */
/*   for (int i=0; i < union_set_size; i++) printf("% i ", union_set[i]); */

  
/*   getchar(); */
/*   return 0; */
/* } */


/* /\* Driver program to test above function *\/ */
/* int main() */
/* { */
/*   int arr1[] = {1, 2, 4, 5, 6}; */
/*   int arr2[] = {2, 3, 5, 7}; */
/*   int m = sizeof(arr1)/sizeof(arr1[0]); */
/*   int n = sizeof(arr2)/sizeof(arr2[0]); */
/*   sorted_array_intersection_print(arr1, arr2, m, n); */
/*   int max_size=0; */
/*   if (m<n) */
/*   { */
/*     max_size = m; */
/*   } */
/*   else */
/*     max_size = n; */

/*   printf("\n"); */
/*   int * intersectionset = (int*) malloc(max_size*sizeof(int)); */
/*   int intersection_size; */
/*   intersection(arr1, arr2, m,  n,  intersectionset, &intersection_size ); */
/*   printf("intersection_size = %i\n", intersection_size); */
/*   for (int i=0; i < intersection_size; i++) printf("% i ", intersectionset[i]); */
/*   getchar(); */
/*   return 0; */
/* } */
