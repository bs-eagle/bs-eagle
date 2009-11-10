#ifndef __MEMORY_MACROSES_H_
#define __MEMORY_MACROSES_H_
/*!
 * \file memory_macroses.h
 * \brief usefull macroses for use in solver
 * \author Borschuk Oleg
 * \date 2006-06-19
 */

//! Allocate double array DEST of length LEN elements and copy SRC to DEST
#define MALLOC_COPY_DOUBLE_ARRAY(DEST,SRC,LEN)      \
  MEMORY_ALLOCATION_TEST ((DEST) = new double[(LEN)]); \
  COPY_DOUBLE_ARRAY ((DEST), (SRC), (LEN))

//! allocate memory for double array A -- array, N -- length, F -- allocator error if != 0 error occur
#ifdef MEM_DEBUG
#define FI_DOUBLE_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new double[(N)];                                                \
  if (!(A)) (F) = 1;                                                    \
  else FI_DOUBLE_ARRAY_SET_ILLEGAL(A,N);
#else // MEM_DEBUG
#define FI_DOUBLE_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new double[(N)];                                                \
  if (!(A)) (F) = 1;
#endif // MEM_DEBUG

//! allocate memory for float array A -- array, N -- length, F -- allocator error if != 0 error occur
#ifdef MEM_DEBUG
#define FI_FLOAT_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new float[(N)];                                                \
  if (!(A)) (F) = 1;                                                    \
  else FI_FLOAT_ARRAY_SET_ILLEGAL(A,N);
#else // MEM_DEBUG
#define FI_FLOAT_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new float[(N)];                                                \
  if (!(A)) (F) = 1;
#endif // MEM_DEBUG

//! allocate memory for T array A -- array, N -- length, F -- allocator error if != 0 error occur
#ifdef MEM_DEBUG
#define FI_T_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new T[(N)];                                                \
  if (!(A)) (F) = 1;                                                    \
  else FI_T_ARRAY_SET_ILLEGAL(A,N);
#else // MEM_DEBUG
#define FI_T_ARRAY_ALLOCATOR(A,N,F)                                \
  (A) = new T[(N)];                                                \
  if (!(A)) (F) = 1;
#endif // MEM_DEBUG

//! Copy double array SRC of length LEN elements to double array DEST of length LEN elements
#define COPY_DOUBLE_ARRAY(DEST,SRC,LEN) \
  memcpy ((DEST), (SRC), (LEN) * sizeof (double))

//! Copy float array SRC of length LEN elements to float array DEST of length LEN elements
#define COPY_FLOAT_ARRAY(DEST,SRC,LEN) \
  memcpy ((DEST), (SRC), (LEN) * sizeof (float))



//! allocate memory for int array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_INT_ARRAY_ALLOCATOR(A,N,F)                                   \
  (A) = new int[(N)];                                                   \
  if (!(A)) (F) = 1;

//! allocate memory for unsigned int array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_UINT_ARRAY_ALLOCATOR(A,N,F)                                  \
  (A) = new unsigned int[(N)];                                          \
  if (!(A)) (F) = 1;


//! allocate memory for char array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_CHAR_ARRAY_ALLOCATOR(A,N,F)                                  \
  (A) = new char[(N)];                                                  \
  if (!(A)) (F) = 1;


//! reallocate memory for double array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_DOUBLE_ARRAY_REALLOCATOR(A,N,F)                              \
  if (A) delete[] (A);                                                  \
  FI_DOUBLE_ARRAY_ALLOCATOR(A,N,F);

//! reallocate memory for float array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_FLOAT_ARRAY_REALLOCATOR(A,N,F)                              \
  if (A) delete[] (A);                                                  \
  FI_FLOAT_ARRAY_ALLOCATOR(A,N,F);

//! reallocate memory for int array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_INT_ARRAY_REALLOCATOR(A,N,F)                                 \
  if (A) delete[] (A);                                                  \
  FI_INT_ARRAY_ALLOCATOR(A,N,F);

//! reallocate memory for T array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_T_ARRAY_REALLOCATOR(A,N,F)                                 \
  if (A) delete[] (A);                                                  \
  FI_T_ARRAY_ALLOCATOR(A,N,F);

//! reallocate memory for unsigned int array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_UINT_ARRAY_REALLOCATOR(A,N,F)                                 \
  if (A) delete[] (A);                                                  \
  FI_UINT_ARRAY_ALLOCATOR(A,N,F);


//! reallocate memory for char array A -- array, N -- length, F -- allocator error if != 0 error occur
#define FI_CHAR_ARRAY_REALLOCATOR(A,N,F)                                \
  if (A) delete[] (A);                                                  \
  FI_CHAR_ARRAY_ALLOCATOR(A,N,F);

//! set values of double array to illegal A -- array, N -- length
#define FI_DOUBLE_ARRAY_SET_ILLEGAL(A,N)                                \
  memset ((A), -1, (N) * sizeof (double));

//! set values of float array to illegal A -- array, N -- length
#define FI_FLOAT_ARRAY_SET_ILLEGAL(A,N)                                \
  memset ((A), -1, (N) * sizeof (float));

//! set values of float array to illegal A -- array, N -- length
#define FI_T_ARRAY_SET_ILLEGAL(A,N)                                \
  memset ((A), 0, (N) * sizeof (T));
//! reallocate array and copy all data from previos
#define FI_DOUBLE_ARRAY_COPYREALLOC(A,NOLD,NNEW,F)                      \
  if (!(A))                                                             \
    {                                                                   \
      FI_DOUBLE_ARRAY_ALLOCATOR((A),(NNEW),(F));                        \
      (NOLD) = (NNEW);                                                  \
    }                                                                   \
  else if ((NNEW) > (NOLD))                                             \
    {                                                                   \
      int n_true = (NOLD);                                              \
      double *d_temp = 0;                                               \
      while (n_true < (NNEW))                                           \
        n_true *= 2;                                                    \
      FI_DOUBLE_ARRAY_ALLOCATOR (d_temp, n_true, (F));                  \
      if (d_temp)                                                       \
        {                                                               \
          memcpy (d_temp, (A), sizeof (double) * (NOLD));               \
          delete [] (A);                                                \
          (A) = d_temp;                                                 \
        }                                                               \
      (NOLD) = n_true;                                                  \
    }

//! reallocate array and copy all data from previos
#define FI_DOUBLE_ARRAY_REALLOC_COPY(A,NOLD,NNEW,F)                     \
  {                                                                     \
      double *d_temp = 0;                                               \
      FI_DOUBLE_ARRAY_ALLOCATOR (d_temp, (NNEW), (F));                  \
      if (d_temp)                                                       \
        {                                                               \
          if ((NNEW) >= (NOLD))                                         \
            {                                                           \
              memcpy (d_temp, (A), sizeof (double) * (NOLD));           \
            }                                                           \
          else                                                          \
            {                                                           \
              memcpy (d_temp, (A), sizeof (double) * (NNEW));           \
            }                                                           \
          delete [] (A);                                                \
          (A) = d_temp;                                                 \
        }                                                               \
   }

//! reallocate array and copy all data from previos
#define FI_INT_ARRAY_COPYREALLOC(A,NOLD,NNEW,F)                         \
  if (!(A))                                                             \
    {                                                                   \
      FI_INT_ARRAY_ALLOCATOR((A),(NNEW),(F));                           \
      (NOLD) = (NNEW);                                                  \
    }                                                                   \
  else if ((NNEW) > (NOLD))                                             \
    {                                                                   \
      int n_true = (NOLD);                                              \
      int *i_temp;                                                      \
      while (n_true < (NNEW))                                           \
        n_true *= 2;                                                    \
      FI_INT_ARRAY_ALLOCATOR (i_temp, n_true, (F));                     \
      if (i_temp)                                                       \
        {                                                               \
          memcpy (i_temp, (A), sizeof (int) * (NOLD));                  \
          delete [] (A);                                                \
          (A) = i_temp;                                                 \
        }                                                               \
      (NOLD) = n_true;                                                  \
    }

//! free array
#define FI_FREE_ARRAY(A)                                                \
  if (A) delete[] (A);                                                  \
  (A) = 0;

//! free
#define FI_FREE(A)                                                      \
  if (A) delete A;                                                      \
  (A) = 0;

//! fill array A -- array, N1 -- start, N2 -- end, V -- value
#define FI_FILL_ARRAY(A,N1,N2,V)                                        \
  for (int iii = (int)(N1); iii < (int)(N2); ++iii)                     \
    {                                                                   \
      (A)[iii] = (V);                                                   \
    }
//! decrement array by value A -- array, N1 -- start, N2 -- end, V -- value
#define FI_DECR_ARRAY(A,N1,N2,V)                                        \
  for (int iii = (int)(N1); iii < (int)(N2); ++iii)                     \
    {                                                                   \
      (A)[iii] -= (V);                                                  \
    }

#endif //__MEMORY_MACROSES_H_
