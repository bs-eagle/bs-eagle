#ifndef __MATRIX_MACROSES_H_
#define __MATRIX_MACROSES_H_
/*!
 * \file matrix_macroses.h
 * \brief usefull macroses for use in solver
 * \author Borschuk Oleg
 * \date 2006-06-19
 */

//! compute matrix x vector product
inline void
mv_product (int n, double *m, double *v, double *r)
{
  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          r[i] += m[j + i * n] * v[j];
        }
    }
}

//! find solution for lower triangular L matrix with 1 on diagonal
// in v, as i understood, we storage vector of rhs
inline void
ulu_find_root_l (int n, double *l, double *v)
{
  int i, j;
  for (i = 1; i < n; ++i)
    {
      for (j = 0; j < i; ++j)
        {
          v[i] -= l[j + i * n] * v[j];
        }
    }
}

//! find solution for upper triangular U matrix
inline void
ulu_find_root_u (int n, double *u, double *v)
{
  int i, j;
  for (i = n - 1; i >= 0; --i)
    {
      for (j = n - 1; j > i; --j)
        {
          v[i] -= u[j + i * n] * v[j];
        }
      v[i] /= u[i + i * n];
    }
}
//! update for solution
inline void
lu_find_root_update (int n, double *m, double *v, double *r)
{
  int i, j;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          r[i] -= m[j + i * n] * v[j];
        }
    }
}

//! calculate B = B - L*U
inline void
lu_upgrade (int n, double *l, double *u, double *b)
{
  int i, j, k;
  double v;
  for (i = 0; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          v = l[j + i * n];
          for (k = 0; k < n; ++k)
            {
              b[k + i * n] -= v * u[k + j * n];
            }
        }
    }
}

//! matrix to matrix product
inline void
mm_prod (int n, double *m, double *r)
{
  int  i;
  int b_sqr = n * n;
  for (i = 0; i < b_sqr; ++i)
    {
      r[i] += m[i];
    }
}

//! find L matrix in LU factorization (solve LU = A for unknown L)
inline void
ulu_seek_l (int n, double *u, double *l, double /* t */)
{
  int i, j, k;
  for (j = 0; j < n; ++j)
    {
      for (i = 0; i < n; ++i)
        {
          l[j + i * n] /= u[j + j * n];
          for (k = j + 1; k < n; ++k)
            {
              l[k + i * n] -= l[j + i * n] * u[k + j * n];
            }
        }
    }
}

//! LU factorization for n x n block
inline void
ulu (int n, double *m, double &/* t */)
{
  int i, j, k;
  for (i = 1; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          for (k = 0; (k < i) && (k < j); ++k)
            {
              m[j + i * n] -= m[k + i * n] * m[j + k * n];
            }
          if (j < i)
            m[j + i * n] /= m[j + j * n];
        }
    }
}

//! find U matrix in LU factorization (solve LU = A for unknown U)
inline void
ulu_seek_u (int n, double *l, double *u)
{
  int i, j, k;
  for (i = 1; i < n; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          for (k = 0; k < i; ++k)
            {
              u[j + i * n] -= l[k + i * n] * u[j + k * n];
            }
        }
    }
}

//!
inline void
vv_prod_m (int n, double *v1, double *v2, double &r)
{
  int i;
  for (i = 0; i < n; ++i)
    {
      r -= v1[i] * v2[i];
    }
}

template <class m_vector_t, class v_vector_t, class r_vector_t>
void
mv_prod_1x1 (const m_vector_t &m, const v_vector_t &v, r_vector_t &r)
{
  typedef typename r_vector_t::value_type r_item_t;
  r[0] += r_item_t(m[0]) * r_item_t(v[0]);
}

template <class m_vector_t, class v_vector_t, class r_vector_t>
void
mv_prod_2x2 (const m_vector_t &m, const v_vector_t &v, r_vector_t &r)
{
  typedef typename r_vector_t::value_type r_item_t;
  r[0] += r_item_t(m[0]) * r_item_t(v[0]) + r_item_t(m[1]) * r_item_t(v[1]);
  r[1] += r_item_t(m[2]) * r_item_t(v[0]) + r_item_t(m[3]) * r_item_t(v[1]);
}

template <class m_vector_t, class v_vector_t, class r_vector_t>
void
mv_prod_3x3 (const m_vector_t &m, const v_vector_t &v, r_vector_t &r)
{
  typedef typename r_vector_t::value_type r_item_t;
  r[0] += r_item_t(m[0]) * r_item_t (v[0]) + r_item_t(m[1]) * r_item_t(v[1]) + r_item_t(m[2]) * r_item_t(v[2]);
  r[1] += r_item_t(m[3]) * r_item_t(v[0]) + r_item_t(m[4]) * r_item_t(v[1]) + r_item_t(m[5]) * r_item_t(v[2]);
  r[2] += r_item_t(m[6]) * r_item_t(v[0]) + r_item_t(m[7]) * r_item_t(v[1]) + r_item_t(m[8]) * r_item_t(v[2]);
}


#ifndef COMPOSITIONAL
//! matrix vector product for 1x1 matrix
#define MV_PROD_1x1(M,V,R)   ((R)[0] += (M)[0] * (V)[0]);
//! matrix vector product for 2x2 matrix
#define MV_PROD_2x2(M,V,R)   ((R)[0] += (M)[0] * (V)[0] + (M)[1] * (V)[1]); ((R)[1] += (M)[2] * (V)[0] + (M)[3] * (V)[1]);
//! matrix vector product for 3x3 matrix
#define MV_PROD_3x3(M,V,R)      \
  ((R)[0] += (M)[0] * (V)[0] + (M)[1] * (V)[1] + (M)[2] * (V)[2]);              \
  ((R)[1] += (M)[3] * (V)[0] + (M)[4] * (V)[1] + (M)[5] * (V)[2]);              \
  ((R)[2] += (M)[6] * (V)[0] + (M)[7] * (V)[1] + (M)[8] * (V)[2]);
//! matrix vector product for 1x1, 2x2 and 3x3 matrixies
#define MV_PROD(F,M,V,R)                          \
  if ((F) == 1)                                   \
    {MV_PROD_1x1 ((M), (V), (R));}                \
  else if ((F) == 2)                              \
    {MV_PROD_2x2 ((M), (V), (R));}                \
  else if ((F) == 3)                              \
    {MV_PROD_3x3 ((M), (V), (R));}
#else
//! matrix vector product for 1x1, 2x2 and 3x3 matrixies
#define MV_PROD(F,M,V,R)  \
   mv_product ((F), (M), (V), (R));
#endif


//! matrix vector product for 2x2 matrix
#define MV_PROD_MULT_2x2(M,V,MULT,R)   \
  ((R)[0] += (MULT) * ((M)[0] * (V)[0] + (M)[1] * (V)[1])); \
  ((R)[1] += (MULT) * ((M)[2] * (V)[0] + (M)[3] * (V)[1]));
//! matrix vector product for 3x3 matrix
#define MV_PROD_MULT_3x3(M,V,MULT,R)      \
  ((R)[0] += (MULT) * ((M)[0] * (V)[0] + (M)[1] * (V)[1] + (M)[2] * (V)[2]));              \
  ((R)[1] += (MULT) * ((M)[3] * (V)[0] + (M)[4] * (V)[1] + (M)[5] * (V)[2]));              \
  ((R)[2] += (MULT) * ((M)[6] * (V)[0] + (M)[7] * (V)[1] + (M)[8] * (V)[2]));
//! matrix vector product for 1x1, 2x2 and 3x3 matrixies
#define MV_PROD_MULT(F,M,V,MULT,R)                \
  if ((F) == 1)                                   \
    {(R)[0] += (MULT) * ((M)[0] * (V)[0]);}       \
  else if ((F) == 2)                              \
    {MV_PROD_MULT_2x2 ((M), (V), (MULT), (R));}   \
  else if ((F) == 3)                              \
    {MV_PROD_MULT_3x3 ((M), (V), (MULT), (R));}

//! matrix vector product minus for 1x1, 2x2 and 3x3 matrixies
#define MV_PROD_MINUS_2x2(M,V,R)                                                \
  ((R)[0] -= (M)[0] * (V)[0] + (M)[1] * (V)[1]);                                \
  ((R)[1] -= (M)[2] * (V)[0] + (M)[3] * (V)[1]);

#define MV_PROD_MINUS_3x3(M,V,R)                                                \
  ((R)[0] -= (M)[0] * (V)[0] + (M)[1] * (V)[1] + (M)[2] * (V)[2]);              \
  ((R)[1] -= (M)[3] * (V)[0] + (M)[4] * (V)[1] + (M)[5] * (V)[2]);              \
  ((R)[2] -= (M)[6] * (V)[0] + (M)[7] * (V)[1] + (M)[8] * (V)[2]);


#define MV_PROD_MINUS(F,M,V,R)                          \
  if ((F) == 1)                                         \
    {(R)[0] -= (M)[0] * (V)[0];}                        \
  else if ((F) == 2)                                    \
    {MV_PROD_MINUS_2x2 ((M), (V), (R));}                \
  else if ((F) == 3)                                    \
    {MV_PROD_MINUS_3x3 ((M), (V), (R));}


#ifndef COMPOSITIONAL
//! LU factorization for 2x2 block
#define uLU_2x2(M)                                      \
  (M)[2] /= (M)[0];                                     \
  (M)[3] -= (M)[2] * (M)[1];

//! LU factorization for 3x3 block
#define uLU_3x3(M,T)                                    \
  (T) = 1.0 / (M)[0];                                   \
  (M)[3] *= (T);                                        \
  (M)[4] -= (M)[3] * (M)[1];                            \
  (M)[5] -= (M)[3] * (M)[2];                            \
  (M)[6] *= (T);                                        \
  (M)[7] -= (M)[6] * (M)[1];                            \
  (M)[7] /= (M)[4];                                     \
  (M)[8] -= (M)[6] * (M)[2] + (M)[7] * (M)[5];

//! LU factorization for 1x1, 2x2 and 3x3 blocks
#define uLU(F,M,T)                                      \
  if ((F) == 2)                                         \
    {uLU_2x2 ((M));}                                    \
  else if ((F) == 3)                                    \
    {uLU_3x3 ((M), (T));}
#else // COMPOSITIONAL
#define uLU(F,M,T)                                      \
  ulu ((F), (M), (T));
#endif // COMPOSITIONAL

#ifndef COMPOSITIONAL
#define uLU_SEEK_U_2x2(L,U)                             \
  (U)[2] -= (U)[0] * (L)[2];                            \
  (U)[3] -= (U)[1] * (L)[2];

#define uLU_SEEK_U_3x3(L,U)                             \
  (U)[3] -= (U)[0] * (L)[3];                            \
  (U)[4] -= (U)[1] * (L)[3];                            \
  (U)[5] -= (U)[2] * (L)[3];                            \
  (U)[6] -= (U)[0] * (L)[6] + (U)[3] * (L)[7];          \
  (U)[7] -= (U)[1] * (L)[6] + (U)[4] * (L)[7];          \
  (U)[8] -= (U)[2] * (L)[6] + (U)[5] * (L)[7];

#define uLU_SEEK_U(F,L,U)                               \
  if ((F) == 2)                                         \
    {uLU_SEEK_U_2x2 ((L), (U));}                        \
  else if ((F) == 3)                                    \
    {uLU_SEEK_U_3x3 ((L), (U));}
#else // COMPOSITIONAL
#define uLU_SEEK_U(F,L,U)                               \
  ulu_seek_u ((F), (L), (U));
#endif // COMPOSITIONAL

#ifndef COMPOSITIONAL
#define uLU_SEEK_L_2x2(U,L,T)                           \
  (T) = 1.0 / (U)[0];                                   \
  (L)[0] *= (T);                                        \
  (L)[2] *= (T);                                        \
  (T) = 1.0 / (U)[3];                                   \
  (L)[1] = ((L)[1] - (L)[0] * (U)[1]) * (T);            \
  (L)[3] = ((L)[3] - (L)[2] * (U)[1]) * (T);

#define uLU_SEEK_L_3x3(U,L,T)                                           \
  (T) = 1.0 / (U)[0];                                                   \
  (L)[0] *= (T);                                                        \
  (L)[3] *= (T);                                                        \
  (L)[6] *= (T);                                                        \
  (T) = 1.0 / (U)[4];                                                   \
  (L)[1] = ((L)[1] - (L)[0] * (U)[1]) * (T);                            \
  (L)[4] = ((L)[4] - (L)[3] * (U)[1]) * (T);                            \
  (L)[7] = ((L)[7] - (L)[6] * (U)[1]) * (T);                            \
  (T) = 1.0 / (U)[8];                                                   \
  (L)[2] = ((L)[2] - (L)[0] * (U)[2] - (L)[1] * (U)[5]) * (T);          \
  (L)[5] = ((L)[5] - (L)[3] * (U)[2] - (L)[4] * (U)[5]) * (T);          \
  (L)[8] = ((L)[8] - (L)[6] * (U)[2] - (L)[7] * (U)[5]) * (T);

#define uLU_SEEK_L(F,U,L,T)                             \
  if ((F) == 1)                                         \
    (L)[0] /= (U)[0];                                   \
  if ((F) == 2)                                         \
    {uLU_SEEK_L_2x2 ((U), (L), (T));}                   \
  else if ((F) == 3)                                    \
    {uLU_SEEK_L_3x3 ((U), (L), (T));}
#else // COMPOSITIONAL
#define uLU_SEEK_L(F,U,L,T)                             \
  ulu_seek_l ((F), (U), (L), (T));
#endif // COMPOSITIONAL

#ifndef COMPOSITIONAL
#define LU_UPGRADE_2x2(L,U,B)                           \
  (B[0]) -= (L)[0] * (U)[0] + (L)[1] * (U)[2];          \
  (B[1]) -= (L)[0] * (U)[1] + (L)[1] * (U)[3];          \
  (B[2]) -= (L)[2] * (U)[0] + (L)[3] * (U)[2];          \
  (B[3]) -= (L)[2] * (U)[1] + (L)[3] * (U)[3];

#define LU_UPGRADE_3x3(L,U,B)                                           \
  (B)[0] -= (L)[0] * (U)[0] + (L)[1] * (U)[3] + (L)[2] * (U)[6];        \
  (B)[1] -= (L)[0] * (U)[1] + (L)[1] * (U)[4] + (L)[2] * (U)[7];        \
  (B)[2] -= (L)[0] * (U)[2] + (L)[1] * (U)[5] + (L)[2] * (U)[8];        \
  (B)[3] -= (L)[3] * (U)[0] + (L)[4] * (U)[3] + (L)[5] * (U)[6];        \
  (B)[4] -= (L)[3] * (U)[1] + (L)[4] * (U)[4] + (L)[5] * (U)[7];        \
  (B)[5] -= (L)[3] * (U)[2] + (L)[4] * (U)[5] + (L)[5] * (U)[8];        \
  (B)[6] -= (L)[6] * (U)[0] + (L)[7] * (U)[3] + (L)[8] * (U)[6];        \
  (B)[7] -= (L)[6] * (U)[1] + (L)[7] * (U)[4] + (L)[8] * (U)[7];        \
  (B)[8] -= (L)[6] * (U)[2] + (L)[7] * (U)[5] + (L)[8] * (U)[8];

#define LU_UPGRADE(F,L,U,B)                             \
  if ((F) == 1)                                         \
    {(B)[0] -= (L)[0] * (U)[0];}                        \
  else if ((F) == 2)                                    \
    {LU_UPGRADE_2x2 ((L), (U), (B));}                   \
  else if ((F) == 3)                                    \
    {LU_UPGRADE_3x3 ((L), (U), (B));}
#else // COMPOSITIONAL
#define LU_UPGRADE(F,L,U,B)                             \
  lu_upgrade ((F), (L), (U), (B));
#endif // COMPOSITIONAL

#ifndef COMPOSITIONAL
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
#define LU_FIND_ROOT_UPDATE_2x2(M,V,R)                  \
  (R)[0] -= (V)[0] * (M)[0] + (V)[1] * (M)[1];          \
  (R)[1] -= (V)[0] * (M)[2] + (V)[1] * (M)[3];

#define LU_FIND_ROOT_UPDATE_3x3(M,V,R)                  \
  (R)[0] -= (V)[0] * (M)[0] + (V)[1] * (M)[1] + (V)[2] * (M)[2];        \
  (R)[1] -= (V)[0] * (M)[3] + (V)[1] * (M)[4] + (V)[2] * (M)[5];        \
  (R)[2] -= (V)[0] * (M)[6] + (V)[1] * (M)[7] + (V)[2] * (M)[8];

#define LU_FIND_ROOT_UPDATE(F,M,V,R)                                    \
  if ((F) == 1)                                                         \
    {(R)[0] -= (V)[0] * (M)[0];}                                        \
  else if ((F) == 2)                                                    \
    {LU_FIND_ROOT_UPDATE_2x2 ((M), (V), (R));}                          \
  else if ((F) == 3)                                                    \
    {LU_FIND_ROOT_UPDATE_3x3 ((M), (V), (R));}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
#else // COMPOSITIONAL
#define LU_FIND_ROOT_UPDATE(F,M,V,R)                                    \
  lu_find_root_update ((F), (M), (V), (R));
#endif // COMPOSITIONAL

#ifndef COMPOSITIONAL
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
#define uLU_FIND_ROOT_L_2x2(L,V)                        \
  (V)[1] -= (V)[0] * (L)[2];

#define uLU_FIND_ROOT_L_3x3(L,V)                        \
  (V)[1] -= (V)[0] * (L)[3];                            \
  (V)[2] -= ((V)[0] * (L)[6] + (V)[1] * (L)[7]);

#define uLU_FIND_ROOT_L(F,L,V)                          \
  if ((F) == 2)                                         \
    {uLU_FIND_ROOT_L_2x2 ((L), (V));}                   \
  else if ((F) == 3)                                    \
    {uLU_FIND_ROOT_L_3x3 ((L), (V));}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
#else
#define uLU_FIND_ROOT_L(F,L,V)                          \
  ulu_find_root_l ((F), (L), (V));
#endif

#ifndef COMPOSITIONAL
#define uLU_FIND_ROOT_U_2x2(U,V)                        \
  (V)[1] /= (U)[3];                                     \
  (V)[0] = ((V)[0] - (V)[1] * (U)[1]) / (U)[0];

#define uLU_FIND_ROOT_U_3x3(U,V)                                        \
  (V)[2] /= (U)[8];                                                     \
  (V)[1] = ((V)[1] - (V)[2] * (U)[5]) / (U)[4];                         \
  (V)[0] = ((V)[0] - (V)[2] * (U)[2] - (V)[1] * (U)[1]) / (U)[0];

#define uLU_FIND_ROOT_U(F,U,V)                          \
  if ((F) == 1)                                         \
    {(V)[0] /= (U)[0];}                                 \
  else if ((F) == 2)                                    \
    {uLU_FIND_ROOT_U_2x2 ((U), (V));}                   \
  else if ((F) == 3)                                    \
    {uLU_FIND_ROOT_U_3x3 ((U), (V));}
#else
#define uLU_FIND_ROOT_U(F,U,V)                          \
  ulu_find_root_u ((F), (U), (V));
#endif

#define uLU_FIND_MATRIX_ROOT_L_2x2(L,M)                 \
  (M)[2] -= (M)[0] * (L)[2];                            \
  (M)[3] -= (M)[1] * (L)[2];

#define uLU_FIND_MATRIX_ROOT_L_3x3(L,M)                                 \
  (M)[3] -= (M)[0] * (L)[3];                                            \
  (M)[6] -= ((M)[0] * (L)[6] + (M)[3] * (L)[7]);                        \
  (M)[4] -= (M)[1] * (L)[3];                                            \
  (M)[7] -= ((M)[1] * (L)[6] + (M)[4] * (L)[7]);                        \
  (M)[5] -= (M)[2] * (L)[3];                                            \
  (M)[8] -= ((M)[2] * (L)[6] + (M)[5] * (L)[7]);

#define uLU_FIND_MATRIX_ROOT_L(F,L,M)                                   \
  if ((F) == 2)                                                         \
    {uLU_FIND_MATRIX_ROOT_L_2x2 ((L), (M));}                            \
  else if ((F) == 3)                                                    \
    {uLU_FIND_MATRIX_ROOT_L_3x3 ((L), (M));}


#define uLU_FIND_MATRIX_ROOT_U_2x2(U,M)                 \
  (M)[2] /= (U)[3];                                     \
  (M)[0] = ((M)[0] - (M)[2] * (U)[1]) / (U)[0];         \
  (M)[3] /= (U)[3];                                     \
  (M)[1] = ((M)[1] - (M)[3] * (U)[1]) / (U)[0];

#define uLU_FIND_MATRIX_ROOT_U_3x3(U,M)                                         \
  (M)[6] /= (U)[8];                                                             \
  (M)[3] = ((M)[3] - (M)[6] * (U)[5]) / (U)[4];                                 \
  (M)[0] = ((M)[0] - (M)[6] * (U)[2] - (M)[3] * (U)[1]) / (U)[0];               \
  (M)[7] /= (U)[8];                                                             \
  (M)[4] = ((M)[4] - (M)[7] * (U)[5]) / (U)[4];                                 \
  (M)[1] = ((M)[1] - (M)[7] * (U)[2] - (M)[4] * (U)[1]) / (U)[0];               \
  (M)[8] /= (U)[8];                                                             \
  (M)[5] = ((M)[5] - (M)[8] * (U)[5]) / (U)[4];                                 \
  (M)[2] = ((M)[2] - (M)[8] * (U)[2] - (M)[5] * (U)[1]) / (U)[0];

#define uLU_FIND_MATRIX_ROOT_U(F,U,M)                                   \
  if ((F) == 1)                                                         \
    (M)[0] /= (U)[0];                                                   \
  else if ((F) == 2)                                                    \
    {uLU_FIND_MATRIX_ROOT_U_2x2 ((U), (M));}                            \
  else if ((F) == 3)                                                    \
    {uLU_FIND_MATRIX_ROOT_U_3x3 ((U), (M));}

#ifndef COMPOSITIONAL
#define VV_PROD_M_2(V1,V2,R)  (R) -= (V1)[0] * (V2)[0] + (V1)[1] * (V2)[1];
#define VV_PROD_M_3(V1,V2,R)  (R) -= (V1)[0] * (V2)[0] + (V1)[1] * (V2)[1] + (V1)[2] * (V2)[2];
#define VV_PROD_M(F,V1,V2,R)                                            \
  if ((F) == 1)                                                         \
    (R) -= (V1)[0] * (V2)[0];                                           \
  else if ((F) == 2)                                                    \
    {VV_PROD_M_2 (V1, V2, R);}                                          \
  else if ((F) == 3)                                                    \
    {VV_PROD_M_3 (V1, V2, R);}
#else // COMPOSITIONAL
#define VV_PROD_M(F,V1,V2,R)                          \
  vv_prod_m ((F), (V1), (V2), (R));
#endif // COMPOSITIONAL


#define INVERSE_MATRIX_2x2(M,DET,D)                     \
  (DET) = ((M)[0] * (M)[3] - (M)[1] * (M)[2]);         \
  (D) = (M)[0];                                         \
  (M)[0] = (M)[3] / (DET);                              \
  (M)[3] = (D) / (DET);                                 \
  (D) = (M)[1];                                         \
  (M)[1] = -(M)[1] / (DET);                             \
  (M)[2] = -(M)[2] / (DET);

#define MV_PROD_PLUS_WELL_2x2(M,V,DT,GAS_MULT,OIL_MULT,R)                       \
  (R)[0] += (GAS_MULT) * (DT) * ((M)[0] * (V)[0] + (M)[1] * (V)[1]);     \
  (R)[1] += (OIL_MULT) * (DT) * ((M)[2] * (V)[0] + (M)[3] * (V)[1]);

#define MV_PROD_PLUS_WELL_3x3(M,V,DT,GAS_MULT,OIL_MULT,R)                          \
  (R)[0] +=              (DT) * ((M)[0] * (V)[0] + (M)[1] * (V)[1] + (M)[2] * (V)[2]);     \
  (R)[1] += (GAS_MULT) * (DT) * ((M)[3] * (V)[0] + (M)[4] * (V)[1] + (M)[5] * (V)[2]);     \
  (R)[2] += (OIL_MULT) * (DT) * ((M)[6] * (V)[0] + (M)[7] * (V)[1] + (M)[8] * (V)[2]);

#define MV_PROD_PLUS_WELL(IS_W,IS_G,IS_O,M,V,DT,GAS_MULT,OIL_MULT,R)              \
  if ((IS_W) && (IS_G) && (IS_O))                                        \
    { MV_PROD_PLUS_WELL_3x3(M,V,DT,GAS_MULT,OIL_MULT,R); }                        \
  else if ((IS_W) && (IS_O))                                             \
    { MV_PROD_PLUS_WELL_2x2(M,V,DT,1.0,OIL_MULT,R); }                             \
  else if ((IS_G) && (IS_O))                                             \
    { MV_PROD_PLUS_WELL_2x2(M,V,DT,GAS_MULT,OIL_MULT,R); }                        \
  else if ((IS_W) || (IS_O))                                             \
    { (R)[0] += (OIL_MULT) * (DT) * (M)[0] * (V)[0]; }                              \
  else if ((IS_G))                                                       \
    { (R)[0] += (GAS_MULT) * (DT) * (M)[0] * (V)[0]; }                 \
  else return -1;

#define MM_PROD_PLUS_WELL_2x2(M1,M2,DT,GAS_MULT,OIL_MULT,R)                           \
  (R)[0] += (GAS_MULT) * (DT) * ((M1)[0] * (M2)[0] + (M1)[1] * (M2)[2]);     \
  (R)[1] += (GAS_MULT) * (DT) * ((M1)[0] * (M2)[1] + (M1)[1] * (M2)[3]);     \
  (R)[2] += (OIL_MULT) * (DT) * ((M1)[2] * (M2)[0] + (M1)[3] * (M2)[2]);                  \
  (R)[3] += (OIL_MULT) * (DT) * ((M1)[2] * (M2)[1] + (M1)[3] * (M2)[3]);

#define MM_PROD_PLUS_WELL_3x3(M1,M2,DT,GAS_MULT,OIL_MULT,R)            \
  (R)[0] += (DT) * ((M1)[0] * (M2)[0] + (M1)[1] * (M2)[3] + (M1)[2] * (M2)[6]);     \
  (R)[1] += (DT) * ((M1)[0] * (M2)[1] + (M1)[1] * (M2)[4] + (M1)[2] * (M2)[7]);     \
  (R)[2] += (DT) * ((M1)[0] * (M2)[2] + (M1)[1] * (M2)[5] + (M1)[2] * (M2)[8]);     \
  (R)[3] += (GAS_MULT) * (DT) * ((M1)[3] * (M2)[0] + (M1)[4] * (M2)[3] + (M1)[5] * (M2)[6]);     \
  (R)[4] += (GAS_MULT) * (DT) * ((M1)[3] * (M2)[1] + (M1)[4] * (M2)[4] + (M1)[5] * (M2)[7]);     \
  (R)[5] += (GAS_MULT) * (DT) * ((M1)[3] * (M2)[2] + (M1)[4] * (M2)[5] + (M1)[5] * (M2)[8]);     \
  (R)[6] += (OIL_MULT) * (DT) * ((M1)[6] * (M2)[0] + (M1)[7] * (M2)[3] + (M1)[8] * (M2)[6]);     \
  (R)[7] += (OIL_MULT) * (DT) * ((M1)[6] * (M2)[1] + (M1)[7] * (M2)[4] + (M1)[8] * (M2)[7]);     \
  (R)[8] += (OIL_MULT) * (DT) * ((M1)[6] * (M2)[2] + (M1)[7] * (M2)[5] + (M1)[8] * (M2)[8]);

#define MM_PROD_PLUS_WELL(IS_W,IS_G,IS_O,M1,M2,DT,GAS_MULT,OIL_MULT,R)            \
  if ((IS_W) && (IS_G) && (IS_O))                                        \
    { MM_PROD_PLUS_WELL_3x3(M1,M2,DT,GAS_MULT,OIL_MULT,R); }                    \
  else if ((IS_W) && (IS_O))                                             \
    { MM_PROD_PLUS_WELL_2x2(M1,M2,DT,1.0,OIL_MULT,R); }                         \
  else if ((IS_G) && (IS_O))                                             \
    { MM_PROD_PLUS_WELL_2x2(M1,M2,DT,GAS_MULT,OIL_MULT,R); }                    \
  else if ((IS_W) || (IS_O))                                             \
    { (R)[0] += (DT) * (M1)[0] * (M2)[0]; }                              \
  else if ((IS_G))                                                       \
    { (R)[0] += (GAS_MULT) * (DT) * (M1)[0] * (M2)[0]; }                 \
  else return -1;

#ifndef COMPOSITIONAL
#define MM_SUM_2x2(M,R)           \
  (R)[0] += (M)[0];               \
  (R)[1] += (M)[1];               \
  (R)[2] += (M)[2];               \
  (R)[3] += (M)[3];

#define MM_SUM_3x3(M,R)           \
  (R)[0] += (M)[0];               \
  (R)[1] += (M)[1];               \
  (R)[2] += (M)[2];               \
  (R)[3] += (M)[3];               \
  (R)[4] += (M)[4];               \
  (R)[5] += (M)[5];               \
  (R)[6] += (M)[6];               \
  (R)[7] += (M)[7];               \
  (R)[8] += (M)[8];

#define MM_SUM(F,M,R)          \
  if ((F) == 3)                \
    { MM_SUM_3x3 ((M), (R)) }  \
  else if ((F) == 2)           \
    { MM_SUM_2x2 ((M), (R)) }  \
  else if ((F) == 1)           \
    (R)[0] += (M)[0];
#else
#define MM_SUM(F,M,R)          \
  mm_prod ((F), (M), (R));
#endif


//! A=A-B*C; A(n*1); B(n*1); C(1)
#define V_MINUS_VS_PROD_1x1(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];
//! A=A-B*C; A(n*1); B(n*1); C(1)
#define V_MINUS_VS_PROD_2x2(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];                                    \
  (A)[1] -= (B)[1] * (C)[0];
//! A=A-B*C; A(n*1); B(n*1); C(1)
#define V_MINUS_VS_PROD_3x3(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];                                    \
  (A)[1] -= (B)[1] * (C)[0];                                    \
  (A)[2] -= (B)[2] * (C)[0];

//! A=A-B*C; A(n*1); B(n*1); C(1)
#define V_MINUS_VS_PROD(F,B,C,A)                                \
  if ((F) == 3)                                                 \
    {V_MINUS_VS_PROD_3x3 (B,C,A)}                               \
  else if ((F) == 2)                                            \
    {V_MINUS_VS_PROD_2x2 (B,C,A)}                               \
  else                                                          \
    {V_MINUS_VS_PROD_1x1 (B,C,A)}

//! A=A-B*C A (n*n); B(n*1); C(1*n)
#define M_MINUS_VV_PROD_1x1(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];

//! A=A-B*C A (n*n); B(n*1); C(1*n)
#define M_MINUS_VV_PROD_2x2(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];                                    \
  (A)[1] -= (B)[0] * (C)[1];                                    \
  (A)[2] -= (B)[1] * (C)[0];                                    \
  (A)[3] -= (B)[1] * (C)[1];

//! A=A-B*C A (n*n); B(n*1); C(1*n)
#define M_MINUS_VV_PROD_3x3(B,C,A)                              \
  (A)[0] -= (B)[0] * (C)[0];                                    \
  (A)[1] -= (B)[0] * (C)[1];                                    \
  (A)[2] -= (B)[0] * (C)[2];                                    \
  (A)[3] -= (B)[1] * (C)[0];                                    \
  (A)[4] -= (B)[1] * (C)[1];                                    \
  (A)[5] -= (B)[1] * (C)[2];                                    \
  (A)[6] -= (B)[2] * (C)[0];                                    \
  (A)[7] -= (B)[2] * (C)[1];                                    \
  (A)[8] -= (B)[2] * (C)[2];


//! A=A-B*C A (n*n); B(n*1); C(1*n)
#define M_MINUS_VV_PROD(F,B,C,A)                                \
  if ((F) == 3)                                                 \
    {M_MINUS_VV_PROD_3x3 (B,C,A)}                               \
  else if ((F) == 2)                                            \
    {M_MINUS_VV_PROD_2x2 (B,C,A)}                               \
  else                                                          \
    {M_MINUS_VV_PROD_1x1 (B,C,A)}







#if 0
#define BLOCK_NORM_2x2(M,N)     \
  (N) = fabs((M)[0]);           \
  if ((N) < fabs((M)[1]))       \
  (N) = fabs((M)[1]);         \
  if ((N) < fabs((M)[2]))       \
  (N) = fabs((M)[2]);         \
  if ((N) < fabs((M)[3]))       \
  (N) = fabs((M)[3]);

#define BLOCK_NORM_3x3(M,N)     \
  (N) = fabs((M)[0]);           \
  if ((N) < fabs((M)[1]))       \
  (N) = fabs((M)[1]);         \
  if ((N) < fabs((M)[2]))       \
  (N) = fabs((M)[2]);         \
  if ((N) < fabs((M)[3]))       \
  (N) = fabs((M)[3]);         \
  if ((N) < fabs((M)[4]))       \
  (N) = fabs((M)[4]);         \
  if ((N) < fabs((M)[5]))       \
  (N) = fabs((M)[5]);         \
  if ((N) < fabs((M)[6]))       \
  (N) = fabs((M)[6]);         \
  if ((N) < fabs((M)[7]))       \
  (N) = fabs((M)[7]);         \
  if ((N) < fabs((M)[8]))       \
  (N) = fabs((M)[8]);


#define BLOCK_NORM(F,M,N)       \
  if ((F) == 3)                 \
    {BLOCK_NORM_3x3((M), (N))}  \
  else if ((F) == 2)            \
    {BLOCK_NORM_2x2((M), (N))}  \
  else if ((F) == 1)            \
    {(N) = fabs ((M)[0]);}
#else //0
#define BLOCK_NORM_2x2(M,N)     \
  (N) = (M)[0] * (M)[0];           \
  (N) += ((M)[1]) * ((M)[1]);         \
  (N) += ((M)[2]) * ((M)[2]);         \
  (N) += ((M)[3]) * ((M)[3]);

#define BLOCK_NORM_3x3(M,N)     \
  (N) = ((M)[0]) * ((M)[0]);           \
  (N) += ((M)[1]) * ((M)[1]);         \
  (N) += ((M)[2]) * ((M)[2]);         \
  (N) += ((M)[3]) * ((M)[3]);         \
  (N) += ((M)[4]) * ((M)[4]);         \
  (N) += ((M)[5]) * ((M)[5]);         \
  (N) += ((M)[6]) * ((M)[6]);         \
  (N) += ((M)[7]) * ((M)[7]);         \
  (N) += ((M)[8]) * ((M)[8]);


#define BLOCK_NORM(F,M,N)       \
  if ((F) == 3)                 \
    {BLOCK_NORM_3x3((M), (N))}  \
  else if ((F) == 2)            \
    {BLOCK_NORM_2x2((M), (N))}  \
  else if ((F) == 1)            \
    {(N) = ((M)[0]) * ((M)[0]);}

#endif

#endif //__MATRIX_MACROSES_H_
