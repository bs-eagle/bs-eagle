#ifndef __MV_FUN_H
#define __MV_FUN_H
/** 
 * @file mv_functions.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-09-10
 */
namespace blue_sky {
template <class vector_v1_t, class vector_v2_t> inline typename vector_v1_t::value_type
mv_vector_inner_product (const vector_v1_t &v1, const vector_v2_t &v2, int /* obsolete */ = 0)
    {
      typename vector_v1_t::value_type sum = 0;
      size_t i = 0;
      size_t n = v1.size ();

      if (v1.size () != v2.size ())
        return 0;

#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
#pragma omp parallel for reduction (+: sum)
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL
      for (i = 0; i < n; ++i)
        {
          sum += v1[i] * v2[i];
        }

      return sum;
    }
/** 
 * @file mv_functions.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-09-10
 */
template <class T, class TI> inline T
mv_vector_inner_product_n (const T *v1, const T *v2, TI n)
  {
    T sum = 0;
    TI i = 0;

#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
#pragma omp parallel for reduction (+: sum)
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL
    for (i = 0; i < n; ++i)
      {
        sum += v1[i] * v2[i];
      }

    return sum;
  }

    /**
    * @brief sum vector
    *
    * @param A -- first vector
    * @param B -- second vector
    * @param alpha -- first scalar
    * @param beta -- second scalar
    * @param RES -- result vector
    * @param I -- index
    * @param N -- vector length

    *
    * @return nothing
    */
    #define SUM_VECTOR(A,B,alpha,beta,RES,I,N)                           \
    for ((I) = 0; (I) < (N); ++(I))                       \
    (RES)[(I)] = alpha*(A)[(I)] + beta*(B)[(I)];

    template <class vector_t>
    inline void sum_vector (vector_t &x, typename vector_t::value_type alpha, vector_t &y, typename vector_t::value_type beta, vector_t &res)
    {
      //int i = 0, cnt = (int)a.size ();

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      BS_ASSERT (x.size ());
      BS_ASSERT (y.size ());
      BS_ASSERT (res.size ());
      BS_ASSERT (x.size () == y.size ());
      BS_ASSERT (x.size () == res.size ());

      size_t i = 0, cnt = x.size ();

      for (i = 0; i < cnt; ++i)
        {
          res[i] = alpha * x[i] + beta * y[i];
        }
    }
    template <class T, class TI>
    inline void sum_vector_n (T *x, T alpha, T *y, T beta, T *res, TI n)
    {
      //int i = 0, cnt = (int)a.size ();

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      TI i = 0;

      for (i = 0; i < n; ++i)
        {
          res[i] = alpha * x[i] + beta * y[i];
        }
    }




    /**
    * @brief scale vector
    *
    * @param A -- vector
    * @param T -- scale factor
    * @param I -- index
    * @param N -- vector length
    *
    * @return nothing
    */
    #define SCALE_VECTOR(A,T,I,N)                           \
    for ((I) = 0; (I) < (N); ++(I))                       \
    (A)[(I)] *= (T);

    // if we got performance decrease we have to specify macroses for seq and mpi separately
    // used only in gmres_solver2
    template <class vector_t>
    inline void scale_vector (vector_t &a, typename vector_t::value_type t)
    {
      int i = 0, cnt = (int)a.size ();

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      for (i = 0; i < cnt; ++i)
        {
          a[i] *= t;
        }
    }

    template <class T, class TI>
    inline void scale_vector_n (T *a, T t, TI n)
    {
      TI i = 0;

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      for (i = 0; i < n; ++i)
        {
          a[i] *= t;
        }
    }

    /**
    * @brief AXPY (X = X + T * Y)
    *
    * @param X -- vector
    * @param T -- scale factor
    * @param Y -- vector
    * @param I -- index
    * @param N -- vectors length
    *
    * @return nothing
    */
    #define AXPY(X,T,Y,I,N)                                 \
    for ((I) = 0; (I) < (N); ++(I))                        \
    (X)[(I)] += (T) * (Y)[(I)];

    // if we get performance decrease we have to specify macroses for seq and mpi separately
    // used only in gmres_solver2
    template <class vector_t>
    inline void axpy (vector_t &x, const vector_t &y, typename vector_t::value_type t)
    {
      BS_ASSERT (x.size () == y.size ());
      BS_ASSERT (x.size ());

      size_t i = 0, cnt = x.size ();
      const size_t unroll_factor = 8;
      size_t cnt2 = cnt - (cnt % unroll_factor);

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      for (i = 0; i < cnt2; i += unroll_factor)
        {
          x[i] += t * y[i];
          x[i + 1] += t * y[i + 1];
          x[i + 2] += t * y[i + 2];
          x[i + 3] += t * y[i + 3];
          x[i + 4] += t * y[i + 4];
          x[i + 5] += t * y[i + 5];
          x[i + 6] += t * y[i + 6];
          x[i + 7] += t * y[i + 7];
        }
      for (; i < cnt; ++i)
        {
          x[i] += t * y[i];
        }
    }
    template <class T, class TI>
    inline void axpy_n (T *x, const T *y, T t, TI n)
    {
      TI i = 0;
      const TI unroll_factor = 8;
      TI cnt2 = n - (n % unroll_factor);

      #ifdef GMRES_SOLVER2_SOLVE_PARALLEL
      #pragma omp parallel for
      #endif
      for (i = 0; i < cnt2; i += unroll_factor)
        {
          x[i] += t * y[i];
          x[i + 1] += t * y[i + 1];
          x[i + 2] += t * y[i + 2];
          x[i + 3] += t * y[i + 3];
          x[i + 4] += t * y[i + 4];
          x[i + 5] += t * y[i + 5];
          x[i + 6] += t * y[i + 6];
          x[i + 7] += t * y[i + 7];
        }
      for (; i < n; ++i)
        {
          x[i] += t * y[i];
        }
    }

    /**
    * @brief AXPY_AYPX // X = Z + P * (X + T * Y)
    *
    * @param X -- vector
    * @param T -- scale factor
    * @param Y -- vector
    * @param P -- scale factor
    * @param Z -- vector
    * @param I -- index
    * @param N -- vectors length
    *
    * @return nothing
    */
    #define AXPY_AYPX(X,T,Y,P,Z,I,N)                             \
    for ((I) = 0; (I) < (N); ++(I))                             \
    (X)[(I)] = (Z)[(I)] + (P) * ((X)[(I)] + (T) * (Y)[(I)]);

    template <class vector_t> inline void
    axpy_aypx (vector_t &x, typename vector_t::value_type t, const vector_t &y, typename vector_t::value_type p, const vector_t &z)
    {
      BS_ASSERT (x.size ());
      BS_ASSERT (y.size ());
      BS_ASSERT (z.size ());
      BS_ASSERT (x.size () == y.size ());
      BS_ASSERT (x.size () == z.size ());

      size_t i = 0, cnt = x.size ();

      for (i = 0; i < cnt; ++i)
        {
          x[i] = z[i] + p * (x[i] + t * y[i]);
        }
    }

    template <class T, class TI> inline void
    axpy_aypx_n (T *x, T t, const T *y, T p, const T *z, TI n)
    {
      TI i = 0;

      for (i = 0; i < n; ++i)
        {
          x[i] = z[i] + p * (x[i] + t * y[i]);
        }
    }
} // end of namespace
#endif //__MV_FUN_H

