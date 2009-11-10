/** 
 * \file lin_solv_macro.h
 * \brief
 * \author 
 * \date 
 * */
#ifndef BS_LIN_SOLV_MACRO_H_
#define BS_LIN_SOLV_MACRO_H_

namespace blue_sky
  {
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

  } // end namespace

#endif  // #ifndef BS_LIN_SOLV_MACRO_H_

