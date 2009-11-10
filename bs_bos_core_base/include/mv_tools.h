#ifndef __MV_TOOLS_H
#define __MV_TOOLS_H

namespace blue_sky
  {
  namespace bos_helper
    {

    namespace helper
      {
      template <typename T> struct is_int
        {
          enum { value = 0};
        };
      template <>						struct is_int <int>
        {
          enum { value = 1};
        };
      template <>						struct is_int <unsigned int>
        {
          enum { value = 1};
        };
      template <>						struct is_int <long>
        {
          enum { value = 1};
        };
      template <>						struct is_int <unsigned long>
        {
          enum { value = 1};
        };
    }

    template <class vector_v1_t, class vector_v2_t> inline typename vector_v1_t::value_type
    mv_vector_inner_product (const vector_v1_t &v1, const vector_v2_t &v2, int /* obsolete */ = 0)
    {
      BOOST_STATIC_ASSERT (helper::is_int <typename vector_v1_t::value_type>::value == 0);
      BOOST_STATIC_ASSERT (helper::is_int <typename vector_v2_t::value_type>::value == 0);

      typename vector_v1_t::value_type sum = 0;
      size_t i = 0;
      size_t n = v1.size ();
      size_t n2 = n - (n % 4);

      BS_ASSERT (v1.size () == v2.size ());

#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
#pragma omp parallel for reduction (+: sum)
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL
      for (i = 0; i < n2; i+=4)
        {
          sum += v1[i + 0] * v2[i + 0];
          sum += v1[i + 1] * v2[i + 1];
          sum += v1[i + 2] * v2[i + 2];
          sum += v1[i + 3] * v2[i + 3];
        }
      for (; i < n; ++i)
        {
          sum += v1[i] * v2[i];
        }

      return sum;
    }

#ifdef _MPI
    template <class T> inline T
    mv_vector_inner_product (const mpi_vector <T> &v1, const mpi_vector <T> &v2, int /* obsolete */ = 0)
    {
      const typename mpi_vector <T>::vector_t &v1_val = v1.get_local_part ();
      const typename mpi_vector <T>::vector_t &v2_val = v2.get_local_part ();

      BS_ASSERT (v1.size () == v2.size ());
      BS_ASSERT (v1.size ());
      BS_ASSERT (v1_val.size () == v2_val.size ());
      BS_ASSERT (v1_val.size ());

      double local_res, res;

      local_res = res = 0.0;
      for (int i = 0, n_local = (int)v1_val.size (); i < n_local; i++)
        local_res += v1_val[i] * v2_val[i];

      MPI_Allreduce (&local_res, &res, 1, mpi_type_t<T>::value, MPI_SUM, MPI_COMM_WORLD);

      return res;
    }
#endif

  } // namespace bos_helper

  template <class strategy_t>
  struct mv_tools
    {
      typedef typename strategy_t::item_t fp_type;
      typedef typename strategy_t::item_array_t item_array_t;

      static inline fp_type
      mv_vector_inner_product2 (const fp_type *v1, const fp_type *v2, const int n)
      {
        int i, istart = 0, iend = n;
        fp_type sum = 0; // double

#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
        int thread_num, n_threads;
        fp_type total_sum = 0; // double
#pragma omp parallel private (i, sum, thread_num, n_threads, istart, iend)
        {
          sum = 0;
          thread_num = omp_get_thread_num ();
          n_threads = omp_get_max_threads ();
          istart =  thread_num * n / n_threads;
          iend = (thread_num + 1) * n / n_threads;
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL

          for (i = istart; i < iend; ++i)
            sum += v1[i] * v2[i];
#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
#pragma omp atomic
          total_sum += sum;
        } //end parallel
        sum = total_sum;
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL
        return sum;
      }

      // r = a + cf1 * b + cf2 * c
      static inline void
      mv_lin_comb_1 (const int n, const fp_type cf1, const fp_type cf2, const fp_type *a,
                     const fp_type *b, const fp_type *c, fp_type *r)
      {
        int i;
#ifdef OTHER_NON_IMPORTANT_PARALLEL
#pragma omp parallel for
#endif //OTHER_NON_IMPORTANT_PARALLEL
        for (i = 0; i < n; ++i)
          r[i] = a[i] + cf1 * b[i] + cf2 * c[i];
      }

      // r = a + cf * b
      static inline void
      mv_lin_comb_2 (const int n, const fp_type cf, const fp_type *a, const fp_type *b, fp_type *r)
      {
        int i;

#ifdef OTHER_NON_IMPORTANT_PARALLEL
#pragma omp parallel for
#endif //OTHER_NON_IMPORTANT_PARALLEL
        for (i = 0; i < n; ++i)
          r[i] = a[i] + cf * b[i];
      }

      static inline void
      mv_set (const int n, const fp_type *a, fp_type *r)
      {
        int i;

        for (i = 0; i < n; ++i)
          r[i] = a[i];
      }

      static inline void
      mv_update_solution (const int n, const int k, const int m, fp_type *h, fp_type *x, fp_type *s, fp_type *v)
      {
        int i, j;
        fp_type *cur_h; //double
        fp_type d; // double
        // Backsolve:
        for (i = k; i >= 0; --i)
          {
            cur_h = h + (m + 1) * i;
            d = s[i] / cur_h[i];
            for (j = i - 1; j >= 0; --j)
              s[j] -= cur_h[j] * d;
            s[i] = d;
          }

        for (j = 0; j <= k; ++j)
          {
            // x += alpha * phat
            mv_lin_comb_2 (n, s[j], x, v + n * j, x);
          }
      }


      static inline void
      mv_vector_print (const fp_type *v, const int n)
      {
        int i;
        for (i = 0; i < n; ++i)
          printf ("%d \t-- %lf\n", i, v[i]);
      }


      static inline void
      mv_vector_print_file (const fp_type *v, const int n, const char *name)
      {
#ifndef UFA_SOLVER
        // TODO: IMPL
#else
        FILE *f;
        int i;
        f = fopen (name, "w");
        for (i = 0; i < n; ++i)
          fprintf (f, "%d \t-- %30.20lf\n", i, v[i]);
        fclose (f);
#endif
      }
    };  // mv_tools

} // namespace blue_sky
#endif //__MV_TOOLS_H

