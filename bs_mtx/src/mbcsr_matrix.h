/**
 * @file mbcsr_matrix.h
 * @brief Block CSR multi matrix implementation
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-27
 */
#ifndef __MBCSR_MATRIX_H
#define __MBCSR_MATRIX_H

#include "mbcsr_matrix_iface.h"

#include <map>
#include <string>

namespace blue_sky
{
  /**
   * @brief interface class for block CSR matrix storage and manipulation
   */

  class BS_API_PLUGIN mbcsr_matrix: public mbcsr_matrix_iface
    {
    public:

      typedef bcsr_matrix_iface                                 bcsr_iface_t;
      typedef smart_ptr<bcsr_iface_t, true>                     sp_bcsr_iface_t;
      typedef mbcsr_matrix                                      this_t;
      typedef matrix_iface                                      base_t;
      typedef std::map<std::string, sp_bcsr_iface_t >            map_t;

      //! destructor
      virtual ~mbcsr_matrix ()
      {
      }

      //-----------------------------------------
      //  matrix_iface METHODS
      //-----------------------------------------

      //! calculate matrix vector product, v -- input vector, r -- output vector
      //! r += A * v
      //! return 0 if success
      virtual int matrix_vector_product (spv_double v, spv_double r) const
      {
        int f = 0;
        map_t::const_iterator   ib, ie, i;

        ib = mat_map.begin ();
        ie = mat_map.end ();

        for (i = ib; i != ie; ++i)
          {
            f = f || i->second->matrix_vector_product (v, r);
          }
        return f;
      }

      //! calculate matrix^t vector product, v -- input vector, r -- output vector
      //! r += A^T * v
      //! return 0 if success
      virtual int matrix_vector_product_t (spv_double v, spv_double r) const
      {
        int f = 0;
        map_t::const_iterator   ib, ie, i;

        ib = mat_map.begin ();
        ie = mat_map.end ();

        for (i = ib; i != ie; ++i)
          {
            f = f || i->second->matrix_vector_product_t (v, r);
          }
        return f;
      }

      //! calculate linear combination r = alpha * Au + beta * v
      //! alpha, beta -- scalar
      //! v, u -- input vector
      //! r -- output vector
      //! return 0 if success
      virtual int calc_lin_comb (t_double alpha, t_double beta, spv_double u, spv_double v, spv_double r) const;

      //! return total amount of allocated memory in bytes
      virtual t_double get_allocated_memory_in_mbytes () const
      {
        t_double f = 0;
        map_t::const_iterator   ib, ie, i;

        ib = mat_map.begin ();
        ie = mat_map.end ();

        for (i = ib; i != ie; ++i)
          {
            f = f + i->second->get_allocated_memory_in_mbytes ();
          }
        return f;
      }

      //! return block size
      virtual t_long get_n_block_size () const
      {
        if (mat_map.size ())
          return mat_map.begin ()->second->get_n_block_size ();
        else
          return 0;
      }

      //! return number of rows in matrix
      virtual t_long get_n_rows () const
      {
        if (mat_map.size ())
          return mat_map.begin ()->second->get_n_rows ();
        else
          return 0;
      }

      //! return number of cols in matrix (return -1 in number of columns is unknown)
      virtual t_long get_n_cols () const
      {
        if (mat_map.size ())
          return mat_map.begin ()->second->get_n_cols ();
        else
          return 0;
      }

      //! return true if number of rows is equal to the number of columns
      virtual bool is_square () const
        {
        if (mat_map.size ())
          return mat_map.begin ()->second->is_square ();
        else
          return false;
        }

      //! initialize vector
      virtual void init_vector (spv_double v) const
        {
        if (mat_map.size ())
          mat_map.begin ()->second->init_vector (v);
        }

      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      // check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const
        {
          int f = 0;
          map_t::const_iterator   ib, ie, i;

          ib = mat_map.begin ();
          ie = mat_map.end ();

          for (i = ib; i != ie; ++i)
            {
              f = f || i->second->internal_check ();
            }
          return f;
        }

      //! return smart pointer to the matrix with name
      virtual sp_bcsr_iface_t get_matrix (const std::string &name)
        {
          map_t::iterator it = mat_map.find (name);

          if (it == mat_map.end ())
            throw "Error";
          return it->second;
        }

      //! add new BCSR matrix to the list
      virtual void add_matrix (const std::string &name, sp_bcsr_iface_t m)
        {
          map_t::iterator it;

          if (!mat_map.size ())
            {
              mat_map.insert (std::pair <std::string, sp_bcsr_iface_t> (name, m));
            }
          else if (mat_map.begin ()->second->get_n_rows () != m->get_n_rows ()
                   || mat_map.begin ()->second->get_n_cols () != m->get_n_cols ()
                   || mat_map.begin ()->second->get_n_block_size () != m->get_n_block_size ())
            {
              // TODO: fix
              throw "Matrix size Error";
            }
          else
            {
              it = mat_map.find (name);

              if (it == mat_map.end ())
                mat_map.insert (std::pair <std::string, sp_bcsr_iface_t> (name, m));
              else
                mat_map[name] = m;
            }
        }

      //! clear all
      virtual void clear ()
        {
          mat_map.clear ();
        }


      /**
       * @brief merge all matrixies together (A = SUM_i (A_i) )and apply filter
       *        filter is a array length n_rows with elements 0 or 1
       *        if 0 in the i position than merged matrix i-th row should include
       *        only diagonal element, the same for the i-th column
       *
       * @param filter -- <INPUT> array length n_rows, elements should be 0 or 1
       *
       * @return merged BCSR matrix
       */
      virtual sp_bcsr_iface_t merge (spv_int filter);

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;
          map_t::const_iterator   ib, ie, i;
          t_long count;

          ib = mat_map.begin ();
          ie = mat_map.end ();

          for (i = ib, count = 0; i != ie; ++i, ++count)
            {
              s << "Matrix: " << count << " -- " << i->first << std::endl;
              s << i->second->py_str ();
            }

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN
    public:

    protected:

      map_t             mat_map;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (mbcsr_matrix);

    };



}//namespace blue_sky

#endif // __MBCSR_MATRIX_H
