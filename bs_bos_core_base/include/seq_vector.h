/**
 * \file seq_vector.h
 * \brief vector of sequences data (inhereted from std::vector (or bs_array). introduced to hide difference between std::vector and mpi_vector
 * \author Sergey Miryanov
 * \date 27.05.2008
 * */
#ifndef BS_SEQ_VECTOR_H_
#define BS_SEQ_VECTOR_H_

#include "aligned_allocator.h"
#include "shared_vector.h"

namespace blue_sky
{
#ifdef USE_STD_VECTOR_IN_SEQ_VECTOR_
  template <class T>
  class seq_vector : public std::vector <T, aligned_allocator <T, 16> >
  {
  public:

    typedef typename std::vector<T, aligned_allocator <T, 16> > base_t;

    typedef typename base_t::iterator iterator;
    typedef typename base_t::const_iterator const_iterator;

    template <typename Y>
    struct array
    {
      typedef seq_vector <Y> type;
    };

    using base_t::swap;

    seq_vector ()
    : base_t ()
    {
    }

    seq_vector (size_t ns)
    : base_t (ns)
    {
    }

    seq_vector (iterator first, iterator end)
    : base_t (first, end)
    {
    }

    template <class mx_t>
    int
    init_by_matrix (const mx_t *mx)
    {
      base_t::assign (mx->n_rows * mx->n_block_size, 0);
      return 0;
    }
  };
#else
  template <class T>
  class seq_vector : public shared_vector <T, aligned_allocator <T, 16> >
  {
  public:

    typedef shared_vector <T, aligned_allocator <T, 16> > base_t;
    typedef typename base_t::iterator                     iterator;
    typedef typename base_t::const_iterator               const_iterator;

    template <typename Y>
    struct array
    {
      typedef seq_vector <Y> type;
    };

    using base_t::swap;

    seq_vector ()
    : base_t ()
    {
    }

    seq_vector (size_t ns)
    : base_t (ns)
    {
    }

    seq_vector (iterator first, iterator end)
    : base_t (first, end)
    {
    }

    template <class mx_t>
    int
    init_by_matrix (const mx_t *mx)
    {
      base_t::assign (mx->n_rows * mx->n_block_size, 0);
      return 0;
    }
  };
#endif

} // namespace blue_sky


#endif  // #ifndef BS_SEQ_VECTOR_H_
