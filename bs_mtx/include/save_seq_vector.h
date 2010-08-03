/**
 * \file save_seq_vector.h
 * \brief save seq_vector in file
 * \author Sergey Miryanov
 * \date 18.06.2008
 * */
#ifndef BS_SAVE_SEQ_VECTOR_H_
#define BS_SAVE_SEQ_VECTOR_H_

#include <cstdio>
#include "locale_keeper.h"
#include "bs_assert.h"

namespace blue_sky
  {
  namespace tools
    {

    struct save_shared_vector
      {
        save_shared_vector (const char *filename)
            : lkeeper ("C", LC_ALL)
        {
          file = fopen (filename, "w");
          BS_ASSERT (file) (filename);
        }

        ~save_shared_vector ()
        {
          if (file)
            {
              fflush (file);
              fclose (file);
            }
        }

        template <class array_t>
        save_shared_vector &
        save (const array_t &array)
        {
          for (size_t i = 0, cnt = array.size (); i < cnt; ++i)
            {
              const typename array_t::value_type &v = array[i];
              save_item (file, v);
            }

          return *this;
        }

        template <class functor_t>
        save_shared_vector &
        save_via_fn (const functor_t &fn)
        {
          for (size_t i = 0, cnt = fn.size_i (); i < cnt; ++i)
            {
              for (size_t j = 0, jcnt = fn.size_j (); j < jcnt; ++j)
                {
                  const typename functor_t::value_type &v = fn.get (i, j);
                  save_item (file, v);
                }
            }

          return *this;
        }

        static void
        save_item (FILE *file, int item)
        {
          fprintf (file, "%d\n", item);
        }
        static void
        save_item (FILE *file, float item)
        {
          fprintf (file, "%f\n", item);
        }
        static void
        save_item (FILE *file, double item)
        {
          fprintf (file, "%10.20e\n", item);
        }

private:

        locale_keeper lkeeper;
        FILE *file;
      };


    template <typename item_array_t>
    struct float_saver_fn
      {
        typedef float value_type;

        float_saver_fn (const item_array_t &array_)
            : array_ (array_)
        {
        }

        size_t
        size () const
          {
            return array_.size ();
          }

        value_type
        operator [] (size_t i) const
          {
            return array_[i];
          }

        const item_array_t &array_;
      };

    template <typename item_array_t>
    float_saver_fn <item_array_t>
    float_saver (const item_array_t &array_)
    {
      return float_saver_fn <item_array_t> (array_);
    }


  } // namespace tools
} // namespace blue_sky


#endif  // #ifndef BS_SAVE_SEQ_VECTOR_H_
