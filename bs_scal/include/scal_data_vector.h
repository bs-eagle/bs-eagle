/**
 * \file scal_data_vector.h
 * \brief vector that "hold" scal data and iterator that allow iterate over these data
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#ifndef BS_SCAL_DATA_VECTOR_H_
#define BS_SCAL_DATA_VECTOR_H_

namespace blue_sky
  {
  //namespace scal
  //  {

    template <typename strategy_t>
    struct data_vector
      {
public:
        typedef typename strategy_t::item_t item_t;

        struct iterator : public std::iterator <std::random_access_iterator_tag, item_t>
          {
public:
            typedef std::iterator <std::random_access_iterator_tag, item_t> base_t;
            typedef typename base_t::value_type                             value_type;
            typedef typename base_t::difference_type                        difference_type;

            iterator()
                : data (0)
                , step (0)
                , position (0)
            {}

            iterator (const item_t *d, int s, int p = 0)
                : data (d)
                , step (s)
                , position (p)
            {
              if (position < 0)
                position = 0;
            }

            value_type operator * () const
              {
                BS_ASSERT (data);
                return data [position];
              }

            iterator operator ++ (int)
            {
              iterator it (data, step, position);
              position += step;
              return it;
            }
            iterator operator ++ ()
            {
              position += step;
              return iterator (data, step, position);
            }

            void operator += (difference_type i)
            {
              position += int(step * i);
            }

            iterator operator + (size_t i) const
              {
                return iterator (data, step, position + i * step);
              }
            iterator operator - (size_t i) const
              {
                return iterator (data, step, position - i * step);
              }

            difference_type operator - (const iterator &it) const
              {
                return (position - it.position) / step;
              }

            bool operator == (const iterator &it) const
              {
                BS_ASSERT (data == it.data);
                return position == it.position;
              }
            bool operator != (const iterator &it) const
              {
                BS_ASSERT (data == it.data);
                return position != it.position;
              }

private:

            const item_t    *data;
            int             step;
            int             position;
          };


        typedef iterator const_iterator;
        typedef item_t value_type;

        data_vector (const value_type *data, int step, int count)
            : data (data)
            , step (step)
            , count (count)
        {
        }

        value_type &operator [] (size_t index)
        {
          BS_ASSERT (data);
          return const_cast <value_type&> (data[index * step]); // TODO: BUG:
        }

        const value_type &operator [] (size_t index) const
          {
            BS_ASSERT (data);
            return data[index * step];
          }

        const_iterator begin () const
          {
            return const_iterator (data, step);
          }
        const_iterator end () const
          {
            return const_iterator (data, step, step * count);
          }

        size_t size () const
          {
            return count;
          }

        value_type front () const
          {
            return data[0];
          }
        value_type back () const
          {
            return data[(count - 1) * step];
          }

private:

        const value_type		*data;
        size_t              step;
        size_t              count;
      };


//  } // namespace scal
} // namespace blue_sky

#endif  // #ifndef BS_SCAL_DATA_VECTOR_H_
