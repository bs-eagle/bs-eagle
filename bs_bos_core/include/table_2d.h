/**
 * \file table_2d.h
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#ifndef BS_TABLE_2D_H_
#define BS_TABLE_2D_H_

#include "seq_vector.h"

namespace blue_sky
  {
  namespace table
    {

    template <typename strategy_t> class data_group;
    template <typename strategy_t> class data_row;
    template <typename strategy_t> class table_2d;
    template <typename strategy_t> class data_row_push_back;

    template <typename strategy_t>
    class data_row
      {
      public:
        typedef typename strategy_t::item_t        item_t;
        typedef typename strategy_t::item_array_t  item_array_t;

        typedef data_row <strategy_t>              this_t;

      private:
        friend class data_group <strategy_t>;

        data_row (item_t *data_ptr, int columns_count)
            : data_ptr (data_ptr)
            , columns_count (columns_count)
        {
        }

      public:

        data_row ()
            : data_ptr (0)
            , columns_count (0)
        {
        }

        data_row (const this_t &row)
        {
          *this = row;
        }
        data_row &operator= (const this_t &row)
        {
          if (this != &row)
            {
              data_ptr = row.data_ptr;
              columns_count = row.columns_count;
            }
          return *this;
        }

        inline item_t &operator[] (int column_index)
        {
          BS_ASSERT (column_index < columns_count);
          BS_ASSERT (column_index >= 0);
          BS_ASSERT (data_ptr);
#ifdef _DEBUG
          static item_t dummy = 0.0;

          if (column_index >= 0 && column_index < columns_count && data_ptr)
            return data_ptr[column_index];

          throw bs_exception ("data_row", "out of range");
          return dummy;
#else
          return data_ptr[column_index];
#endif
        }
        inline const item_t &operator[] (int column_index) const
          {
            BS_ASSERT (column_index < columns_count);
            BS_ASSERT (column_index >= 0);
            BS_ASSERT (data_ptr);
#ifdef _DEBUG
            static item_t dummy = 0.0;

            if (column_index >= 0 && column_index < columns_count && data_ptr)
              return data_ptr[column_index];

            throw bs_exception ("data_row", "out of range");
            return dummy;
#else
            return data_ptr[column_index];
#endif
          }

      private:

        item_t    *data_ptr;
        int       columns_count;

      };

    template <typename strategy_t>
    class data_group
      {
      public:
        typedef typename strategy_t::item_t        item_t;
        typedef typename strategy_t::item_array_t  item_array_t;

        typedef data_group<strategy_t>             this_t;
        typedef data_row<strategy_t>               data_row_t;

      private:

        friend class table_2d <strategy_t>;

      public:
        data_group (item_t *data_ptr = 0, int rows_count = 0, int columns_count = 0)
            : data_ptr (data_ptr)
            , rows_count (rows_count)
            , columns_count (columns_count)
        {
          row_size = columns_count;
        }

        inline void      set_rows_count (int count)
        {
          rows_count = count;
        }
        inline void      set_data_ptr (item_t *ptr)
        {
          data_ptr = ptr;
        }

      public:

        data_group (const this_t & group)
        {
          *this = group;
        }
        this_t &operator= (const this_t &group)
        {
          if (this != &group)
            {
              data_ptr      = group.data_ptr;
              rows_count    = group.rows_count;
              columns_count = group.columns_count;
              row_size      = group.row_size;
            }

          return *this;
        }

        inline data_row_t  get_row (int row_index) const
          {
#ifdef _DEBUG
            if (row_index >= 0 && row_index < rows_count)
              return data_row_t (data_ptr + row_index * row_size, columns_count);

            throw bs_exception ("data_row", "out of range");
#else
            return data_row_t (data_ptr + row_index * row_size, columns_count);
#endif
          }

        inline int       get_rows_count () const
          {
            return rows_count;
          }
        inline int       get_columns_count () const
          {
            return columns_count;
          }

        inline void      set_columns_count (int count)
        {
          columns_count = count;
          row_size = columns_count;
        }

      private:

        item_t    *data_ptr;
        int       rows_count;
        int       columns_count;
        int       row_size;

      };

    template <typename strategy_t>
    class data_row_push_back
      {
      public:
        typedef typename strategy_t::item_t        item_t;
        typedef typename strategy_t::item_array_t  item_array_t;

        typedef data_row_push_back <strategy_t>    this_t;

      private:
        friend class table_2d <strategy_t>;

        data_row_push_back (item_t *data_ptr, int columns_count)
            : data_ptr (data_ptr)
            , columns_count (columns_count)
        {

        }

      public:

        this_t &push_back (int column_index, item_t value)
        {
          BS_ASSERT (data_ptr);
          BS_ASSERT (column_index < columns_count);
          BS_ASSERT (column_index >= 0);

          if (column_index >= 0 && column_index < columns_count)
            data_ptr[column_index] = value;

          return *this;
        }

      private:


        item_t  *data_ptr;
        int     columns_count;

      };

    template <typename strategy_t>
    class table_2d
      {
      public:

        typedef typename strategy_t::item_t        item_t;
        typedef typename strategy_t::item_array_t  item_array_t;

        typedef typename strategy_t::item_t        index_t;

        typedef data_group <strategy_t>            data_group_t;
        typedef data_row_push_back <strategy_t>    data_row_push_back_t;

        table_2d (int group_count)
            : rows_count (0)
        {
          for (int i = 0; i < group_count; ++i)
            {
              group_list.push_back (data_group_t (0, 0, 0));
            }
        }

        inline void init_dependent (int dependent_rows_count)
        {
          BS_ASSERT (dependent_rows_count);

          int column_count = 0;
          for (int i = 1, cnt = (int)group_list.size (); i < cnt; i++)
            {
              column_count += group_list[i].get_columns_count ();
            }

          BS_ASSERT (column_count);
          BS_ASSERT (data.size ());

          if (data.empty ())
            {
              // TODO: LOG
              BS_ASSERT (false);
              throw bs_exception ("", "data is empty");
            }

          data_t new_data;
          new_data.resize (data.size () + column_count * dependent_rows_count);

          memset (&new_data[0], 0, sizeof (typename data_t::value_type) * new_data.size ());
          memcpy (&new_data[0], &data[0], sizeof (typename data_t::value_type) * data.size ());

          data.swap (new_data);
          item_t *data_ptr = &data[0];

          int rows_count_ = rows_count;
          for (int i = 0, cnt = (int)group_list.size (); i < cnt; i++)
            {
              data_group_t &group = group_list[i];

              group.set_rows_count (rows_count_);
              group.set_data_ptr (data_ptr);

              data_ptr += (rows_count_ * group.get_columns_count ());
              rows_count_ = dependent_rows_count;
            }

#ifdef _DEBUG
          item_t *last_data = &data[0] + data.size ();
          BS_ASSERT (data_ptr == last_data);
#endif
        }
        inline void clear ()
        {
          data.clear ();

          for (int i = 0, cnt = (int)group_list.size (); i < cnt; i++)
            {
              data_group_t &group = group_list[i];

              group.set_rows_count (0);
              group.set_data_ptr (0);
            }
        }

        inline data_group_t &get_data_group (int group_index)
        {
          BS_ASSERT (group_index < (int)group_list.size ());
          BS_ASSERT (group_index >= 0);
#ifdef _DEBUG
          static data_group_t dummy(0, 0, 0);

          if (group_index >= 0 && group_index < (int)group_list.size ()) 
            return group_list[group_index];

          throw bs_exception ("data_row", "out of range");
#else
          return group_list[group_index];
#endif
        }
        inline const data_group_t &get_data_group (int group_index) const
          {
            BS_ASSERT (group_index < (int)group_list.size ());
            BS_ASSERT (group_index >= 0);
#ifdef _DEBUG
            static data_group_t dummy(0, 0, 0);

            if (group_index >= 0 && group_index < (int)group_list.size ()) 
              return group_list[group_index];
            
            throw bs_exception ("data_row", "out of range");
#else
            return group_list[group_index];
#endif
          }
        inline data_group_t get_initial_group ()
        {
          data_group_t group (group_list[0]);
          item_t *data_ptr = &data[0];

          group.set_rows_count (rows_count);
          group.set_data_ptr (data_ptr);

          return group;
        }

        inline int get_groups_count () const
          {
            return (int)group_list.size ();
          }

        inline int get_data_rows_count () const
          {
            return (int)rows_count;
          }

        inline data_row_push_back_t add_row ()
        {
          BS_ASSERT (group_list.size ());

          rows_count += 1;
          data.resize (rows_count * group_list[0].get_columns_count ());

          return data_row_push_back_t (&data[0] + data.size () - group_list[0].get_columns_count (), group_list[0].get_columns_count ());
        }

      private:

        typedef item_array_t              data_t;
        typedef seq_vector <data_group_t> group_list_t;

        data_t        data;
        group_list_t  group_list;
        index_t           rows_count;
      };


  } // namespace table
} // namespace blue_sky


#endif // #ifndef BS_TABLE_2D_H_
