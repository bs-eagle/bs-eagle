/**
 *       \file  py_table_2d_debug.h
 *      \brief  Python wrappers for table_2d
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.05.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef PY_BS_TABLE_2D_H_
#define PY_BS_TABLE_2D_H_

#ifdef _DEBUG

#include "table_2d.h"

namespace blue_sky
  {
  namespace python
    {

    class py_data_row
      {
      public:
        typedef table::data_row  data_row_t;

        py_data_row ()
        {

        }

        py_data_row (const data_row_t &r)
            : row (r)
        {

        }

        typename strategy_t::item_t get (int column_index) const
          {
            return row[column_index];
          }

      private:

        data_row_t row;

      };

    class py_data_group
      {
      public:
        typedef table::data_group  data_group_t;

        py_data_group ()
            : group (0, 0, 0)
        {

        }

        py_data_group (const data_group_t &g)
            : group (g)
        {

        }

        int get_rows_count () const
          {
            return group.get_rows_count ();
          }
        int get_columns_count () const
          {
            return group.get_columns_count ();
          }

        py_data_row get_row (int row_index) const
          {
            return py_data_row (group.get_row (row_index));
          }

      private:


        data_group_t group;

      };

    class py_table_2d
      {
      public:
        typedef table::table_2d  table_2d_t;

        py_table_2d ()
            : data (0)
        {

        }

        py_table_2d (const table_2d_t &d)
            : data (d)
        {

        }

        int get_groups_count () const
          {
            return data.get_groups_count ();
          }

        py_data_group get_data_group (int group_index)
        {
          return py_data_group (data.get_data_group (group_index));
        }

      private:

        table_2d_t data;

      };

    void
    py_export_table_2d ();


  } // namespace python
} // namespace blue_sky


#endif  // #ifdef _DEBUG


#endif // #ifdef PY_BS_TABLE_2D_H_
