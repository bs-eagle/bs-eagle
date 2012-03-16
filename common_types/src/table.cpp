/** 
 * @file table.cpp
 * @brief implementation of the table interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif
#include <memory.h>
#include <iomanip>

#include "bs_kernel.h"
#include "table.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>



using namespace std;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{
  table::table (bs_type_ctor_param) 
        : table_iface ()
    {
    }

  table::table (const table& rhs) 
        : bs_refcounter (), table_iface ()
    {
      *this = rhs;
    }

  //! copy 
  int 
  table::copy (const sp_table_t a)
    {
      clear ();
      t_long n_cols = a->get_n_cols ();
      init (a->get_n_rows (), n_cols);
      for (t_long i = 0; i < n_cols; i++)
        {
          set_col_values (i, a->get_col_values (i));
          set_col_name (i, a->get_col_name (i));
        }
      return 0;  
    } 
    
  int 
  table::init (const t_long n_rows, const t_long n_cols)
    {
      values.clear ();
      values.resize (n_cols);
      for (t_long i = 0; i < n_cols; ++i)
        {
          if (n_rows < 1)
            values[i].reserve (20);
          else
            values[i].reserve (2 * n_rows);

          values[i].resize (n_rows);
        }
      col_names.clear ();
      col_names.resize (n_cols);
      return 0;
    }

  void 
  table::add_row (const t_long row_index)
    {
      table_t::iterator i, e;
      t_double v;
      for (i = values.begin (), e = values.end (); i != e; ++i)
        {
          if (row_index == 0)
            {
              v = (*i)[0];
              i->insert (i->begin (), v);
            }
          else if (row_index >= (t_long)i->size ())
            {
              if (i->size () == 0)
                v = 0;
              else
                v = i->back ();
              i->resize (i->size () + 1, v);
            }
          else
            {
              v = 0.5 * ((*i)[row_index - 1] + (*i)[row_index]);
              i->insert (i->begin () + row_index, v);
            }
        }
    }
  
  void 
  table::remove_row (const t_long row_index)
    {
      table_t::iterator i, e;

      for (i = values.begin (), e = values.end (); i != e; ++i)
        {
            if (row_index >= (t_long)i->size ())
                continue;
            i->erase(i->begin() + row_index);
        }
    }
  
  spv_double 
  table::convert_to_array (const t_long n_rows, const t_long n_cols) const 
    {
      //BS_ASSERT (n_rows * n_cols <= values.size ());
      spv_double data = BS_KERNEL.create_object (v_double::bs_type ());
      data->resize (n_rows * n_cols);
      t_double *data_array = &(*data)[0];
      
      for (t_long i = 0; i < n_rows; ++i)
        for (t_long j = 0; j < n_cols; ++j)
          data_array[i * n_cols + j] = get_value (i, j);
          
      return data;    
    } 
     
  void 
  table::convert_from_array (const t_long n_rows, const t_long n_cols, spv_double data) 
    {
      BS_ASSERT (n_rows * n_cols == (int)(data->size ()));
      t_double *data_array = &(*data)[0];
      
      init(0, n_cols);
      
      for (t_long j = 0; j < n_cols; ++j)
        values[j].reserve(n_rows);
        
      for (t_long i = 0; i < n_rows; ++i)
        for (t_long j = 0; j < n_cols; ++j)
          values[j].push_back (data_array[i * n_cols + j]);
    }    

  void 
  table::convert_from_tf_array (const t_long n_rows, const t_long n_cols, t_float *data)
    {
      //BS_ASSERT (n_rows * n_cols == data->size ());
      
      init(0, n_cols);
      
      for (t_long j = 0; j < n_cols; ++j)
        values[j].reserve(n_rows);
        
      for (t_long i = 0; i < n_rows; ++i)
        for (t_long j = 0; j < n_cols; ++j)
          values[j].push_back (data[i * n_cols + j]);
    }

  void 
  table::save (toa_t &ar) const
    {
      //using namespace boost::serialization;
      ar & col_names & values;
    }

  void 
  table::load (tia_t &ar)
    {
      ar & col_names & values;
    }

  table::sp_table_t 
  table::check_serial () const
    {
      std::ostringstream oss;
      std::istringstream iss;

      boost::archive::text_oarchive oar(oss);

      sp_table_t sp_table = BS_KERNEL.create_object ("table");

      save (oar);
      iss.str (oss.str ());
      boost::archive::text_iarchive iar(iss);
      sp_table->load (iar);
      return sp_table;

    }
    
#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  table::to_str () const
    {
      std::ostringstream oss;

      boost::archive::text_oarchive oar(oss);

      save (oar);
      return oss.str ();
    }
  void 
  table::from_str (const std::string &s)
    {
      std::istringstream iss;

      iss.str (s);
      boost::archive::text_iarchive iar(iss);
      load (iar);
    }

  void 
  table::set_col_values (const t_long col, spv_double val)
    {
      if (col < 0 || col >= get_n_cols ()
          || get_n_rows () != (t_long)val->size ())
        {
          // TODO: print error message
          throw "Invalid parameters in set_col_values"; 
        }
      memcpy (&(values[col])[0], &(*val)[0], sizeof (t_double) * get_n_rows ()); 
    }
  
  spv_double 
  table::get_col_values (const t_long col) const
    {
      if (col < 0 || col >= get_n_cols ())
        {
          // TODO: print error message
          throw "Invalid parameters in get_col_values"; 
        }
      spv_double a = BS_KERNEL.create_object (v_double::bs_type ());

      a->resize (get_n_rows ());
      memcpy (&(*a)[0], &(values[col])[0], sizeof (t_double) * get_n_rows ()); 
      return a;
    }

  std::string 
  table::py_str () const
    {
      std::stringstream s;
      std::stringstream s_line;
      t_long i, j;
      const t_long fw = 10; 
      t_long nc = get_n_cols ();
      t_long nr = get_n_rows ();
      
      s << "Table: [" << nr << ", " << nc << "]\n";
      s << "========================\n";
      for (i = 0; i < nc; ++i)
        {
          s << "|" << std::setw (fw) << col_names[i];
          s_line << "|" ;
          for (j = 0; j < fw; ++j)
            {
              s_line << "-";
            }
        }
      s_line << "|" << std::endl;
      s << "|" << std::endl;
      s <<  s_line.str ();
      
      for (j = 0; j < nr; ++j)
        {
          for (i = 0; i < nc; ++i)
            {
              s << "|" << std::setw (10) << std::setprecision (4) << values[i][j];
            }
          s << "|" << std::endl;
        }
      s <<  s_line.str ();
      
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (table);
  BLUE_SKY_TYPE_STD_COPY (table);

  BLUE_SKY_TYPE_IMPL(table,  table_iface, "table", "Table storage", "realization of table storage");

}  // blue_sky namespace
