/** 
 * @file vartype_table.cpp
 * @brief implementation of the vartype_table interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#include <memory.h>
#include <iomanip>

#include "bs_kernel.h"
#include "vartype_table.h"


using namespace std;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{
  template <class var_type_t> 
  vartype_table <var_type_t>::vartype_table (bs_type_ctor_param) 
        : vartype_table_iface<var_type_t> ()
    {
    }

  template <class var_type_t>
  vartype_table <var_type_t>::vartype_table (const vartype_table& rhs) 
        : bs_refcounter (), vartype_table_iface<var_type_t> ()
    {
      *this = rhs;
    }
  
  template <class var_type_t>
  int 
  vartype_table<var_type_t>::init (const t_long n_cols)
    {
      values.clear ();
      values.resize (n_cols);
      col_names.clear ();
      col_names.resize (n_cols);
      return 0;
    }
#if 0
#ifdef BSPY_EXPORTING_PLUGIN
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
#endif 
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(vartype_table, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(vartype_table, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (vartype_table<t_float >), 1,  (vartype_table_iface <t_float>), "float_type_table", "Property storage", "realization of property storage", false);

}  // blue_sky namespace
