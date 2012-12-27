/** 
 * @file vartype_table.h
 * @brief implementation of the table interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifndef VARTYPE_TABLE_X6PKPWC9

#define VARTYPE_TABLE_X6PKPWC9

#include <string>
#include <sstream>
#include <vector>

#include "vartype_table_iface.h"

namespace blue_sky
{
  template <class var_type_t>
  class BS_API_PLUGIN vartype_table : public vartype_table_iface <var_type_t>
    {
    
    public: 

      typedef typename vartype_table_iface<var_type_t>::vector_t      vector_t;
      typedef bs_array <var_type_t, BS_EAGLE_ARRAY_TRAITS>            var_type_array_t;
      typedef smart_ptr <var_type_array_t, true >            sp_var_type_array_t;
      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      ~vartype_table ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:

      /** 
       * @brief Initialize or reinitialize table by <n_rows> <n_cols> 
       * 
       * @param n_cols -- <INPUT> number of cols in table
       * 
       * @return 0 if success
       */
      virtual int init (const t_long n_cols);

      /** 
       * @brief set <name> of the column <col>
       * 
       * @param col  -- <INPUT> column index
       * @param name -- <INPUT> new name
       * 
       * @return 0 if success
       */
      virtual int set_col_name (const t_long col, const std::wstring &name)
        {
          col_names[col] = name;
          return 0;
        }

      /** 
       * @brief return name of the <col> column
       * 
       * @param col -- <INPUT> column index
       * 
       * @return name
       */
      virtual std::wstring get_col_name (const t_long col) const
        {
          return col_names[col];
        }

      /** 
       * @brief return pointer to the column <col> data 
       * 
       * @param col -- <INPUT> column index
       * 
       * @return pointer
       */
      virtual var_type_t *get_col_ptr (const t_long col)
        {
          return values[col].size () ? &(values[col])[0] : 0;
        }

      /** 
       * @brief return pointer to the column <col> data 
       * 
       * @param col -- <INPUT> column index
       * 
       * @return 
       */
      virtual vector_t & get_col_vector (const t_long col)
        {
          return values[col];
        }  
      /** 
       * @brief add new column and return reference to the column <col> data 
       * 
       * @param name -- <INPUT> new name       
       * 
       * @return last added column index
       */
      virtual t_long add_col_vector (const t_long col, const std::wstring &name, sp_var_type_array_t new_vector)
        {
          var_type_array_t &new_vector_ = *new_vector;
          t_long size = values.size ();
          BS_ASSERT (col >= 0 && col < size);
          
          t_long new_vec_size = new_vector_.size ();
          values[col].resize (new_vec_size);
          for (t_long i = 0; i < new_vec_size; ++i)
            {
              values[col][i] = new_vector_[i];
            } 
          set_col_name (col, name);
          
          return 0;
        }
        
      /** 
       * @brief get numver of rows in table
       */
      virtual t_long get_n_rows (const t_long col) const
        {
          return values[col].size ();
        }
      
      /** 
       * @brief get number of columns in table
       */
      virtual t_long get_n_cols () const
        {
          return values.size ();
        }

      /** 
       * @brief clear table data 
       */
      virtual void clear ()
        {
          values.clear ();
          col_names.clear ();
        }

      virtual void remove_col_vector (t_long col)
      {
        BS_ASSERT (col >= 0 && col < static_cast <t_long> (values.size ()));
        values[col].clear ();
      }
#if 0
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

      /** 
       * @brief set column <col> values by numpy array
       * 
       * @param col     -- <INPUT> column index
       * @param values  -- <INPUT> numpy array
       */
      virtual void set_col_values (const t_long col, spv_double values);

      /** 
       * @brief return numpy array for given column <col>
       * 
       * @param col -- <INPUT> column index
       * 
       * @return numpy array
       */
      virtual spv_double get_col_values (const t_long col) const;
#endif //BSPY_EXPORTING_PLUGIN
#endif // 0      
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      std::vector<std::wstring> col_names;
      std::vector<std::vector<var_type_t> > values;

      BLUE_SKY_TYPE_DECL (vartype_table);
    };

}; //end of blue_sky namespace


#endif /* end of include guard: TABLE_X6PKPWC9 */
