/** 
 * @file table.h
 * @brief implementation of the table interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifndef TABLE_X6PKPWC9

#define TABLE_X6PKPWC9

#include <string>
#include <sstream>
#include <vector>

#include "table_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN table : public table_iface
    {
    
    public: 

      typedef std::vector <t_double>                  vector_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      ~table ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------
    public:

      /** 
       * @brief Initialize or reinitialize table by <n_rows> <n_cols> 
       * 
       * @param n_rows -- <INPUT> number of rows in table
       * @param n_cols -- <INPUT> number of cols in table
       * 
       * @return 0 if success
       */
      virtual int init (const t_long n_rows, const t_long n_cols);

      /** 
       * @brief set <name> of the column <col>
       * 
       * @param col  -- <INPUT> column index
       * @param name -- <INPUT> new name
       * 
       * @return 0 if success
       */
      virtual int set_col_name (const t_long col, const std::string &name)
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
      virtual std::string get_col_name (const t_long col) const
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
      virtual t_double *get_col_ptr (const t_long col)
        {
          return &(values[col])[0];
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
      virtual t_long add_col_vector (const std::string &name)
        {
          if (values.empty ())
            return -1;
          t_long n_rows = (t_long) values.front ().size ();
          t_long size = (t_long) values.size ();
          
          values.resize (size + 1); 
          col_names.resize (size + 1);
          
          values[size].resize (n_rows);
          set_col_name (size, name);
          return size;
        }
        
      /** 
       * @brief get numver of rows in table
       */
      virtual t_long get_n_rows () const
        {
          return (t_long) values[0].size ();
        }
      
      /** 
       * @brief get number of columns in table
       */
      virtual t_long get_n_cols () const
        {
          return (t_long) values.size ();
        }

      /** 
       * @brief clear table data 
       */
      virtual void clear ()
        {
          values.clear ();
          col_names.clear ();
        }

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
      
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      std::vector<std::string> col_names;
      std::vector<std::vector<t_double> > values;

      BLUE_SKY_TYPE_DECL (table);
    };

}; //end of blue_sky namespace


#endif /* end of include guard: TABLE_X6PKPWC9 */

