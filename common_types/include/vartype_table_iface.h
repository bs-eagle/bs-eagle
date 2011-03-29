/** 
 * @file vartype_table_iface.h
 * @brief 
 * @author ildar badykov
 * @date 2009-08-27
 */
#ifndef __VARTYPE_TABLE_IFACE_H
#define __VARTYPE_TABLE_IFACE_H
#include <string>

#include "bs_object_base.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{

  /**
  * \brief properties
  */
  template <class var_type_t>
  class vartype_table_iface : public objbase
    {
      public:
        typedef std::vector <var_type_t>                       vector_t;
        typedef bs_array <var_type_t, bs_nparray>              var_type_array_t;    
        typedef smart_ptr <var_type_array_t, true >            sp_var_type_array_t;
 
      public:
        /** 
         * @brief destructor
         */
        virtual ~vartype_table_iface ()
          {}
      /** 
       * @brief Initialize or reinitialize table by <n_cols> 
       * 
       * @param n_cols -- <INPUT> number of cols in table
       * 
       * @return 0 if success
       */
      virtual int init (const t_long n_cols) = 0;

      /** 
       * @brief set <name> of the column <col>
       * 
       * @param col  -- <INPUT> column index
       * @param name -- <INPUT> new name
       * 
       * @return 0 if success
       */
      virtual int set_col_name (const t_long col, const std::string &name) = 0;

      /** 
       * @brief return name of the <col> column
       * 
       * @param col -- <INPUT> column index
       * 
       * @return name
       */
      virtual std::string get_col_name (const t_long col) const = 0;

      /** 
       * @brief return pointer to the column <col> data 
       * 
       * @param col -- <INPUT> column index
       * 
       * @return pointer
       */
      virtual var_type_t *get_col_ptr (const t_long col) = 0;


      /** 
       * @brief return reference to the column <col> data 
       * 
       * @param col -- <INPUT> column index
       * 
       * @return pointer
       */
      virtual vector_t & get_col_vector (const t_long col) = 0;

      /** 
       * @brief add new column and return reference to the column <col> data 
       * 
       * @param name -- <INPUT> new name       
       * 
       * @return pointer
       */
      virtual t_long add_col_vector (const t_long col, const std::string &name, sp_var_type_array_t new_vector) = 0;
      
      /** 
       * @brief get number of rows in table in column col
       */
      virtual t_long get_n_rows (const t_long col) const = 0;
      
      /** 
       * @brief get number of columns in table
       */
      virtual t_long get_n_cols () const = 0;
      
      /** 
       * @brief clear table data 
       */
      virtual void clear () = 0;
#if 0
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

      /** 
       * @brief set column <col> values by numpy array
       * 
       * @param col     -- <INPUT> column index
       * @param values  -- <INPUT> numpy array
       */
      virtual void set_col_values (const t_long col, spv_double values) = 0;

      /** 
       * @brief return numpy array for given column <col>
       * 
       * @param col -- <INPUT> column index
       * 
       * @return numpy array
       */
      virtual spv_double get_col_values (const t_long col) const = 0;
#endif //BSPY_EXPORTING_PLUGIN
#endif // 0
 
    };
}

#endif //__VARTYPE_TABLE_IFACE_H
 