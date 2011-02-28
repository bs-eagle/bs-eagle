/** 
 * @file table_iface.h
 * @brief Table interface class
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-25
 */
#ifndef TABLE_IFACE_EUVNSA26

#define TABLE_IFACE_EUVNSA26
#include <string>

#include "bs_object_base.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "conf.h"
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{
class table_iface : public objbase
  {
    public:
      /** 
       * @brief destructor
       */
      virtual ~table_iface ()
        {}

      /** 
       * @brief Initialize or reinitialize table by <n_rows> <n_cols> 
       * 
       * @param n_rows -- <INPUT> number of rows in table
       * @param n_cols -- <INPUT> number of cols in table
       * 
       * @return 0 if success
       */
      virtual int init (const t_long n_rows, const t_long n_cols) = 0;

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
      virtual t_double *get_col_ptr (const t_long col) = 0;

      /** 
       * @brief get numver of rows in table
       */
      virtual t_long get_n_rows () const = 0;
      
      /** 
       * @brief get number of columns in table
       */
      virtual t_long get_n_cols () const = 0;

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
};

}  // end of bluesky name space

#endif /* end of include guard: TABLE_IFACE_EUVNSA26 */
