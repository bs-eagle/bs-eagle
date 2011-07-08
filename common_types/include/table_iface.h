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
      typedef std::vector <t_double>                  vector_t;
      typedef BS_SP (table_iface)                     sp_table_iface;

    public:
      /** 
       * @brief destructor
       */
      virtual ~table_iface ()
        {}
      
      /*!
        \brief copy 
      */
      virtual int copy (const sp_table_iface a) = 0;
      
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
      virtual t_long add_col_vector (const std::string &name) = 0;

      virtual void add_col_vector (t_long col, std::string const &name, spv_double array) = 0;
      
      /** 
       * @brief get numver of rows in table
       */
      virtual t_long get_n_rows (t_long col = 0) const = 0;
      
      /** 
       * @brief get number of columns in table
       */
      virtual t_long get_n_cols () const = 0;
      
      /** 
       * @brief clear table data 
       */
      virtual void clear () = 0;

      /** 
       * @brief add new row to the table at index #row_index
       * 
       * @param row_index -- <INPUT> given index of new row
       *                        if 0 insert row at first
       *                        if > n_rows add row to the end
       */
      virtual void add_row (const t_long row_index) = 0;

      /** 
       * @brief return value at the given #row and #col
       * 
       * @param row     -- <INPUT> given row
       * @param col     -- <INPUT> given column
       * 
       */
      virtual t_double get_value (const t_long row, const t_long col) const = 0;

      /** 
       * @brief set new value at given #row and #col
       * 
       * @param row     -- <INPUT> given row
       * @param col     -- <INPUT> given column
       * @param val     -- <INPUT> given value
       */
      virtual void set_value (const t_long row, const t_long col, const t_double val) = 0;
      
      /*!
        \brief convert table to array 
      */
      virtual spv_double convert_to_array (const t_long n_rows, const t_long n_cols) const = 0;
      

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
