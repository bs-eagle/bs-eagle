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

#include "bs_serialize_decl.h"
#include "table_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN table : public table_iface
    {
    
    public: 

      typedef std::vector <t_double>                  vector_t;
      typedef std::vector <vector_t>                  table_t;
      typedef BS_SP (table_iface)                     sp_table_t;

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
      /*!
        \brief copy 
      */
      virtual int copy (const sp_table_t a);

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
      virtual int set_col_name (const t_long col, const std::wstring &name)
        {
          col_names[col] = name;
          return 0;
        }

      /** 
       * @brief remove row from bs.comm.table
       * 
       * @param row_index -- row index
       */
      virtual void remove_row (const t_long row_index);
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
       * @brief return col names
       */
      virtual std::vector<std::wstring> get_col_names () const
        {
          return col_names;
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
      virtual t_long add_col_vector (const std::wstring &name)
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

      virtual void 
      add_col_vector (t_long col, const std::wstring &name, spv_double new_vector)
        {
          t_double *new_vector_ = new_vector->data ();
          t_long size = values.size ();
          BS_ASSERT (col >= 0 && col < size);
          
          t_long new_vec_size = new_vector->size ();
          values[col].resize (new_vec_size);
          for (t_long i = 0; i < new_vec_size; ++i)
            {
              values[col][i] = new_vector_[i];
            } 
          set_col_name (col, name);
        }

        
      /** 
       * @brief get numver of rows in table
       */
      virtual t_long get_n_rows (t_long col = 0) const
        {
          if (values.size ())
            return (t_long) values[col].size ();
          else
            return 0;
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

      /** 
       * @brief add new row to the table at index #row_index
       * 
       * @param row_index -- <INPUT> given index of new row
       *                        if 0 insert row at first
       *                        if > n_rows add row to the end
       */
      virtual void add_row (const t_long row_index);

      /** 
       * @brief return value at the given #row and #col
       * 
       * @param row     -- <INPUT> given row
       * @param col     -- <INPUT> given column
       * 
       */
      virtual t_double get_value (const t_long row, const t_long col) const
        {
          return (values[col])[row];
        }

      /** 
       * @brief return values 
       */
      virtual table_t get_values () const 
        {
          return values;
        }


      /** 
       * @brief set new value at given #row and #col
       * 
       * @param row     -- <INPUT> given row
       * @param col     -- <INPUT> given column
       * @param val     -- <INPUT> given value
       */
      virtual void set_value (const t_long row, const t_long col, const t_double val)
        {
          (values[col])[row] = val;
        }
      
      /*!
        \brief convert table from array 
      */
      virtual void convert_from_array (const t_long n_rows, const t_long n_cols, spv_double data);

      virtual void convert_from_tf_array (const t_long n_rows, const t_long n_cols, t_float *data);
      
      /*!
        \brief convert table to array 
      */
      virtual spv_double convert_to_array (const t_long n_rows, const t_long n_cols) const;

      virtual spv_float convert_to_array_f(const t_long n_rows, const t_long n_cols) const;
      /** 
       * @brief push back new row
       * 
       * @param v  -- <INPUT> row values
       */
      virtual void push_back (const std::vector<t_double> &v)
        {
          if (v.size () != values.size ())
            return;
          t_long i, n;
          n = v.size ();
          for (i = 0; i < n; ++i)
            {
              values[i].push_back (v[i]);
            }
        }
      virtual void push_back_f (const std::vector<t_float> &v) {
        this->push_back(std::vector< t_double >(v.begin(), v.end()));
      }

      virtual sp_table_t check_serial () const;

      /** 
       * @brief pack(serialize) all information of class to text string 
       * 
       * @return string
       */
      virtual std::string to_str () const; 

      /** 
       * @brief Reastore all class information from input text string
       * 
       * @param s -- <INPUT> string
       */
      virtual void from_str (const std::string &s);

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;
#endif //BSPY_EXPORTING_PLUGIN

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
      
      void remove_col(const t_long col);
      void remove_col(const std::wstring& col_name);
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      std::vector<std::wstring> col_names;
      table_t values;

      BLUE_SKY_TYPE_DECL (table);
      friend class bs_serialize;
    };

}; //end of blue_sky namespace


#endif /* end of include guard: TABLE_X6PKPWC9 */

