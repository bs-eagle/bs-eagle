#ifndef ROCKTAB_TABLE_H
#define ROCKTAB_TABLE_H

#include "interpolation_macro.h"
#include "conf.h"

#include "bs_serialize_decl.h"

namespace blue_sky
  {
  /**
   * @brief table class
   */
  
  class BS_API_PLUGIN bs_table
    {

      //-------------------------
      // METHODS
      //-------------------------
    public:
      // constructor
      bs_table ();
      // Destructor
      ~bs_table ();
      // returns number of rows
      t_int get_num_rows();
      // returns number of columns
      t_int get_num_cols();

      // set number of rows in table
      void set_num_rows (const t_int num_rows);
      // adds new column to the table
      void add_column (const std::string &_name);
      // erases all the elements of a list
      void clear();

      // returns name of the column with given index
      const std::wstring &get_col_name (const t_int col_num);

      //returns name of current column
      void set_col_name(const t_int col_num, const std::wstring &new_name);

      //returns current data array from column
      std::vector<t_float> &get_column(const t_int col_num);



      //-------------------------
      // VARIABLES
      //-------------------------
    protected:
      struct column
        {
          std::wstring name;               //!< column name
          std::vector<t_float> data;       //!< data array stored in one column
        };

      std::vector<column> columns;      //!< list of columns
      t_int n_rows;                       //!< number of rows

      friend class blue_sky::bs_serialize;
    };

  
  class BS_API_PLUGIN rocktab_table  : public bs_table 
    {
    public:

      //! ctor
      rocktab_table();
      //! copy ctor
      //rocktab_table(const rocktab_table&);
      //! dtor
      ~rocktab_table();

	  rocktab_table (const rocktab_table &rhs);
	  rocktab_table & operator = (const rocktab_table & rhs);

      //! interpolate method

      //! \brief interpolate values in rocktab table
      template <typename T>
      inline 
      void interpolate(T pressure,    //!< pressure value 
                       T *phi_m,      //!< porosity multiplier                     
                       T *dp_phi_m,   //!< pressure derivation of dp_phi_m
                       T *truns_m,    //!< transmissibility multiplier
                       T *dp_truns_m  //!< pressure derivation of transmissibility
                       )
      {
        std::vector<t_float> &p = this->get_column (0);
        std::vector<t_float> &phi = this->get_column (1);
        std::vector<t_float> &truns = this->get_column (2);

        if (this->n_rows < 1)
          return;

        if (pressure <= p[0])
          {
            *phi_m = phi[0];
            *dp_phi_m = 0.0;
            *truns_m = truns[0];
            *dp_truns_m = 0.0;
          }
        else if (pressure >= p[this->n_rows - 1])
          {
            *phi_m = phi[this->n_rows - 1];
            *dp_phi_m = 0.0;
            *truns_m = truns[this->n_rows - 1];
            *dp_truns_m = 0.0;
          }
        else
          {
            t_int iu, im, il;
            t_double denom;

            BINARY_SEARCH (p, pressure, this->n_rows, iu, im, il);
            im = (il - 1);
            iu = il;
            denom = (t_double) (1.0 / (p[iu] - p[im]));
            *dp_phi_m = (phi[iu] - phi[im]) * denom;
            *dp_truns_m = (truns[iu] - truns[im]) * denom;
            *phi_m = phi[im] + *dp_phi_m * (pressure - p[im]);
            *truns_m = truns[im] + *dp_truns_m * (pressure - p[im]);
          }
      }
    };
}

#endif // ROCKTAB_TABLE_H
