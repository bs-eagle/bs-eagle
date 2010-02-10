#ifndef ROCKTAB_TABLE_H
#define ROCKTAB_TABLE_H

#include "interpolation_macro.h"

namespace blue_sky
  {
  /**
   * @brief table class
   */
  template <class strategy_t>
  class BS_API_PLUGIN bs_table
    {

      typedef typename strategy_t::item_t   item_t;
      typedef typename strategy_t::index_t  index_t;
      //-------------------------
      // METHODS
      //-------------------------
    public:
      // constructor
      bs_table ();
      // Destructor
      ~bs_table ();
      // returns number of rows
      index_t get_num_rows();
      // returns number of columns
      index_t get_num_cols();

      // set number of rows in table
      void set_num_rows (const index_t num_rows);
      // adds new column to the table
      void add_column (const std::string &_name);
      // erases all the elements of a list
      void clear();

      // returns name of the column with given index
      const std::string &get_col_name (const index_t col_num);

      //returns name of current column
      void set_col_name(const index_t col_num, const std::string &new_name);

      //returns current data array from column
      std::vector<item_t> &get_column(const index_t col_num);



      //-------------------------
      // VARIABLES
      //-------------------------
    protected:
      struct column
        {
          std::string name;               //!< column name
          std::vector<item_t> data;       //!< data array stored in one column
        };

      std::vector<column> columns;      //!< list of columns
      index_t n_rows;                       //!< number of rows
    };

  template <class strategy_t>
  class BS_API_PLUGIN rocktab_table  : public bs_table <strategy_t>
    {
    public:
      typedef typename strategy_t::index_t  index_t;
      typedef typename strategy_t::item_t   item_t;

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
        std::vector<item_t> &p = this->get_column (0);
        std::vector<item_t> &phi = this->get_column (1);
        std::vector<item_t> &truns = this->get_column (2);

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
            index_t iu, im, il;
            item_t denom;

            BINARY_SEARCH (p, pressure, this->n_rows, iu, im, il);
            im = (il - 1);
            iu = il;
            denom = (item_t) (1.0 / (p[iu] - p[im]));
            *dp_phi_m = (phi[iu] - phi[im]) * denom;
            *dp_truns_m = (truns[iu] - truns[im]) * denom;
            *phi_m = phi[im] + *dp_phi_m * (pressure - p[im]);
            *truns_m = truns[im] + *dp_truns_m * (pressure - p[im]);
          }
      }
    };
}

#endif // ROCKTAB_TABLE_H
