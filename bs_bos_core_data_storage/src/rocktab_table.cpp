/**
 * @file bs_table.cpp
 * @brief table implementation
 * @author
 * @date 2007-09-28
 */
#include "bs_bos_core_data_storage_stdafx.h"
#include "rocktab_table.h"

namespace blue_sky
  {

  /**
   * @brief constructor
   */
  bs_table::bs_table()
  {
    n_rows = 0;
  }

  /**
   * @brief destructor
   */
  bs_table::~bs_table()
  {
    clear();
  }

  /**
   * @brief get number of rows
   *
   * @return number of rows
   */
  strategy_t::index_t
  bs_table::get_num_rows ()
  {
    return n_rows;
  }

  /**
   * @brief get number of columns
   *
   * @return number of columns
   */
  strategy_t::index_t
  bs_table::get_num_cols ()
  {
    return (int)columns.size ();
  }

  /**
   * @brief set number of rows
   *
   * @param num_rows -- new number of rows
   */
  void
  bs_table::set_num_rows (const index_t num_rows)
  {
    std::vector<column>::iterator itr, end_itr = columns.end ();

    for (itr = columns.begin (); itr != end_itr; ++itr)
      {
        itr->data.resize(num_rows);
      }
    n_rows  = num_rows;
  }

  /**
   * @brief add columnt with given name
   *
   * @param _name -- name of new column
   */
  void
  bs_table::add_column(const std::string &_name)
  {
    column last_column;

    columns.push_back (last_column);
    columns.back ().name = _name;
    columns.back ().data.resize (n_rows);
  }

  /**
   * @brief clear
   */
  void
  bs_table::clear()
  {
    columns.clear();
  }

  /**
   * @brief get column name
   *
   * @param col_num -- column index
   *
   * @return column name
   */
  const std::string &
  bs_table::get_col_name(const index_t col_num)
  {
    return columns[col_num].name;
  }

  /**
   * @brief set column name
   *
   * @param col_num -- column index
   * @param new_name -- name
   */
  void
  bs_table::set_col_name (const index_t col_num, const std::string &new_name)
  {
    columns[col_num].name = new_name;
  }

  /**
   * @brief get column data
   *
   * @param col_num -- column index
   *
   * @return column data
   */
  std::vector<strategy_t::item_t> &
  bs_table::get_column (const index_t col_num)
  {
    return columns[col_num].data;
  }

  /**
   * @brief constructor
   */
  rocktab_table::rocktab_table ()
  {
    this->add_column ("pressure");
    this->add_column ("poro volume multiplier");
    this->add_column ("transmissibility multiplier");
  }

  rocktab_table & 
	  rocktab_table::operator= (const rocktab_table & rhs)
  {
	  if (this == &rhs)
		  return *this;

	  bs_table::operator=(rhs);
	  return *this; 
  }

  rocktab_table::rocktab_table (const rocktab_table & rhs)
  {
	  *this = rhs;
  }

  /**
   * @brief destructor
   */
  rocktab_table::~rocktab_table ()
  {
  }

} // namespace blue_sky
