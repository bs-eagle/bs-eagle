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
  template <class strategy_t>
  bs_table<strategy_t>::bs_table()
  {
    n_rows = 0;
  }

  /**
   * @brief destructor
   */
  template <class strategy_t>
  bs_table <strategy_t>::~bs_table()
  {
    clear();
  }

  /**
   * @brief get number of rows
   *
   * @return number of rows
   */
  template <class strategy_t>
  typename strategy_t::index_t
  bs_table <strategy_t> ::get_num_rows ()
  {
    return n_rows;
  }

  /**
   * @brief get number of columns
   *
   * @return number of columns
   */
  template <class strategy_t>
  typename strategy_t::index_t
  bs_table<strategy_t>::get_num_cols ()
  {
    return (int)columns.size ();
  }

  /**
   * @brief set number of rows
   *
   * @param num_rows -- new number of rows
   */
  template <class strategy_t>
  void
  bs_table<strategy_t>::set_num_rows (const index_t num_rows)
  {
    typename std::vector<column>::iterator itr, end_itr = columns.end ();

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
  template <class strategy_t>
  void
  bs_table<strategy_t>::add_column(const std::string &_name)
  {
    column last_column;

    columns.push_back (last_column);
    columns.back ().name = _name;
    columns.back ().data.resize (n_rows);
  }

  /**
   * @brief clear
   */
  template <class strategy_t>
  void
  bs_table<strategy_t>::clear()
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
  template <class strategy_t>
  const std::string &
  bs_table <strategy_t>::get_col_name(const index_t col_num)
  {
    return columns[col_num].name;
  }

  /**
   * @brief set column name
   *
   * @param col_num -- column index
   * @param new_name -- name
   */
  template <class strategy_t>
  void
  bs_table <strategy_t>::set_col_name (const index_t col_num, const std::string &new_name)
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
  template <class strategy_t>
  std::vector<typename strategy_t::item_t> &
  bs_table <strategy_t>::get_column (const index_t col_num)
  {
    return columns[col_num].data;
  }




  /**
   * @brief constructor
   */
  template <class strategy_t>
  rocktab_table <strategy_t>::rocktab_table ()
  {
    this->add_column ("pressure");
    this->add_column ("poro volume multiplier");
    this->add_column ("transmissibility multiplier");
  }

  /**
   * @brief destructor
   */
  template <class strategy_t>
  rocktab_table <strategy_t>::~rocktab_table ()
  {
  }

  template class bs_table <base_strategy_di>;
  template class bs_table <base_strategy_fi>;
  template class bs_table <base_strategy_mixi>;

  template class rocktab_table <base_strategy_di>;
  template class rocktab_table <base_strategy_fi>;
  template class rocktab_table <base_strategy_mixi>;

} // namespace blue_sky
