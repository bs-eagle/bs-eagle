/**
 * */
#ifndef BS_EAGLE_ADD_CONNECTION_TO_LIST_H_
#define BS_EAGLE_ADD_CONNECTION_TO_LIST_H_

#include "sort_connections.h"

namespace blue_sky {
namespace wells {
namespace detail {

  template <typename index_t, typename connection_list_t, typename well_t>
  typename connection_list_t::value_type
  add_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block, connection_list_t &connection_list, well_t *well)
  {
    if (n_block < 0)
      {
        bs_throw_exception (boost::format ("Invalid connection n_block value (ijk: {%d, %d, %d}, n_block: %d, well: %s)") 
          % i_coord % j_coord % k_coord
          % n_block
          % well->name ());
      }

    typename connection_list_t::value_type connection = BS_KERNEL.create_object (connection_list_t::value_type::element_type::bs_type (), true);
    if (!connection)
        bs_throw_exception (boost::format ("Can't create connection (well: %s)") % well->name ());

    connection->set_coord (i_coord, j_coord, k_coord, n_block);
    connection_list.push_back (connection);
    std::sort (connection_list.begin (), connection_list.end (), sort_connection_list <typename connection_list_t::value_type> ());
    //connection_map_.insert (std::make_pair (connection->n_block (), (index_t)primary_connection_list_.size () - 1));

    return connection;
  }

} // namespace detail 
} // namespace wells 
} // namespace blue_sky

#endif // #ifndef BS_EAGLE_ADD_CONNECTION_TO_LIST_H_

