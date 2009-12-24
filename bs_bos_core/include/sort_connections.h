/**
 * */
#ifndef BS_EAGLE_SORT_CONNECTIONS_H_
#define BS_EAGLE_SORT_CONNECTIONS_H_

namespace blue_sky {
namespace wells {
namespace detail {

  /**
   * \class sort_connection_list
   * \brief Sorts connection list in ASC order
   * */
  template <typename connection_t>
  struct sort_connection_list : std::binary_function <const connection_t &, const connection_t &, bool>
  {
    /**
     * \brief  Sorts connection list in ASC order
     * \param  lhs
     * \param  rhs
     * \return True if n_block of lhs < of n_block of rhs
     * */
    bool
    operator () (const connection_t &lhs, const connection_t &rhs) const
    {
      return lhs->n_block () < rhs->n_block ();
    }
  };

} // namespace detail
} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_SORT_CONNECTIONS_H_

