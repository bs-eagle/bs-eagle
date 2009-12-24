/**
 * */
#ifndef BS_EAGLE_WELL_TOOLS_H_
#define BS_EAGLE_WELL_TOOLS_H_

namespace blue_sky {
namespace wells {
namespace detail {

  template <typename strategy_t, typename well_t, typename connection_t>
  void
  fill_rows (const well_t *well, typename strategy_t::index_array_t &rows)
  {
    typedef default_connection_iterator_impl <strategy_t, well_t, connection_t> iterator_t;
    iterator_t it (well, begin_iterator_tag), e (well, end_iterator_tag);
    for (; !well->is_shut () && it != e; ++it)
      {
        if (!it->is_shut ())
          {
            typename strategy_t::index_t n_block = it->n_block ();
            typename strategy_t::index_t k = rows[n_block + 1] != 0;
            rows[n_block + 1] += well->open_connections_count_ - k;
          }
      }
  }


} // namespace detail
} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_WELL_TOOLS_H_

