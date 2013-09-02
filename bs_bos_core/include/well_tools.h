/**
 * */
#ifndef BS_EAGLE_WELL_TOOLS_H_
#define BS_EAGLE_WELL_TOOLS_H_

namespace blue_sky {
namespace wells {
namespace detail {

  template <typename well_t, typename connection_t>
  void
  fill_rows (const well_t *well, v_long &rows)
  {
    typedef default_connection_iterator_impl <well_t, connection_t> iterator_t;
    iterator_t it (well, begin_iterator_tag), e (well, end_iterator_tag);
    for (; !well->is_shut () && it != e; ++it)
      {
        if (!it->is_shut ())
          {
            t_long n_block = it->n_block ();
            t_long k = rows[n_block + 1] != 0;
            rows[n_block + 1] += well->open_connections_count_ - k;
          }
      }
  }

  // FIXME: dereference
  template <typename well_t, typename connection_t>
  void
  fill_jacobian (const well_t                 *well, 
    double                                    dt, 
    t_long                                    block_size,
    const spv_long                            &rows,
    spv_long                                  &cols,
    spv_float                                 &values, 
    stdv_long                                 &markers)
  {
    typedef default_connection_iterator_impl <well_t, connection_t> iterator_t;
    iterator_t it (well, begin_iterator_tag), e (well, end_iterator_tag);

    t_long b_sqr = block_size * block_size;
    for (; !well->is_shut () && it != e; ++it)
      {
        const smart_ptr <connection_t> &rw_con = *it;
        if (!rw_con->is_shut ())
          {
            t_long n_block = rw_con->n_block ();
            t_long l = (*rows)[n_block];

            if (markers[n_block] == 0)
              {
                markers[n_block] = 1;
              }

            iterator_t wr_it (well, begin_iterator_tag);
            for (; wr_it != e; ++wr_it)
              {
                t_long index = l;
                const smart_ptr <connection_t> &wr_con = *wr_it;
                if (!wr_con->is_shut ())
                  {
                    if (rw_con->n_block () == wr_con->n_block ())
                      {
                        BS_ASSERT ((*cols)[l] == -1) ((*cols)[l]) (n_block);
                        (*cols)[l] = n_block;

                        BS_ASSERT (l * b_sqr < (t_long)values->size ()) (l) (values->size ());
                        eliminate (well, &(*values)[l * b_sqr], rw_con, wr_con, dt, block_size);
                      }
                    else
                      {
                        BS_ASSERT ((*cols)[l + markers[n_block]] == -1) ((*cols)[l + markers[n_block]]) (wr_con->n_block ());
                        index = l + markers[n_block];
                        (*cols)[index] = wr_con->n_block ();
                        markers[n_block]++;

                        BS_ASSERT (index * b_sqr < (t_long)values->size ()) (index) (b_sqr) (values->size ());
                        eliminate (well, &(*values)[index * b_sqr], rw_con, wr_con, dt, block_size);
                      }
                  }
              }
          }
      }
  }


} // namespace detail
} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_WELL_TOOLS_H_

