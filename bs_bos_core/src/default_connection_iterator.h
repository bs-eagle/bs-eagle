/**
 * */
#ifndef BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_
#define BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_

#include "connection_iterator.h"
#include "default_well.h"

#include <boost/shared_ptr.hpp>

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  struct default_connection_iterator : connection_iterator <strategy_t>::impl
  {
    typedef typename connection_iterator <strategy_t>::impl   base_t;
    typedef connection <strategy_t>                           connection_t;
    typedef default_connection <strategy_t>                   default_connection_t;
    typedef default_well <strategy_t>                         default_well_t;

    typedef smart_ptr <connection_t>                          sp_connection_t;
    typedef smart_ptr <default_connection_t>                  sp_default_connection_t;
    typedef smart_ptr <default_well_t>                        sp_default_well_t;

  public:

    default_connection_iterator (const sp_default_well_t &well, 
      typename strategy_t::index_t connection_idx)
    : well_ (well)
    , connection_idx_ (connection_idx)
    , connection_count_ (well_->get_connections_count ())
    {
    }

    virtual sp_connection_t
    get () const
    {
#ifdef _DEBUG
      if (connection_count_ != well_->get_connections_count ())
        {
          bs_throw_exception (boost::format ("Iterator no more valid (well: %s)") 
            % well_->name ());
        }
      if (connection_idx_ >= connection_count_)
        {
          bs_throw_exception (boost::format ("Index out of range (idx: %ld, count: %ld, well: %s)") 
            % connection_idx_ % connection_count_ % well_->name ())
        }
#endif

      return well_->get_connection (connection_idx_);
    }

    virtual void
    increment () 
    {
#ifdef _DEBUG
      if (connection_count_ != well_->get_connections_count ())
        {
          bs_throw_exception (boost::format ("Iterator no more valid (well: %s)") 
            % well_->name ());
        }
      if (connection_idx_ >= connection_count_)
        {
          bs_throw_exception (boost::format ("Index out of range (idx: %ld, count: %ld, well: %s)") 
            % connection_idx_ % connection_count_ % well_->name ())
        }
#endif
      ++connection_idx_;
    }

    bool
    is_equal (const boost::shared_ptr <base_t> &rhs) const
    {
      return connection_idx_ == rhs->position ();
    }

    typename strategy_t::index_t
    position () const
    {
      return connection_idx_;
    }

  private:
    sp_default_well_t               well_;
    typename strategy_t::index_t    connection_idx_;
    typename strategy_t::index_t    connection_count_;
  };


} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_

