/**
 * */
#ifndef BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_
#define BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_

#include "connection_iterator.h"
#include "default_well.h"

#include <boost/shared_ptr.hpp>

namespace blue_sky {
namespace wells {

  enum iterator_tag
  {
    begin_iterator_tag,
    end_iterator_tag,
  };

  template <typename well_t, typename connection_t>
  struct default_connection_iterator_impl 
  {
    typedef default_connection_iterator_impl <well_t, connection_t>     this_t;

    typedef smart_ptr <connection_t>                  sp_connection_t;
    typedef smart_ptr <well_t>                        sp_well_t;

  public:

    default_connection_iterator_impl (const well_t *well, 
      iterator_tag tag)
    : well_ (const_cast <well_t *> (well))
    , list_idx_ (0)
    {
      if (tag == begin_iterator_tag)
        {
          connection_idx_[0] = 0;
          connection_idx_[1] = 0;
        }
      else
        {
          connection_idx_[0] = well_->primary_connection_list_.size ();
          connection_idx_[1] = well_->secondary_connection_list_.size ();
        }

      connection_list_[0] = &well_->primary_connection_list_;
      connection_list_[1] = &well_->secondary_connection_list_;

#ifdef _DEBUG
      connection_count_[0] = well_->primary_connection_list_.size ();
      connection_count_[1] = well_->secondary_connection_list_.size ();
#endif
    }

    inline sp_connection_t
    operator* () const
    {
#ifdef _DEBUG
      if (static_cast <size_t> (connection_count_[list_idx_]) != connection_list_[list_idx_]->size ())
        {
          bs_throw_exception (boost::format ("Iterator no more valid (list: %ld, well: %s)")
            % list_idx_
            % well_->name ());
        }
      if (static_cast <size_t> (connection_idx_[list_idx_]) >= connection_list_[list_idx_]->size ())
        {
          bs_throw_exception (boost::format ("Index out of range (list: %ls, idx: %ld, count: %ld, well: %s)")
            % list_idx_
            % connection_idx_[list_idx_] % connection_count_[list_idx_] % well_->name ())
        }
#endif

      return connection_list_[list_idx_]->operator[] (connection_idx_[list_idx_]);
    }

    inline sp_connection_t
    operator-> () const
    {
      return operator* ();
    }

    inline this_t &
    operator++ ()
    {
      if (list_idx_ > 1)
        {
          //bs_throw_exception (boost::format ("List index out of range (list: %ld, well: %s)")
          //  % list_idx_ % well_->name ());

          return *this;
        }

#ifdef _DEBUG
      if (static_cast <size_t> (connection_count_[list_idx_]) != connection_list_[list_idx_]->size ())
        {
          bs_throw_exception (boost::format ("Iterator no more valid (list: %ld, well: %s)")
            % list_idx_
            % well_->name ());
        }
      if (static_cast <size_t> (connection_idx_[list_idx_]) >= connection_list_[list_idx_]->size ())
        {
          bs_throw_exception (boost::format ("Index out of range (list: %ls, idx: %ld, count: %ld, well: %s)")
            % list_idx_ % connection_idx_[list_idx_] % connection_count_[list_idx_] % well_->name ())
        }
#endif
      ++connection_idx_[list_idx_];
      list_idx_ += static_cast <size_t> (connection_idx_[list_idx_]) >= connection_list_[list_idx_]->size ();

      return *this;
    }

    inline bool
    operator== (const this_t &rhs) const
    {
      return position () == rhs.position ();
    }
    inline bool
    operator!= (const this_t &rhs) const
    {
      return !operator== (rhs);
    }

    typename strategy_t::index_t
    position () const
    {
      return connection_idx_[0] + connection_idx_[1];
    }

  private:
    well_t                                      *well_;
    typename strategy_t::index_t                list_idx_;
    typename strategy_t::index_t                connection_idx_[2];
    typename well_t::connection_list_t          *connection_list_[2];

#ifdef _DEBUG
    typename strategy_t::index_t                connection_count_[2];
#endif
  };

  template <typename well_t, typename connection_t>
  struct default_connection_iterator : connection_iterator::impl
  {
    typedef connection_iterator::impl   base_t;
    typedef connection base_connection_t;
    typedef smart_ptr <base_connection_t>                     sp_base_connection_t;

  public:

    default_connection_iterator (const well_t *well, 
      iterator_tag tag)
    : impl_ (well, tag)
    {
    }

    virtual sp_base_connection_t
    operator* () const
    {
      return impl_.operator* ();
    }

    virtual sp_base_connection_t
    operator-> () const
    {
      return impl_.operator-> ();
    }

    virtual base_t &
    operator++ ()
    {
      impl_.operator++ ();
      return *this;
    }

    virtual bool
    operator== (const boost::shared_ptr <base_t> &rhs) const
    {
      return impl_.position () == rhs->position ();
    }
    virtual bool
    operator!= (const boost::shared_ptr <base_t> &rhs) const
    {
      return !operator== (rhs);
    }

    virtual typename strategy_t::index_t
    position () const
    {
      return impl_.position ();
    }

  private:
    default_connection_iterator_impl <well_t, connection_t> impl_;
  };


} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_EAGLE_DEFAULT_CONNECTION_ITERATOR_H_

