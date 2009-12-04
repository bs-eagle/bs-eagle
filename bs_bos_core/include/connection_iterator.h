/**
 *       \file  connection_iterator.h
 *      \brief  Compound iterator for well connections
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_EAGLE_CONNECTION_ITERATOR_H_
#define BS_EAGLE_CONNECTION_ITERATOR_H_

#include "well_connection.h"
#include <boost/shared_ptr.hpp>

namespace blue_sky {

  template <typename strategy_t>
  struct connection_iterator : 
    std::iterator <
      std::forward_iterator_tag, 
      smart_ptr <wells::connection <strategy_t> > 
    >
  {
    typedef connection_iterator <strategy_t>  this_t;
    typedef wells::connection <strategy_t>    connection_t;

    typedef smart_ptr <connection_t>          sp_connection_t;

  public:

    /**
     * \brief  Returns object which iterator points
     * \return Base connection object
     * */
    sp_connection_t
    operator* () const
    {
      BS_ASSERT (impl_);
      return impl_->get ();
    }

    /**
     * \brief  Returns object which iterator points
     * \return Base connection object
     * */
    sp_connection_t
    operator-> () const
    {
      BS_ASSERT (impl_);
      return impl_->get ();
    }

    /**
     * \brief  Moves iterator forward, only 
     *         preincrement operation supports
     * \return This object
     * */
    this_t &
    operator++ ()
    {
      BS_ASSERT (impl_);
      impl_->increment ();
      return *this;
    }

    /**
     * \brief  Checks two iterators 
     * \return True if two iterators are equal
     * */
    bool
    operator== (const this_t &rhs) const
    {
      BS_ASSERT (impl_);
      return impl_->is_equal (rhs.impl_);
    }

    /**
     * \brief  Checks two iterators
     * \return True if two iterators are not equal
     * */
    bool
    operator!= (const this_t &rhs) const
    {
      BS_ASSERT (impl_);
      return !impl_->is_equal (rhs.impl_);
    }

    struct impl
    {
      virtual ~impl () {}

      virtual sp_connection_t
      get () const = 0;

      virtual void
      increment () = 0;

      virtual bool
      is_equal (const boost::shared_ptr <impl> &rhs) const = 0;

      virtual typename strategy_t::index_t
      position () const = 0;
    };

    connection_iterator (impl *impl_)
    : impl_ (impl_)
    {
    }

  protected:
    boost::shared_ptr <impl> impl_;
  };


} // namespace blue_sky

#endif // #ifndef BS_EAGLE_CONNECTION_ITERATOR_H_

