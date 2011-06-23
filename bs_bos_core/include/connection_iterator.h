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

  struct connection_iterator : 
    std::iterator <
      std::forward_iterator_tag, 
      smart_ptr <wells::connection, true>
    >
  {
    typedef connection_iterator this_t;
    typedef wells::connection connection_t;

    typedef smart_ptr <connection_t, true> sp_connection_t;

  public:

    /**
     * \brief  Returns object which iterator points
     * \return Base connection object
     * */
    sp_connection_t
    operator* () const
    {
      BS_ASSERT (impl_);
      return impl_->operator * ();
    }

    /**
     * \brief  Returns object which iterator points
     * \return Base connection object
     * */
    sp_connection_t
    operator-> () const
    {
      BS_ASSERT (impl_);
      return impl_->operator * ();
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
      impl_->operator++ ();
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
      return impl_->operator== (rhs.impl_);
    }

    /**
     * \brief  Checks two iterators
     * \return True if two iterators are not equal
     * */
    bool
    operator!= (const this_t &rhs) const
    {
      BS_ASSERT (impl_);
      return !impl_->operator== (rhs.impl_);
    }

    /**
     * \brief  Returns position of iterator, useful
     *         for some auxulliary functions
     * \return Iterator position
     * */
    typename strategy_t::index_t
    position () const
    {
      BS_ASSERT (impl_);
      return impl_->position ();
    }

    struct impl
    {
      virtual ~impl () {}

      virtual sp_connection_t
      operator* () const = 0;

      virtual sp_connection_t
      operator-> () const = 0;

      virtual impl &
      operator++ () = 0;

      virtual bool
      operator== (const boost::shared_ptr <impl> &rhs) const = 0;

      virtual bool
      operator!= (const boost::shared_ptr <impl> &rhs) const = 0;

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

