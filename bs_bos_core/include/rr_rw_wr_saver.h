/**
 *       \file  rr_rw_wr_saver.h
 *      \brief  Helpers for save rr, rw and wr arrays from connections
 *              as one plain file
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_RR_RW_WR_SAVER_H_
#define BS_RR_RW_WR_SAVER_H_

#include "facility_manager.h"


namespace blue_sky
  {
  namespace tools
    {

    template <typename strategy_t, typename array_type, typename accessor_t>
    struct rr_rw_wr_saver
      {
        typedef typename strategy_t::item_t                   item_t;
        typedef typename strategy_t::index_t                  index_t;
        typedef facility_manager <strategy_t>                 facility_manager_t;
        typedef wells::connection <strategy_t>                connection_t;
        typedef typename facility_manager_t::well_iterator_t  well_iterator_t;
        typedef item_t                                        value_type;
        typedef array_type                                    array_t;

        typedef smart_ptr <connection_t, true>                sp_connection_t;

        rr_rw_wr_saver (well_iterator_t wb, well_iterator_t we, size_t connection_count)
            : wb_ (wb)
            , we_ (we)
            , array_size_ (array_t::static_size)
            , total_connection_count_ (connection_count)
            , accessor_ (accessor_t ())
            , current_connection_ (0)
            , current_connection_count_ (wb == we ? 0 : (*wb)->get_connection_list ().size ())
            , current_item_ (0)
            , total_item_ (0)
        {
        }

        size_t
        size () const
          {
            return array_size_ * total_connection_count_;
          }

        value_type
        operator [] (size_t i) const
          {
            if (current_item_ >= array_size_)
              {
                current_item_ = 0;
                current_connection_ ++;
              }

            if (current_connection_ >= current_connection_count_)
              {
                if (wb_ != we_)
                  {
                    ++wb_;
                    if (wb_ != we_)
                      {
                        current_connection_ = 0;
                        current_connection_count_ = (*wb_)->get_connection_list ().size ();
                      }
                  }
              }

            total_item_++;
            BS_ASSERT (total_item_ <= size ()) (total_item_) (size ());
            return accessor_ ((*wb_)->get_connection_list ()[current_connection_], current_item_++);
          }

        mutable well_iterator_t   wb_;
        well_iterator_t           we_;
        size_t                    array_size_;
        size_t                    total_connection_count_;
        accessor_t                accessor_;

        mutable size_t            current_connection_;
        mutable size_t            current_connection_count_;
        mutable size_t            current_item_;
        mutable size_t            total_item_;
      };

    template <typename strategy_t, typename accessor_t>
    struct connection_member_saver
      {
        typedef typename strategy_t::item_t     item_t;
        typedef typename strategy_t::index_t    index_t;

        typedef facility_manager <strategy_t>   facility_manager_t;
        typedef wells::connection <strategy_t>  connection_t;
        typedef typename facility_manager_t::well_iterator_t  well_iterator_t;

        typedef item_t                          value_type;

        typedef smart_ptr <connection_t, true>  sp_connection_t;

        connection_member_saver (well_iterator_t wb, well_iterator_t we, size_t array_size)
            : wb_ (wb)
            , we_ (we)
            , array_size_ (array_size)
            , accessor_ (accessor_t ())
            , current_connection_ (0)
            , current_connection_count_ (wb == we ? 0 : (*wb)->get_connection_list ().size ())
            , last_value_ (0)
        {
        }

        size_t
        size () const
          {
            return array_size_;
          }

        value_type
        operator[] (size_t i) const
          {
            if (wb_ == we_)
              {
                BS_ASSERT (wb_ != we_);
                return last_value_;
              }

            if (current_connection_ >= current_connection_count_)
              {
                ++wb_;
                if (wb_ != we_)
                  {
                    current_connection_ = 0;
                    current_connection_count_ = (*wb_)->get_connection_list ().size ();
                  }
              }

            if (wb_ == we_)
              {
                BS_ASSERT (wb_ != we_);
                return last_value_;
              }

            last_value_ = accessor_ ((*wb_)->get_connection_list ()[current_connection_++]);
            return last_value_;
          }

        mutable well_iterator_t       wb_;
        well_iterator_t               we_;
        size_t                        array_size_;
        accessor_t                    accessor_;

        mutable size_t                current_connection_;
        mutable size_t                current_connection_count_;
        mutable item_t                last_value_;
      };

    template <typename strategy_t, typename function_t>
    struct well_member_saver
      {
        typedef typename strategy_t::item_t                   item_t;
        typedef facility_manager <strategy_t>                 facility_manager_t;
        typedef typename facility_manager_t::well_iterator_t  well_iterator_t;

        typedef item_t                                        value_type;
        typedef well <strategy_t>                             well_t;

        //typedef item_t (well_t::*accessor_t) () const;

        typedef function_t                                    accessor_t;

        well_member_saver (well_iterator_t wb, well_iterator_t we)
            : wb_ (wb)
            , we_ (we)
            , accessor_ (accessor_t ())
            , last_value_ (0)
        {
        }

        size_t
        size () const
          {
            return std::distance (wb_, we_);
          }

        value_type
        operator[] (size_t i) const
          {
            if (wb_ == we_)
              {
                BS_ASSERT (wb_ != we_);
                return last_value_;
              }

            last_value_ = accessor_ (*wb_);
            ++wb_;
            return last_value_;
          }

        mutable well_iterator_t   wb_;
        well_iterator_t           we_;
        accessor_t                accessor_;
        mutable item_t            last_value_;
      };

  } // namespace tools
} // namespace blue_sky



#endif  // #ifndef BS_RR_RW_WR_SAVER_H_

