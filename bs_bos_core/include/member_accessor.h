/**
 * \file member_accessor.h
 * \brief helper to adapt vector of calc_model_data to do save via tools::save_seq_vector (only for debug purpose)
 * \author Sergey Miryanov
 * \date 31.10.2008
 * */
#ifndef BS_TOOLS_MEMBER_ACCESSOR_H_
#define BS_TOOLS_MEMBER_ACCESSOR_H_

namespace blue_sky
  {
  namespace tools
    {

    template <typename data_array_t, typename item_t>
    struct member_accessor
      {
        typedef item_t value_type;

        member_accessor (const data_array_t &data_, const item_t *data_begin_, const item_t *data_offset_, size_t size_)
            : data_ (data_)
            , data_begin_ (data_begin_)
            , data_offset_ (data_offset_)
            , size_ (size_)
        {
        }

        size_t
        size_i () const
          {
            return data_.size ();
          }

        size_t
        size_j () const
          {
            return size_;
          }

        item_t
        get (size_t i, size_t j) const
          {
            const item_t *data = (const item_t *)&data_[i];
            size_t count = data_offset_ - data_begin_;
            return data [count + j];
          }

private:
        const data_array_t  &data_;
        const item_t        *data_begin_;
        const item_t        *data_offset_;
        size_t              size_;
      };


// use like
//
//#define SAVE_BOOST_ARRAY(name)																															\
//	tools::save_seq_vector (BOOST_PP_CAT (BOOST_PP_STRINGIZE(name), ".bs.txt"))								\
//		.save_via_fn (member_accessor (data, (item_t *)&data[0], &data[0].name[0], data[0].name.size ()));
//
//#define SAVE_ITEM(name)																																			\
//	tools::save_seq_vector (BOOST_PP_CAT (BOOST_PP_STRINGIZE(name), ".bs.txt"))								\
//		.save_via_fn (member_accessor (data, (item_t *)&data[0], &data[0].name, 1));


  } // namespace tools
} // namespace blue_sky

#endif  // #ifndef BS_TOOLS_MEMBER_ACCESSOR_H_

