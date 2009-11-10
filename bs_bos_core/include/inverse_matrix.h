///**
// * \file inverse_matrix.h
// * \brief calc invers of matrix
// * \author Sergey Miryanov
// * \date 25.06.2008
// * */
//#ifndef BS_CALC_INVERSE_MATRIX_H_
//#define BS_CALC_INVERSE_MATRIX_H_
//
//namespace blue_sky {
//namespace tools {
//
//    template <class array_t, class index_t>
//    struct inverse_block
//    {
//      typedef typename array_t::value_type item_t;
//
//      inverse_block (array_t &dst, const array_t &src_, index_t bs_)
//        : dst (dst)
//        , block_size (bs_)
//        , b_sqr (bs_ * bs_)
//      {
//        // TODO: !!
//        array_t &src__ = const_cast <array_t&> (src_);
//        src.assign (src__.begin (), src__.end ());
//      }
//
//      void inverse (index_t index)
//      {
//        identity_matrix (index);
//
//        for (index_t i = 0; i < block_size; ++i)
//          {
//            process_lower_row (index, i);
//          }
//        for (index_t i = block_size - 1; i >= 0; --i)
//          {
//            process_upper_row (index, i);
//          }
//        for (index_t i = 0; i < block_size; ++i)
//          {
//            process_diagonal (index, i);
//          }
//      }
//
//    private:
//      void process_lower_row (index_t index, index_t row_idx)
//      {
//        index_t row_idx_ = (index * b_sqr) + (row_idx * block_size);
//
//        item_t mul0 = src[row_idx_];
//        for (index_t i = row_idx + 1; i < block_size; ++i)
//          {
//            index_t idx = index * b_sqr + i * block_size;
//            item_t mul = src[idx];
//            for (index_t j = 0; j < block_size; ++j)
//              {
//                src[idx + j] = src[idx + j] * mul0 - src[row_idx_ + j] * mul;
//                dst[idx + j] = dst[idx + j] * mul0 - dst[row_idx_ + j] * mul;
//              }
//          }
//      }
//
//      void process_upper_row (index_t index, index_t row_index)
//      {
//        index_t row_idx = (index *b_sqr) + (row_index * block_size);
//
//        item_t mul0 = src[row_idx + row_index];
//        for (index_t i = row_index - 1; i >= 0; --i)
//          {
//            index_t idx = index * b_sqr + i * block_size;
//            item_t mul = src[idx + row_index];
//            for (index_t j = block_size - 1; j >= 0; --j)
//              {
//                src[idx + j] = src[idx + j] * mul0 - src[row_idx + j] * mul;
//                dst[idx + j] = dst[idx + j] * mul0 - dst[row_idx + j] * mul;
//              }
//          }
//      }
//
//      void process_diagonal (index_t index, index_t row_index)
//      {
//        index_t row_idx = (index * b_sqr) + (row_index * block_size);
//        item_t div = src[row_idx + row_index];
//
//        for (index_t i = 0; i < block_size; ++i)
//          {
//            dst[row_idx + i] /= div;
//          }
//      }
//
//      void identity_matrix (index_t index)
//      {
//        for (index_t i = 0; i < b_sqr; ++i)
//          {
//            dst[index * b_sqr + i] = 0;
//          }
//        for (index_t i = 0; i < block_size; ++i)
//          {
//            dst[index * b_sqr + i * block_size + i] = 1.0;
//          }
//      }
//
//      array_t   &dst;
//      array_t   src;
//      index_t   block_size;
//      index_t   b_sqr;
//    };
//
//} // namespace tools
//} // namespace blue_sky
//
//
//#endif  // #ifndef BS_CALC_INVERSE_MATRIX_H_
