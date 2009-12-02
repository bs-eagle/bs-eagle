/**
 * \file scale_arrays_placement_strategies.h
 * \brief strategies impl for placement of scale arrays
 * \author Sergey Miryanov
 * \date 25.05.2008
 * */
#ifndef BS_SCALE_ARRAYS_PLACEMENT_STRATEGIES_H_
#define BS_SCALE_ARRAYS_PLACEMENT_STRATEGIES_H_

#include "scal_data_placement_info.h"

namespace blue_sky
  {
  namespace scal
    {
    namespace data_placement
      {

      struct separate_scale_vectors_t
        {
          typedef shared_vector <float> vector_t;

          enum array_name
          {
            socr,
            scr,
            su,
            sl,
            pcp,
          };

          template <typename src_vector_t>
          static void
          place_data (array_name name, vector_t &dst, const src_vector_t &src, scale_array_placement_info &placement_info)
          {
            if (src.empty ())
              return ;

            if (dst.size () <= 1)
              {
                copy (dst, src);

                if (name == socr)
                  placement_info.socr_step = 1, placement_info.socr_offset = 0;
                else if (name == scr)
                  placement_info.scr_step = 1, placement_info.scr_offset = 0;
                else if (name == su)
                  placement_info.su_step = 1, placement_info.su_offset = 0;
                else if (name == sl)
                  placement_info.sl_step = 1, placement_info.sl_offset = 0;
                else if (name == pcp)
                  placement_info.pcp_step = 1, placement_info.pcp_offset = 0;

                placement_info.size = (int)src.size ();
                return ;
              }

            BS_ASSERT ((int)src.size () == placement_info.size);
            if (name == socr)
              {
                if (placement_info.socr_step != 0)
                  {
                    BS_ASSERT (placement_info.socr_offset == 0);
                    set_stored (dst, src, placement_info.socr_offset, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);

                    if (placement_info.pcp_step != 0)
                      copy (dst, dst, placement_info.pcp_offset, placement_info.size);
                    if (placement_info.sl_step != 0)
                      copy (dst, dst, placement_info.sl_offset, placement_info.size);
                    if (placement_info.su_step != 0)
                      copy (dst, dst, placement_info.su_offset, placement_info.size);
                    if (placement_info.scr_step != 0)
                      copy (dst, dst, placement_info.scr_offset, placement_info.size);

                    placement_info.socr_offset = 0;

                    copy_src (dst, src, placement_info.socr_offset, placement_info.size);
                    placement_info.socr_step = 1;
                  }
              }
            else if (name == scr)
              {
                if (placement_info.scr_step != 0)
                  {
                    BS_ASSERT (placement_info.scr_offset <= placement_info.size);
                    set_stored (dst, src, placement_info.scr_offset, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);

                    if (placement_info.pcp_step != 0)
                      copy (dst, dst, placement_info.pcp_offset, placement_info.size);
                    if (placement_info.sl_step != 0)
                      copy (dst, dst, placement_info.sl_offset, placement_info.size);
                    if (placement_info.su_step != 0)
                      copy (dst, dst, placement_info.su_offset, placement_info.size);

                    placement_info.scr_offset = 0;
                    if (placement_info.socr_step != 0)
                      placement_info.scr_offset = 0;

                    copy_src (dst, src, placement_info.scr_offset, placement_info.size);
                    placement_info.scr_step = 1;
                  }
              }
            else if (name == su)
              {
                if (placement_info.su_step != 0)
                  {
                    BS_ASSERT (placement_info.su_offset <= 2 * placement_info.size);
                    set_stored (dst, src, placement_info.su_offset, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);

                    if (placement_info.pcp_step != 0)
                      copy (dst, dst, placement_info.pcp_offset, placement_info.size);
                    if (placement_info.sl_step != 0)
                      copy (dst, dst, placement_info.sl_offset, placement_info.size);

                    placement_info.su_offset = 0;
                    if (placement_info.socr_step != 0)
                      placement_info.su_offset += placement_info.size;
                    if (placement_info.scr_step != 0)
                      placement_info.su_offset += placement_info.size;

                    copy_src (dst, src, placement_info.su_offset, placement_info.size);
                    placement_info.su_step = 1;
                  }
              }
            else if (name == sl)
              {
                if (placement_info.sl_step != 0)
                  {
                    BS_ASSERT (placement_info.sl_offset <= 3 * placement_info.size);
                    set_stored (dst, src, placement_info.sl_offset, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);

                    if (placement_info.pcp_step != 0)
                      copy (dst, dst, placement_info.pcp_offset, placement_info.size);

                    placement_info.sl_offset = 0;
                    if (placement_info.socr_step != 0)
                      placement_info.sl_offset += placement_info.size;
                    if (placement_info.scr_step != 0)
                      placement_info.sl_offset += placement_info.size;
                    if (placement_info.su_step != 0)
                      placement_info.sl_offset += placement_info.size;

                    copy_src (dst, src, placement_info.sl_offset, placement_info.size);
                    placement_info.sl_step = 1;
                  }
              }
            else if (name == pcp)
              {
                if (placement_info.pcp_step != 0)
                  {
                    BS_ASSERT (placement_info.pcp_offset <= 4 * placement_info.size);
                    set_stored (dst, src, placement_info.pcp_offset, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);

                    placement_info.pcp_offset = 0;
                    if (placement_info.socr_step != 0)
                      placement_info.pcp_offset += placement_info.size;
                    if (placement_info.scr_step != 0)
                      placement_info.pcp_offset += placement_info.size;
                    if (placement_info.su_step != 0)
                      placement_info.pcp_offset += placement_info.size;
                    if (placement_info.sl_step != 0)
                      placement_info.pcp_offset += placement_info.size;

                    copy_src (dst, src, placement_info.pcp_offset, placement_info.size);
                    placement_info.pcp_step = 1;
                  }
              }
          }

          static void
          remove_data (array_name name, vector_t &dst, scale_array_placement_info &placement_info)
          {
            if (name == pcp)
              {
                if (placement_info.pcp_step != 0)
                  {
                    vector_t new_dst;

                    new_dst.resize (dst.size () - placement_info.size);

                    void *p_dst = &new_dst[0];
                    void *p_src = &dst[0];

                    memcpy (p_dst, p_src, placement_info.pcp_offset * sizeof (vector_t::value_type));
                    if (placement_info.pcp_offset < 4 * placement_info.size)
                      {
                        p_dst = &new_dst[placement_info.pcp_offset];
                        p_src = &dst[placement_info.pcp_offset + placement_info.size];
                        size_t size = (dst.size () - placement_info.size - placement_info.pcp_offset);

                        BS_ASSERT ((int)size >= 0) (dst.size ()) (placement_info.size) (placement_info.pcp_offset);
                        memcpy (p_dst, p_src, size * sizeof (vector_t::value_type));
                      }

                    dst.swap (new_dst);
                  }
              }
            else
              {
                BS_ASSERT (false && "NOT IMPL YET");
              }
          }

private:

          //static void
          //copy (vector_t &dst, const vector_t &src)
          //{
          //  dst.clear ();
          //  dst.assign (src.begin (), src.end ());
          //}
          template <typename src_t>
          static void
          copy (vector_t &dst, const src_t &src)
          {
            dst.clear ();
            for (size_t i = 0, cnt = src.size (); i < cnt; ++i)
              dst.push_back (src[i]);
          }

          //static void
          //copy (vector_t &dst, const vector_t &src, int &offset, int size)
          //{
          //  memcpy (&dst[offset + size], &src[offset], sizeof (vector_t::value_type) * size);
          //  offset += size;
          //}

          //static void
          //copy_src (vector_t &dst, const vector_t &src, int offset, int size)
          //{
          //  memcpy (&dst[offset], &src[0], sizeof (vector_t::value_type) * size);
          //}

          //static void
          //set_stored (vector_t &dst, const vector_t &src, int offset, int size)
          //{
          //  BS_ASSERT ((int)src.size () == size);

          //  memcpy (&dst[offset], &src[0], sizeof (vector_t::value_type) * size);
          //}

          template <typename src_t>
          static void
          copy (vector_t &dst, const src_t &src, int &offset, int size)
          {
            for (size_t i = 0, cnt = size; i < cnt; ++i)
              {
                dst[offset + size + i] = src[offset + i];
              }

            offset += size;
          }

          template <typename src_t>
          static void
          copy_src (vector_t &dst, const src_t &src, int offset, int size)
          {
            for (size_t i = 0, cnt = size; i < cnt; ++i)
              {
                dst[offset + i] = src[i];
              }
          }

          template <typename src_t>
          static void
          set_stored(vector_t &dst, const src_t &src, int offset, int size)
          {
            BS_ASSERT ((int)src.size () == size);

            for (size_t i = 0, cnt = size; i < cnt; ++i)
              {
                dst[offset + i] = src[i];
              }
          }

        };


      struct struct_like_scale_arrays_t
        {
          typedef shared_vector <float> vector_t;

          enum field_name
          {
            socr,
            scr,
            su,
            sl,
            pcp,

            field_count,
          };

          template <typename src_vector_t>
          static void
          place_data (field_name name, vector_t &dst, const src_vector_t &src, scale_array_placement_info &placement_info)
          {
            if (src.empty ())
              return ;

            if (dst.size () <= 1)
              {
                dst.clear ();
                dst.insert (dst.end (), src.begin (), src.end ());

                if (name == socr)
                  placement_info.socr_step = 1, placement_info.socr_offset = 0;
                else if (name == scr)
                  placement_info.scr_step = 1, placement_info.scr_offset = 0;
                else if (name == su)
                  placement_info.su_step = 1, placement_info.su_offset = 0;
                else if (name == sl)
                  placement_info.sl_step = 1, placement_info.sl_offset = 0;
                else if (name == pcp)
                  placement_info.pcp_step = 1, placement_info.pcp_offset = 0;

                placement_info.size = (int)src.size ();
                return ;
              }

            BS_ASSERT ((int)src.size () == placement_info.size);
            if (name == socr)
              {
                if (placement_info.socr_step)
                  {
                    set_stored (dst, src, placement_info.socr_offset, placement_info.socr_step, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);
                    int i = 0;
                    i += placement_info.scr_step != 0;
                    i += placement_info.su_step != 0;
                    i += placement_info.sl_step != 0;
                    i += placement_info.pcp_step != 0;

                    copy (dst, src, 0, i, placement_info.size);

                    placement_info.socr_step = i;
                    placement_info.socr_offset = 0;
                  }
              }
            else if (name == scr)
              {
                if (placement_info.scr_step)
                  {
                    set_stored (dst, src, placement_info.scr_offset, placement_info.scr_step, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);
                    int i = 0, j = 0;
                    j += placement_info.socr_step != 0;
                    i += placement_info.su_step != 0;
                    i += placement_info.sl_step != 0;
                    i += placement_info.pcp_step != 0;

                    copy (dst, src, j, i, placement_info.size);

                    placement_info.scr_step = j + i;
                    placement_info.scr_offset = 0;
                  }
              }
            else if (name == su)
              {
                if (placement_info.su_step)
                  {
                    set_stored (dst, src, placement_info.su_offset, placement_info.su_step, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);
                    int i = 0, j = 0;
                    j += placement_info.socr_step != 0;
                    j += placement_info.scr_step != 0;
                    i += placement_info.sl_step != 0;
                    i += placement_info.pcp_step != 0;
                    copy (dst, src, j, i, placement_info.size);

                    placement_info.su_step = j + i;
                    placement_info.su_offset = 0;
                  }
              }
            else if (name == sl)
              {
                if (placement_info.sl_step)
                  {
                    set_stored (dst, src, placement_info.sl_offset, placement_info.sl_step, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);
                    int i = 0, j = 0;
                    j += placement_info.socr_step != 0;
                    j += placement_info.scr_step != 0;
                    j += placement_info.su_step != 0;
                    i += placement_info.pcp_step != 0;
                    copy (dst, src, j, i, placement_info.size);

                    placement_info.sl_step = j + i;
                    placement_info.sl_offset = 0;
                  }
              }
            else if (name == pcp)
              {
                if (placement_info.pcp_step)
                  {
                    set_stored (dst, src, placement_info.sl_offset, placement_info.sl_step, placement_info.size);
                  }
                else
                  {
                    dst.resize (dst.size () + placement_info.size);
                    int j = 0;
                    j += placement_info.socr_step != 0;
                    j += placement_info.scr_step != 0;
                    j += placement_info.su_step != 0;
                    j += placement_info.sl_step != 0;
                    copy (dst, src, j, 0, placement_info.size);

                    placement_info.pcp_step = j;
                    placement_info.pcp_offset = 0;
                  }
              }
          }

          static void
          copy (vector_t &dst, const vector_t &src, int before, int after, int size)
          {
            int count = before + after;
            for (int i = size - 1; i >= 0; --i)
              {
                for (int j = before + after - 1; j >= before + 1; --j)
                  dst [i * (count + 1) + j + 1] = dst [i * count + j];

                dst [i * (count + 1) + before] = src [i];

                for (int j = before - 1; j >= 0; --j)
                  dst [i * (count + 1) + j] = dst [i * count + j];
              }
          }

          static void
          set_stored (vector_t &dst, const vector_t &src, int offset, int step, int size)
          {
            vector_t::value_type *d = &dst[offset];
            for (int si = 0, di = 0; si < size; ++si, di += step)
              {
                d[di] = src[si];
              }
          }

          static void
          copy (vector_t &dst, const shared_vector <double> &src, int before, int after, int size)
          {
            int count = before + after;
            for (int i = size - 1; i >= 0; --i)
              {
                for (int j = before + after - 1; j >= before + 1; --j)
                  dst [i * (count + 1) + j + 1] = dst [i * count + j];

                dst [i * (count + 1) + before] = src [i];

                for (int j = before - 1; j >= 0; --j)
                  dst [i * (count + 1) + j] = dst [i * count + j];
              }
          }

          static void
          set_stored (vector_t &dst, const shared_vector <double> &src, int offset, int step, int size)
          {
            vector_t::value_type *d = &dst[offset];
            for (int si = 0, di = 0; si < size; ++si, di += step)
              {
                d[di] = src[si];
              }
          }

        };


    } // namespace data_placement
  } // namespace scal
} // namespace blue_sky


#endif  // #ifndef BS_SCALE_ARRAYS_PLACEMENT_STRATEGIES_H_
