/**
 * \file scal_save_data.h
 * \brief save scal data into file
 * \author Sergey Miryanov
 * \date 28.08.2008
 * */
#ifndef BS_SCAL_SAVE_DATA_H_
#define BS_SCAL_SAVE_DATA_H_

#ifdef _DEBUG

namespace blue_sky
  {
  namespace scal
    {
    namespace data_placement
      {

      struct scal_data_saver
        {
          typedef scal_data_saver this_t;

          scal_data_saver (const std::string &filename)
              : file_ (0)
              , lkeeper_ ("C", LC_ALL)
          {
            file_ = fopen (filename.c_str (), "w");
            BS_ASSERT (file_) (filename);
          }

          ~scal_data_saver ()
          {
            if (file_)
              fclose (file_);
          }

          this_t &
          save_section_name (const std::string &name)
          {
            if (name.length ())
              fprintf (file_, "%s\n", name.c_str ());

            return *this;
          }

          this_t &
          save_section_end ()
          {
            fprintf (file_, "/\n");
            return *this;
          }

          this_t &
          save_item (const float &item)
          {
            fprintf (file_, "\t%10.20f", item);
            return *this;
          }

          this_t &
          save_item (const double &item)
          {
            fprintf (file_, "\t%10.20lf", item);
            return *this;
          }

          this_t &
          save_endl ()
          {
            fprintf (file_, "\n");
            return *this;
          }

          FILE            *file_;
          locale_keeper   lkeeper_;
        };

      template <typename strategy_t>
      void
      save_spof_raw_data (const typename scal_2p_data_holder<strategy_t>::item_array_t &raw_data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);

        file.save_section_name ("SPOF");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.Sp_count == region.So_count) (region.Sp_count) (region.So_count);
            BS_ASSERT (region.sp_offset == region.so_offset) (region.sp_offset) (region.so_offset);

            for (size_t j = 0, jcnt = region.Sp_count * 5; j < jcnt; j+=5)
              {
                file
                .save_item (raw_data[region.sp_offset + j + 0])
                .save_item (raw_data[region.sp_offset + j + 2])
                .save_item (raw_data[region.sp_offset + j + 3])
                .save_item (raw_data[region.sp_offset + j + 4])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_spof_data (scal_2p_data_holder<strategy_t> &data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);
        file.save_section_name ("SPOF");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.Sp_count == region.So_count) (region.Sp_count) (region.So_count);
            BS_ASSERT (region.sp_offset == region.so_offset) (region.sp_offset) (region.so_offset);

            const scal_region<strategy_t> &region_data = data.get_region ((typename strategy_t::index_t)i);

            BS_ASSERT (region_data.Sp.size () == region_data.So.size ()) (region_data.Sp.size ()) (region_data.So.size ());
            BS_ASSERT (region_data.Krp.size () == region_data.Krop.size ()) (region_data.Krp.size ()) (region_data.Krop.size ());
            BS_ASSERT (region_data.Sp.size () == region_data.Pcp.size ()) (region_data.Sp.size ()) (region_data.Pcp.size ());
            BS_ASSERT (region_data.Sp.size () == region_data.Krp.size ()) (region_data.Sp.size ()) (region_data.Krp.size ());

            for (size_t j = 0, jcnt = region_data.Sp.size (); j < jcnt; ++j)
              {
                file
                .save_item (region_data.Sp[j])
                .save_item (region_data.Krp[j])
                .save_item (region_data.Krop[j])
                .save_item (region_data.Pcp[j])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_spfn_raw_data (const typename scal_2p_data_holder<strategy_t>::item_array_t &raw_data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, scal_data_saver &file)
      {
        file.save_section_name ("SPFN");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.Sp_count != -1) (region.Sp_count);

            for (size_t j = 0, jcnt = region.Sp_count * 3; j < jcnt; j+=3)
              {
                file
                .save_item (raw_data[region.sp_offset + j + 0])
                .save_item (raw_data[region.sp_offset + j + 1])
                .save_item (raw_data[region.sp_offset + j + 2])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_sof3_raw_data (const typename scal_2p_data_holder<strategy_t>::item_array_t &raw_data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, scal_data_saver &file)
      {
        file.save_section_name ("SOF3");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.So_count != -1) (region.So_count);

            for (size_t j = 0, jcnt = region.So_count * 3; j < jcnt; j+=3)
              {
                file
                .save_item (raw_data[region.so_offset + j + 0])
                .save_item (raw_data[region.so_offset + j + 1])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_spfn_data (scal_2p_data_holder<strategy_t> &data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, scal_data_saver &file)
      {
        file.save_section_name ("SPFN");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.Sp_count != -1) (region.Sp_count);

            const scal_region<strategy_t> &region_data = data.get_region ((typename strategy_t::index_t)i);

            BS_ASSERT (region_data.Sp.size () == region_data.Krp.size ()) (region_data.Sp.size ()) (region_data.Krp.size ());
            BS_ASSERT (region_data.Sp.size () == region_data.Pcp.size ()) (region_data.Sp.size ()) (region_data.Pcp.size ());

            for (size_t j = 0, jcnt = region_data.Sp.size (); j < jcnt; ++j)
              {
                file
                .save_item (region_data.Sp[j])
                .save_item (region_data.Krp[j])
                .save_item (region_data.Pcp[j])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_sof3_data (scal_2p_data_holder<strategy_t> &data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, scal_data_saver &file)
      {
        file.save_section_name ("SOF3");

        for (size_t i = 0, cnt = region_info.size (); i < cnt; ++i)
          {
            const scal_region_info<strategy_t> &region = region_info[i];

            BS_ASSERT (region.So_count != -1) (region.So_count);

            const scal_region<strategy_t> &region_data = data.get_region ((typename strategy_t::index_t)i);

            BS_ASSERT (region_data.So.size () == region_data.Krop.size ()) (region_data.So.size ()) (region_data.Krop.size ());

            for (size_t j = 0, jcnt = region_data.So.size (); j < jcnt; ++j)
              {
                file
                .save_item (region_data.So[j])
                .save_item (region_data.Krop[j])
                .save_endl ();
              }

            file.save_section_end ();
          }
      }

      template <typename strategy_t>
      void
      save_spfn_sof3_raw_data (const typename scal_2p_data_holder<strategy_t>::item_array_t &raw_data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);

        save_spfn_raw_data <strategy_t> (raw_data, region_info, file);
        save_sof3_raw_data <strategy_t> (raw_data, region_info, file);
      }
      template <typename strategy_t>
      void
      save_spfn_sof3_data (scal_2p_data_holder<strategy_t> &data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);

        save_spfn_data <strategy_t> (data, region_info, file);
        save_sof3_data <strategy_t> (data, region_info, file);
      }

      template <typename strategy_t>
      void
      save_sof3_spfn_raw_data (const typename scal_2p_data_holder<strategy_t>::item_array_t &raw_data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);

        save_sof3_raw_data <strategy_t> (raw_data, region_info, file);
        save_spfn_raw_data <strategy_t> (raw_data, region_info, file);
      }
      template <typename strategy_t>
      void
      save_sof3_spfn_data (scal_2p_data_holder<strategy_t> &data, const typename scal_2p_data_holder<strategy_t>::region_vector_t &region_info, const std::string &filename)
      {
        scal_data_saver file (filename);

        save_sof3_data <strategy_t> (data, region_info, file);
        save_spfn_data <strategy_t> (data, region_info, file);
      }

    } // namespace data_placement
  } // namespace scal
} // namespace blue_sky

#endif  // #ifdef _DEBUG
#endif  // #ifndef BS_SCAL_SAVE_DATA_H_

