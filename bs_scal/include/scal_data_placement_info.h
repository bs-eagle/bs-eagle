/**
 * \file scal_data_placement_info.h
 * \brief information about scal data placement, scal::data_vectors are created based on this info
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#ifndef BS_SCAL_DATA_PLACEMENT_INFO_H_
#define BS_SCAL_DATA_PLACEMENT_INFO_H_

namespace blue_sky
  {
  namespace scal
    {
    namespace data_placement
      {

      enum scal_data_placement_type
      {
        scal_data_placement_null,
        separate_vectors,
        spkrp_sokrop_pcp,
        spof,
        spfn_sof3,
        sof3_spfn,
      };

      struct scal_placement_info
        {
          typedef unsigned char byte;

          byte  sp_step;
          byte  so_step;
          byte  krp_step;
          byte  krop_step;
          byte  pcp_step;
          int   sp_offset;                                ///< offset of Sp from begin of memory region
          int   so_offset;                                ///< offset of So from begin of memory region
          int   krp_offset;                               ///< offset of Krp from begin of memory region
          int   krop_offset;                              ///< offset of Krop from begin of memory region
          int   pcp_offset;                               ///< offset of Pcp from begin of memory region
          scal_data_placement_type  type;

          scal_placement_info ()
          {
            sp_step     = 0;
            so_step     = 0;
            krp_step    = 0;
            krop_step   = 0;
            pcp_step    = 0;

            sp_offset   = 0;
            so_offset   = 0;
            krp_offset  = 0;
            krop_offset = 0;
            pcp_offset  = 0;

            type        = scal_data_placement_null;
          }
        };

      struct scale_array_placement_info
        {
          typedef unsigned char byte;

          byte socr_step;
          byte scr_step;
          byte su_step;
          byte sl_step;
          byte pcp_step;

          int socr_offset;
          int scr_offset;
          int su_offset;
          int sl_offset;
          int pcp_offset;

          int size;

          scale_array_placement_info ()
          {
            socr_step   = 0;
            scr_step    = 0;
            su_step     = 0;
            sl_step     = 0;
            pcp_step    = 0;

            socr_offset = 0;
            scr_offset  = 0;
            su_offset   = 0;
            sl_offset   = 0;
            pcp_offset  = 0;

            size        = 0;
          }
        };


    } // namespace data_placement
  } // namespace scal
} // namespace blue_sky


#endif  // #ifndef BS_SCAL_DATA_PLACEMENT_INFO_H_
