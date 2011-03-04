/** 
 * \file scal_region_info.h
 * \brief scal_region_info definition
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_SCAL_SCAL_REGION_INFO_H_
#define BS_SCAL_SCAL_REGION_INFO_H_

namespace blue_sky {

  /**
   * \brief describe a scal region
   * */
  template <typename item_t>
  struct scal_region_info
  {
    //typedef typename strategy_t::item_t item_t;
    int So_count;                                   ///< count of So and Krop components of region
    int Sp_count;                                   ///< count of Sp and Krp components of region
    int Krp_min_greater_zero;                       ///< index of item that first greater than zero in Krp
    int Krop_min_greater_zero;                      ///< index of item that first greater than zero in Krop
    size_t so_offset;                               ///< offset for so, krop components from begin of allocated memory region
    size_t sp_offset;                               ///< offset for sp, krp and pcp components from begin of allocated memory region

    item_t spr;                                     ///< phase critical saturation (this is the highest saturation value for which the relative permeability is zero)
    item_t sorp;                                    ///< oil critical saturation (this is the highest oil saturation value for which the oil relative permeability is zero)
    item_t kpr;
    item_t krorp;
    item_t pcp_max;

    size_t size () const
    {
      return 2 * So_count + 3 * Sp_count;
    }
  };

}


#endif  // #ifndef BS_SCAL_SCAL_REGION_INFO_H_
