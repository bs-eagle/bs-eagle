#include "bs_bos_core_data_storage_stdafx.h"
#include "prvd_table.h"

namespace blue_sky
  {
  val_vs_depth::val_vs_depth() {}

  val_vs_depth::~val_vs_depth()
  {
    depth.clear();
    values.clear();
  }

  size_t val_vs_depth::get_table_len () const
    {
      BS_ASSERT(depth.size() == values.size());
      return depth.size();
    }

  void val_vs_depth::set_table_len (unsigned len)
  {
    if (len < 1)
      throw bs_exception("val_vs_depth","try to set length < 1");
    depth.clear();
    values.clear();
    for (;len > 0;--len)
      {
        depth.push_back(-1);
        values.push_back(-1);
      }
  }

  val_vs_depth::d_vec &val_vs_depth::tdepth()
  {
    return depth;
  }

  val_vs_depth::d_vec &val_vs_depth::tvalues()
  {
    return values;
  }

  /*!
  	\brief check monotonic of pressure and depth
  	\param region_number -- pvt region number for error printing
  	\return 0 if success, <0 if error occur
  */
  int val_vs_depth::check_monotonic (int region_number) const
    {
      BS_ASSERT(depth.size() == values.size());

      if (depth.size() < 2)
        return 0;

      for (unsigned i = 0; i < depth.size(); ++i)
        {
          if ((depth[i] - depth[i-1] < EPS) || (values[i] - values[i-1] < EPS))
            {
              BOSERR (section::check_data, level::error) 
                << "Input error: values "
                << (i - 1) << " and " << i
                << " for depth (keyword PRVD or RSVD) in region"
                <<  region_number << " are nonmonotonic or equal\n"
                << bs_end;
              return -1;
            }
        }
      return 0;
    }
}
