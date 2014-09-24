#include "bs_common.h"
#include BS_FORCE_PLUGIN_IMPORT()
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT()

#include "prvd_table.h"

#ifdef EPS
#undef EPS
#endif // EPS
#define EPS 1.e-10

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


  /*inline*/ double 
  val_vs_depth::interpolate_linear (double val) const
    {
      BS_ASSERT(depth.size() == values.size());
      if (!depth.size())
        return 0;

      //if value less than smallest depth - extrapolate
      //it as constant
      if (depth[0] >= val)
        return values[0];

      //if value more than largest depth - extrapolate
      //it as constant
      if (depth[depth.size() - 1] <= val)
        return values[values.size() - 1];

      //binary search algorithm
      size_t first = 0, last = depth.size() - 1, mid = 0;
      while (first < last)
        {
          mid = (last + first) / 2;
          if ((last - first) == 1)
            {
              if (fabs (depth[last] - depth[first]) < EPS)
                return values[first];
              return values[first] + (val - depth[first]) /
                     (depth[last] - depth[first]) *
                     (values[last] - values[first]);
            }
          else if (depth[mid] <= val)
            {
              first = mid + 1;
              continue;
            }
          else
            {
              last = mid - 1;
              continue;
            }
        }
      return -1;
    }

  BS_SP( table_iface)
  val_vs_depth::get_table_data () const
  {
    BS_SP( table_iface) table;
    table = BS_KERNEL.create_object ("table");
    table->init (depth.size(), 2);
    for (unsigned i = 0; i < depth.size(); ++i)
      {
        table->set_value(i, 0, depth[i]);
        table->set_value(i, 1, values[i]);
      }
    return table;
  }
}
