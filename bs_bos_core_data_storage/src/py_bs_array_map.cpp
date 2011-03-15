#include "bs_bos_core_data_storage_stdafx.h"

#include "py_bs_array_map.h"
namespace blue_sky
  {
  namespace python
    {
    void py_export_array_maps ()
    {
      using namespace boost::python;
      //base_exporter <bs_array_map <t_int, t_int>, bs_array_map_exporter>::export_class ("bs_array_map_i"); 
      //base_exporter <bs_array_map <t_int, t_double>, bs_array_map_exporter>::export_class ("bs_array_map_fp"); 
    }
  }
}
