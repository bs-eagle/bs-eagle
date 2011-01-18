#include "bs_bos_core_data_storage_stdafx.h"

#include "py_bs_array_map.h"
namespace blue_sky
  {
  namespace python
    {
    void py_export_array_maps ()
    {
      using namespace boost::python;
      base_exporter <bs_array_map <base_strategy_did::i_type_t, base_strategy_did::i_type_t>, bs_array_map_exporter>::export_class ("bs_array_map_i"); 
      base_exporter <bs_array_map <base_strategy_fif::i_type_t, base_strategy_fif::fp_type_t>, bs_array_map_exporter>::export_class ("bs_array_map_fp"); 
    }
  }
}
