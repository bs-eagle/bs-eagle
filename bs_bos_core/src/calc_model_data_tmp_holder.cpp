/**
 *
 * */
#include "stdafx.h"
#include "calc_model.h"

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky
  {

  template <typename strategy_t>
  void
  calc_model_data_tmp_holder<strategy_t>::save (const sp_calc_model_t &cm)
  {
    pressure.assign (cm->pressure.begin (), cm->pressure.end ());
    saturation_3p.assign (cm->saturation_3p.begin (), cm->saturation_3p.end ());

    if (cm->is_gas ())
      {
        gas_oil_ratio.assign (cm->gas_oil_ratio.begin (), cm->gas_oil_ratio.end ());
        main_var.assign (cm->main_variable.begin (), cm->main_variable.end ());
      }
  }

  template <typename strategy_t>
  void
  calc_model_data_tmp_holder<strategy_t>::restore (sp_calc_model_t &cm)
  {
    cm->pressure.assign (pressure.begin (), pressure.end ());
    cm->saturation_3p.assign (saturation_3p.begin (), saturation_3p.end ());

    if (cm->is_gas ())
      {
        cm->gas_oil_ratio.assign (gas_oil_ratio.begin (), gas_oil_ratio.end ());
        cm->main_variable.assign (main_var.begin (), main_var.end ());
      }
  }

  template struct calc_model_data_tmp_holder<base_strategy_fi>;
  template struct calc_model_data_tmp_holder<base_strategy_di>;
  template struct calc_model_data_tmp_holder<base_strategy_mixi>;

} // namespace blue_sky

