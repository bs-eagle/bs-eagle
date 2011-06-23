/**
 *
 * */
#include "stdafx.h"
#include "calc_model.h"

namespace blue_sky
  {

  void
  calc_model_data_tmp_holder::save (const BS_SP (calc_model) &cm)
  {
    pressure.assign (cm->pressure.begin (), cm->pressure.end ());
    saturation_3p.assign (cm->saturation_3p.begin (), cm->saturation_3p.end ());

    if (cm->is_gas ())
      {
        gas_oil_ratio.assign (cm->gas_oil_ratio.begin (), cm->gas_oil_ratio.end ());
        main_var.assign (cm->main_variable.begin (), cm->main_variable.end ());
      }
  }

  void
  calc_model_data_tmp_holder::restore (BS_SP (calc_model) &cm)
  {
    cm->pressure.assign (pressure.begin (), pressure.end ());
    cm->saturation_3p.assign (saturation_3p.begin (), saturation_3p.end ());

    if (cm->is_gas ())
      {
        cm->gas_oil_ratio.assign (gas_oil_ratio.begin (), gas_oil_ratio.end ());
        cm->main_variable.assign (main_var.begin (), main_var.end ());
      }
  }

} // namespace blue_sky

