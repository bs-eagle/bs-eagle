/**
 * \file py_scal_wrapper.cpp
 * \brief
 * \author
 * \date
 * */
#include "bs_scal_stdafx.h"

#include "py_scal_wrapper.h"
#include "scal_save_data.h"

#include "export_python_wrapper.h"

namespace blue_sky
  {
  namespace python
    {

    template <typename strategy_t>
    py_scal_data_holder<strategy_t>::py_scal_data_holder ()
        : base_t (wrapped_t::bs_type ())
    {
    }

    template <typename strategy_t>
    py_scal_data_holder<strategy_t>::py_scal_data_holder (sp_scal_data_holder_t holder)
        : base_t (holder)
    {
    }

#ifdef _DEBUG
    template <typename strategy_t>
    void
    py_scal_data_holder<strategy_t>::save_raw_data (const std::string &filename)
    {
      sp_scal_data_holder_t spx = get_spx <wrapped_t> ();
      BS_ASSERT (spx) (filename);

      const scal::data_placement::scal_placement_info &info = spx->get_placement_info ();
      const typename wrapped_t::item_array_t &raw_data              = spx->get_raw_data ();
      const typename wrapped_t::region_vector_t &region_info         = spx->get_region_info ();

      using namespace scal;
      using namespace data_placement;
      switch (info.type)
        {
        case spof:
          save_spof_raw_data <strategy_t> (raw_data, region_info, filename);
          break;
        case spfn_sof3:
          save_spfn_sof3_raw_data <strategy_t> (raw_data, region_info, filename);
          break;
        case sof3_spfn:
          save_sof3_spfn_raw_data <strategy_t> (raw_data, region_info, filename);
          break;
        default:
          bs_throw_exception ("Unknown type");
        }
    }

    template <typename strategy_t>
    void
    py_scal_data_holder<strategy_t>::save_data (const std::string &filename)
    {
      sp_scal_data_holder_t spx (get_spx <wrapped_t> ());
      BS_ASSERT (spx) (filename);

      const scal::data_placement::scal_placement_info &info = spx->get_placement_info ();
      const typename wrapped_t::region_vector_t &region_info         = spx->get_region_info ();

      using namespace scal;
      using namespace data_placement;
      switch (info.type)
        {
        case spof:
          save_spof_data <strategy_t> (*spx, region_info, filename);
          break;
        case spfn_sof3:
          save_spfn_sof3_data <strategy_t> (*spx, region_info, filename);
          break;
        case sof3_spfn:
          save_sof3_spfn_data <strategy_t> (*spx, region_info, filename);
          break;
        default:
          bs_throw_exception ("Unknown type");
        }
    }
#endif

    PY_EXPORTER (scale_array_holder_exporter, default_exporter)
    PY_EXPORTER_END;

    PY_EXPORTER (scal_data_holder_exporter, default_exporter)
    PY_EXPORTER_END;

    PY_EXPORTER (jfunction_exporter, default_exporter)
      .add_property ("st_phase",  &T::st_phase)
      .add_property ("alpha",     &T::alpha)
      .add_property ("beta",      &T::beta)
      .add_property ("valid",     &T::valid)
    PY_EXPORTER_END;

    PY_EXPORTER (scal_3p_exporter, default_exporter)
      .add_property ("water_scale",      &T::get_water_scale)
      .add_property ("gas_scale",        &T::get_gas_scale)
      .add_property ("water_data",       &T::get_water_data)
      .add_property ("gas_data",         &T::get_gas_data)
      .add_property ("water_jfunction",  &T::get_water_jfunction, &T::set_water_jfunction)
      .add_property ("gas_jfunction",    &T::get_gas_jfunction, &T::set_gas_jfunction)
    PY_EXPORTER_END;

    void py_export_scal ()
    {
      using namespace boost::python;

      enum_ <PHASE_ENUM> ("phase_enum")
        .value ("water", PHASE_WATER)
        .value ("gas", PHASE_GAS)
        .value ("oil", PHASE_OIL)
        ;

      enum_ <RPO_MODEL_ENUM> ("rpo_model_enum")
        .value ("default_model", RPO_DEFAULT_MODEL)
        .value ("stone1_model", STONE1_MODEL)
        .value ("stone2_model", STONE2_MODEL)
        ;

      enum_ <JFUNC_PERM_TYPE_ENUM> ("jfunc_perm_type_enum")
        .value ("xy", JFUNC_PERM_XY)
        .value ("x", JFUNC_PERM_X)
        .value ("y", JFUNC_PERM_Y)
        .value ("z", JFUNC_PERM_Z)
        ;

      strategy_exporter::export_base <jfunction, jfunction_exporter> ("jfunction");
      strategy_exporter::export_base <scale_array_holder, scale_array_holder_exporter> ("scale_arrays");
      strategy_exporter::export_base <scal_2p_data_holder, scal_data_holder_exporter> ("scal_data");
      strategy_exporter::export_base <scal_3p, scal_3p_exporter> ("scal_3p");
    }
  }
}

