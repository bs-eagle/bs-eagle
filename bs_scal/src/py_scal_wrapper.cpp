/**
 * \file py_scal_wrapper.cpp
 * \brief
 * \author
 * \date
 * */
#include "bs_scal_stdafx.h"

#include "py_scal_wrapper.h"
#include "scal_save_data.h"

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

    template <typename strategy_t>
    void py_export_scale_array_holder (const char *name)
    {
      using namespace boost::python;

      class_ < py_scale_array_holder<strategy_t> > (name)
      .def ("insert_socr", &py_scale_array_holder<strategy_t>::insert_socr)
      .def ("insert_scr", &py_scale_array_holder<strategy_t>::insert_scr)
      .def ("insert_sl", &py_scale_array_holder<strategy_t>::insert_sl)
      .def ("insert_su", &py_scale_array_holder<strategy_t>::insert_su)
      ;
    }

    template <typename strategy_t>
    void py_export_scal_data_holder (const char *name)
    {
      using namespace boost::python;

#ifdef _DEBUG
      class_ < py_scal_data_holder<strategy_t> > (name)
      .def ("save_raw_data", &py_scal_data_holder<strategy_t>::save_raw_data)
      .def ("save_data", &py_scal_data_holder<strategy_t>::save_data)
#endif
      ;
    }

    template <typename strategy_t>
    void py_export_scal_3p (const char *name)
    {
      using namespace boost::python;

      class_ < py_scal_3p<strategy_t> > (name)
      .def ("get_relative_perm",  &py_scal_3p<strategy_t>::get_relative_perm)
      .def ("get_capillary",      &py_scal_3p<strategy_t>::get_capillary)
      .def ("get_water_scale",    &py_scal_3p<strategy_t>::get_water_scale)
      .def ("get_gas_scale",      &py_scal_3p<strategy_t>::get_gas_scale)
      .def ("get_water_data",     &py_scal_3p<strategy_t>::get_water_data)
      .def ("get_gas_data",       &py_scal_3p<strategy_t>::get_gas_data)
      //.def ("set_water_jfunction", &py_scal_3p::set_water_jfunction)
      //.def ("set_gas_jfunction", &py_scal_3p::set_gas_jfunction)
      //.def ("set_n_phases", &py_scal_3p::set_n_phases)
      //.def ("set_phases", &py_scal_3p::set_phases)
      //.def ("set_rpo_model", &py_scal_3p::set_rpo_model)
      //.def ("process_model", &py_scal_3p<strategy_t>::process_model)
      //.def ("set_test_model", &py_scal_3p<strategy_t>::set_test_model)
      //.def ("set_phase_d", &py_scal_3p::set_phase_d)
      //.def ("set_sat_d", &py_scal_3p::set_sat_d)
      ;

    }

    //template <typename strategy_t>
    //void py_export_test_model (const char *name)
    //{
    //  using namespace boost::python;

    //  class_ < py_test_model<strategy_t> > (name)
    //  .def ("init", &py_test_model<strategy_t>::init)
    //  ;
    //}

    template <typename strategy_t>
    void py_export_jfunction (const char *name)
    {
      using namespace boost::python;

      class_ <py_jfunction <strategy_t>, bases <py_objbase> > (name)
      .add_property ("st_phase", &py_jfunction<strategy_t>::get_st_phase, &py_jfunction<strategy_t>::set_st_phase)
      .add_property ("alpha", &py_jfunction<strategy_t>::get_alpha, &py_jfunction<strategy_t>::set_alpha)
      .add_property ("beta", &py_jfunction<strategy_t>::get_beta, &py_jfunction<strategy_t>::set_beta)
      .add_property ("valid", &py_jfunction<strategy_t>::valid)
      .def ("init", &py_jfunction<strategy_t>::init)
      ;
    }

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

      py_export_jfunction <base_strategy_fi> ("jfunction_fi");
      py_export_jfunction <base_strategy_di> ("jfunction_di");

      py_export_scale_array_holder <base_strategy_fi> ("scale_arrays_fi");
      py_export_scale_array_holder <base_strategy_di> ("scale_arrays_di");

      py_export_scal_data_holder <base_strategy_fi> ("scal_data_fi");
      py_export_scal_data_holder <base_strategy_di> ("scal_data_di");

      py_export_scal_3p <base_strategy_fi> ("scal_3p_fi");
      py_export_scal_3p <base_strategy_di> ("scal_3p_di");
    }

  }
}

