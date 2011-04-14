/**
 *       \file  calc_model_data_accessors.h
 *      \brief  Accessors for calc_model data, can check invalid access to 
 *              data at compile time.
 *    \details  Uses Boost::Preprocessor to generate for each accessor 
 *              a number of template classes that specialized for valid
 *              phase values.
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Deprecated due performance reasons
 * */
#ifndef BS_CALC_MODEL_DATA_ACCESSORS_H_
#define BS_CALC_MODEL_DATA_ACCESSORS_H_

namespace blue_sky
  {

/////////////////////////////////////////////////////////////////////////////
#define DEF_FUN_I(r, name, phase)                                                                                   \
  template <>                                                                                                       \
  struct BOOST_PP_CAT (hx_, name) <phase>                                                                           \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d)                 \
    {                                                                                                               \
      BS_ASSERT (phase_d[phase] != -1) (phase);                                                                     \
      return data.name [phase_d[phase]];                                                                            \
    }                                                                                                               \
  };

#define DEF_FUN_2_III(r, name, phase1, phase2)                                                                      \
  template <>                                                                                                       \
  struct BOOST_PP_CAT (hx_, name) <phase1, phase2>                                                                  \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d, int n_phases)   \
    {                                                                                                               \
      BS_ASSERT (phase_d[phase1] != -1) (phase1);                                                                   \
      BS_ASSERT (phase_d[phase2] != -1) (phase2);                                                                   \
      return data.name [phase_d[phase1] * n_phases + phase_d[phase2]];                                              \
    }                                                                                                               \
  };

#define DEF_FUN_2_II(r, iter_data, phase1)                                                                          \
  DEF_FUN_2_III (r, BOOST_PP_TUPLE_ELEM (2, 0, iter_data), phase1, BOOST_PP_TUPLE_ELEM (2, 1, iter_data))

#define DEF_FUN_2_I(name, ts1, tp1, ts2, tp2)  BOOST_PP_CAT(DEF_FUN_2_I_, ts2) (name, ts1, tp1, tp2)

#define DEF_FUN_2_I_1(name, ts1, tp1, tp2)                                                                          \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (1, 0, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))

#define DEF_FUN_2_I_2(name, ts1, tp1, tp2)                                                                          \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (2, 0, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))   \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (2, 1, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))

#define DEF_FUN_2_I_3(name, ts1, tp1, tp2)                                                                          \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (3, 0, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))   \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (3, 1, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))   \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_2_II, (name, BOOST_PP_TUPLE_ELEM (3, 2, tp2)), BOOST_PP_TUPLE_TO_SEQ (ts1, tp1))

#define DEF_FUN_PAIR_I(r, iter_data, phase)                                                                         \
  template <>                                                                                                       \
  struct BOOST_PP_CAT (hx_, BOOST_PP_TUPLE_ELEM (3, 0, iter_data)) <phase>                                          \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d)                  \
    {                                                                                                               \
      return data.BOOST_PP_TUPLE_ELEM (3, 1, iter_data) [phase] * data.BOOST_PP_TUPLE_ELEM (3, 2, iter_data) [phase]; \
    }                                                                                                               \
  };


#define DEF_FUN(name, ts, tp)                                                                                       \
  template <FI_PHASE_ENUM phase_value>                                                                              \
  struct BOOST_PP_CAT (hx_, name)                                                                                   \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d)                 \
    {                                                                                                               \
      class inclomplete_type_1;                                                                                     \
      inclomplete_type_1 invalid_phase_value;                                                                       \
      return t_double ();                                                                                           \
    }                                                                                                               \
  };                                                                                                                \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_I, name, BOOST_PP_TUPLE_TO_SEQ (ts, tp))

#define DEF_FUN_2(name, ts1, tp1, ts2, tp2)                                                                         \
  template <FI_PHASE_ENUM phase1, FI_PHASE_ENUM phase2>                                                             \
  struct BOOST_PP_CAT (hx_, name)                                                                                   \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d, int n_phases)   \
    {                                                                                                               \
      class inclomplete_type_1;                                                                                     \
      inclomplete_type_1 invalid_phase_value;                                                                       \
      return t_double ();                                                                                           \
    }                                                                                                               \
  };                                                                                                                \
  DEF_FUN_2_I (name, ts1, tp1, ts2, tp2)

#define DEF_FUN_PAIR(name, name1, name2, ts, tp)                                                                    \
  template <FI_PHASE_ENUM phase>                                                                                    \
  struct BOOST_PP_CAT (hx_, name)                                                                                   \
  {                                                                                                                 \
    static t_double                                                                              \
    get (const calc_model_data &data, const boost::array <t_long, FI_PHASE_TOT> &phase_d)                 \
    {                                                                                                               \
      class inclomplete_type_1;                                                                                     \
      inclomplete_type_1 invalid_phase_value;                                                                       \
      return t_double ();                                                                                           \
    }                                                                                                               \
  };                                                                                                                \
  BOOST_PP_SEQ_FOR_EACH (DEF_FUN_PAIR_I, (name, name1, name2), BOOST_PP_TUPLE_TO_SEQ (ts, tp))

  DEF_FUN   (cap_pressure, 2, (FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (s_deriv_cap_pressure, 2, (FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (relative_perm, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN_2 (s_deriv_relative_perm, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS), 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (invers_fvf, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (p_deriv_invers_fvf, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (invers_viscosity, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (p_deriv_invers_viscosity, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (invers_visc_fvf, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (p_deriv_invers_visc_fvf, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (density, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (p_deriv_density, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (mobility, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN   (p_deriv_mobility,3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN_2 (s_deriv_mobility, 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS), 3, (FI_PHASE_OIL, FI_PHASE_WATER, FI_PHASE_GAS))

  DEF_FUN_PAIR (s_deriv_invers_fvf, p_deriv_invers_fvf, s_deriv_cap_pressure, 2, (FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN_PAIR (s_deriv_invers_viscosity, p_deriv_invers_viscosity, s_deriv_cap_pressure, 2, (FI_PHASE_WATER, FI_PHASE_GAS))
  DEF_FUN_PAIR (s_deriv_invers_visc_fvf, p_deriv_invers_visc_fvf, s_deriv_cap_pressure, 2, (FI_PHASE_WATER, FI_PHASE_GAS))

/////////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG

#define CAP_PRESSURE(data, phase_d, phase)                        hx_cap_pressure<phase>::get (data, phase_d)
#define S_DERIV_CAP_PRESSURE(data, phase_d, phase)                hx_s_deriv_cap_pressure<phase>::get (data, phase_d)
#define RELATIVE_PERM(data, phase_d, phase)                       hx_relative_perm<phase>::get (data, phase_d)
#define S_DERIV_RELATIVE_PERM(data, phase_d, n_phases, phase1, phase2)  hx_s_deriv_relative_perm<phase1, phase2>::get (data, phase_d, n_phases)
#define GAS_OIL_RATIO(data)                                       params.gas_oil_ratio[params.n_block]
#define P_DERIV_GAS_OIL_RATIO(data)                               data.p_deriv_gas_oil_ratio
#define INVERS_FVF(data, phase_d, phase)                          hx_invers_fvf<phase>::get (data, phase_d)
#define P_DERIV_INVERS_FVF(data, phase_d, phase)                  hx_p_deriv_invers_fvf<phase>::get (data, phase_d)
#define GOR_DERIV_INVERS_FVF(data)                                data.gor_deriv_invers_fvf
#define INVERS_VISCOSITY(data, phase_d, phase)                    hx_invers_viscosity<phase>::get (data, phase_d)
#define P_DERIV_INVERS_VISCOSITY(data, phase_d, phase)            hx_p_deriv_invers_viscosity<phase>::get (data, phase_d)
#define GOR_DERIV_INVERS_VISCOSITY(data)                          data.gor_deriv_invers_viscosity
#define INVERS_VISC_FVF(data, phase_d, phase)                     hx_invers_visc_fvf<phase>::get (data, phase_d)
#define P_DERIV_INVERS_VISC_FVF(data, phase_d, phase)             hx_p_deriv_invers_visc_fvf<phase>::get (data, phase_d)
#define GOR_DERIV_INVERS_VISC_FVF(data)                           data.gor_deriv_invers_visc_fvf
#define DENSITY(data, phase_d, phase)                             hx_density<phase>::get (data, phase_d)
#define P_DERIV_DENSITY(data, phase_d, phase)                     hx_p_deriv_density<phase>::get (data, phase_d)
#define GOR_DERIV_DENSITY(data)                                   data.gor_deriv_density
#define MOBILITY(data, phase_d, phase)                            hx_mobility<phase>::get (data, phase_d)
#define P_DERIV_MOBILITY(data, phase_d, phase)                    hx_p_deriv_mobility<phase>::get (data, phase_d)
#define S_DERIV_MOBILITY(data, phase_d, n_phases, phase1, phase2) hx_s_deriv_mobility<phase1, phase2>::get (data, phase_d, n_phases)
#define S_DERIV_INVERS_FVF(data, phase_d, phase)                  hx_s_deriv_invers_fvf<phase>::get (data, phase_d)
#define S_DERIV_INVERS_VISCOSITY(data, phase_d, phase)            hx_s_deriv_invers_viscosity<phase>::get (data, phase_d)
#define S_DERIV_INVERS_VISC_FVF(data, phase_d, phase)             hx_s_deriv_invers_visc_fvf<phase>::get (data, phase_d)

#else

#define CAP_PRESSURE(data, phase_d, phase)                        data.cap_pressure[phase_d[phase]]
#define S_DERIV_CAP_PRESSURE(data, phase_d, phase)                data.s_deriv_cap_pressure[phase_d[phase]]
#define RELATIVE_PERM(data, phase_d, phase)                       data.relative_perm[phase_d[phase]]
#define S_DERIV_RELATIVE_PERM(data, phase_d, n_phases, phase1, phase2)      data.s_deriv_relative_perm[phase_d[phase1] * n_phases + phase_d[phase2]]
#define GAS_OIL_RATIO(data)                                       params.gas_oil_ratio[params.n_block]
#define P_DERIV_GAS_OIL_RATIO(data)                               data.p_deriv_gas_oil_ratio
#define INVERS_FVF(data, phase_d, phase)                          data.invers_fvf[phase_d[phase]]
#define P_DERIV_INVERS_FVF(data, phase_d, phase)                  data.p_deriv_invers_fvf[phase_d[phase]]
#define GOR_DERIV_INVERS_FVF(data)                                data.gor_deriv_invers_fvf
#define INVERS_VISCOSITY(data, phase_d, phase)                    data.invers_viscosity[phase_d[phase]]
#define P_DERIV_INVERS_VISCOSITY(data, phase_d, phase)            data.p_deriv_invers_viscosity[phase_d[phase]]
#define GOR_DERIV_INVERS_VISCOSITY(data)                          data.gor_deriv_invers_viscosity
#define INVERS_VISC_FVF(data, phase_d, phase)                     data.invers_visc_fvf[phase_d[phase]]
#define P_DERIV_INVERS_VISC_FVF(data, phase_d, phase)             data.p_deriv_invers_visc_fvf[phase_d[phase]]
#define GOR_DERIV_INVERS_VISC_FVF(data)                           data.gor_deriv_invers_visc_fvf
#define DENSITY(data, phase_d, phase)                             data.density[phase_d[phase]]
#define P_DERIV_DENSITY(data, phase_d, phase)                     data.p_deriv_density[phase_d[phase]]
#define GOR_DERIV_DENSITY(data)                                   data.gor_deriv_density
#define MOBILITY(data, phase_d, phase)                            data.mobility[phase_d[phase]]
#define P_DERIV_MOBILITY(data, phase_d, phase)                    data.p_deriv_mobility[phase_d[phase]]
#define S_DERIV_MOBILITY(data, phase_d, n_phases, phase1, phase2) data.s_deriv_mobility[phase_d[phase1] * n_phases + phase_d[phase2]]
#define S_DERIV_INVERS_FVF(data, phase_d, phase)                  data.p_deriv_invers_fvf[phase_d[phase]] * data.s_deriv_cap_pressure[phase_d[phase]]
#define S_DERIV_INVERS_VISCOSITY(data, phase_d, phase)            data.p_deriv_invers_viscosity[phase_d[phase]] * data.s_deriv_cap_pressure[phase_d[phase]]
#define S_DERIV_INVERS_VISC_FVF(data, phase_d, phase)             data.p_deriv_invers_visc_fvf[phase_d[phase]] * data.s_deriv_cap_pressure[phase_d[phase]]

#endif
/////////////////////////////////////////////////////////////////////////////

#define RELATIVE_PERM_W RELATIVE_PERM(data, params.phase_d, FI_PHASE_WATER)
#define RELATIVE_PERM_G RELATIVE_PERM(data, params.phase_d, FI_PHASE_GAS)
#define RELATIVE_PERM_O RELATIVE_PERM(data, params.phase_d, FI_PHASE_OIL)

#define S_DERIV_RELATIVE_PERM_WW S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_WATER)
#define S_DERIV_RELATIVE_PERM_WG S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_GAS)
#define S_DERIV_RELATIVE_PERM_WO S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_OIL)
#define S_DERIV_RELATIVE_PERM_GW S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_WATER)
#define S_DERIV_RELATIVE_PERM_GG S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_GAS)
#define S_DERIV_RELATIVE_PERM_GO S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_OIL)
#define S_DERIV_RELATIVE_PERM_OW S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_WATER)
#define S_DERIV_RELATIVE_PERM_OG S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_GAS)
#define S_DERIV_RELATIVE_PERM_OO S_DERIV_RELATIVE_PERM(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_OIL)

#define INVERS_FVF_W INVERS_FVF(data, params.phase_d, FI_PHASE_WATER)
#define INVERS_FVF_G INVERS_FVF(data, params.phase_d, FI_PHASE_GAS)
#define INVERS_FVF_O INVERS_FVF(data, params.phase_d, FI_PHASE_OIL)

#define INVERS_VISCOSITY_W INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_WATER)
#define INVERS_VISCOSITY_G INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_GAS)
#define INVERS_VISCOSITY_O INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_OIL)

#define INVERS_VISC_FVF_W INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_WATER)
#define INVERS_VISC_FVF_G INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_GAS)
#define INVERS_VISC_FVF_O INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_OIL)

#define P_DERIV_INVERS_FVF_W P_DERIV_INVERS_FVF(data, params.phase_d, FI_PHASE_WATER)
#define P_DERIV_INVERS_FVF_O P_DERIV_INVERS_FVF(data, params.phase_d, FI_PHASE_OIL)

#define P_DERIV_INVERS_VISCOSITY_W P_DERIV_INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_WATER)
#define P_DERIV_INVERS_VISCOSITY_O P_DERIV_INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_OIL)

#define S_DERIV_INVERS_VISCOSITY_W S_DERIV_INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_WATER)
#define S_DERIV_INVERS_VISCOSITY_G S_DERIV_INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_GAS)
#define S_DERIV_INVERS_VISCOSITY_O S_DERIV_INVERS_VISCOSITY(data, params.phase_d, FI_PHASE_OIL)

#define S_DERIV_INVERS_VISC_FVF_W S_DERIV_INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_WATER)
#define S_DERIV_INVERS_VISC_FVF_G S_DERIV_INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_GAS)
#define S_DERIV_INVERS_VISC_FVF_O S_DERIV_INVERS_VISC_FVF(data, params.phase_d, FI_PHASE_OIL)

#define P_DERIV_MOBILITY_W P_DERIV_MOBILITY(data, params.phase_d, FI_PHASE_WATER)
#define P_DERIV_MOBILITY_G P_DERIV_MOBILITY(data, params.phase_d, FI_PHASE_GAS)
#define P_DERIV_MOBILITY_O P_DERIV_MOBILITY(data, params.phase_d, FI_PHASE_OIL)

#define S_DERIV_MOBILITY_OW S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_WATER)
#define S_DERIV_MOBILITY_OG S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_GAS)
#define S_DERIV_MOBILITY_OO S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_OIL)

#define S_DERIV_MOBILITY_WW S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_WATER)
#define S_DERIV_MOBILITY_WG S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_GAS)
#define S_DERIV_MOBILITY_WO S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_OIL)

#define S_DERIV_MOBILITY_GW S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_WATER)
#define S_DERIV_MOBILITY_GG S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_GAS)
#define S_DERIV_MOBILITY_GO S_DERIV_MOBILITY(data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_OIL)

#define MOBILITY_W MOBILITY(data, params.phase_d, FI_PHASE_WATER)
#define MOBILITY_G MOBILITY(data, params.phase_d, FI_PHASE_GAS)
#define MOBILITY_O MOBILITY(data, params.phase_d, FI_PHASE_OIL)

#define S_DERIV_CAP_PRESSURE_W S_DERIV_CAP_PRESSURE(data, params.phase_d, FI_PHASE_WATER)
#define S_DERIV_CAP_PRESSURE_G S_DERIV_CAP_PRESSURE(data, params.phase_d, FI_PHASE_GAS)

} // namespace blue_sky


#endif  // #ifndef BS_CALC_MODEL_DATA_ACCESSORS_H_

