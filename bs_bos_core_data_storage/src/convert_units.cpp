/*!
  \file convert_units.cpp
  \brief Contains class convert_units

  Class #convert_units methods.
*/
#include "bs_bos_core_data_storage_stdafx.h"

#include "convert_units.h"
//#include "localization.h"

namespace blue_sky
  {
  
  /*!
    \brief Constructor -- make a new object of class convert_units
    \param units_in_i
    \param units_out_i
  */
  convert_units::convert_units (int units_in_i, int units_out_i)
  {
    set_input_units (units_in_i);
    set_output_units (units_out_i);
  }

  /*!
    \brief set_input_units -- Set input units.
           Return non-zero if unsupported units.
    \param units_in_i
    \return Return non-zero if unsupported units.
  */
  int
  convert_units::set_input_units (int units_in_i)
  {
    switch (units_in_i)
      {
      case UNITS_SI:
      case UNITS_METRIC:
      case UNITS_FIELD:
      case UNITS_LAB:
      case UNITS_INTERNAL:
        units_in = units_in_i;
        input_constants.set_units(units_in_i);
        break;
      default:
        // unsupported:
        units_in = UNITS_INPUT_DEFAULT;
        input_constants.set_units(units_in_i);
        return -1;
      }

    return 0;
  }

  /*!
    \brief set_output_units -- Set output units.\n
           Return non-zero if unsupported units.
    \param units_out_i
    \return Return non-zero if unsupported units.
  */
  int
  convert_units::set_output_units (int units_out_i)
  {
    switch (units_out_i)
      {
      case UNITS_SI:
      case UNITS_METRIC:
      case UNITS_FIELD:
      case UNITS_LAB:
      case UNITS_INTERNAL:
        units_out = units_out_i;
        output_constants.set_units(units_out_i);
        break;
      default:
        // unsupported:
        units_out = UNITS_OUTPUT_DEFAULT;
        output_constants.set_units(units_out_i);
        return -1;
      }
    return YS_SUCCESS;
  }

  /*!
    \brief lenght_mult -- Return multiplier for length conversion
    \return multiplier for length conversion
  */
  double
  convert_units::length_mult (void) const
    {
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL, */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 3.280839895013123, 0., 3.280839895013123,},
        { /* UNITS_FIELD:  */ 0., 0.3048, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 0.3048, 1., 0., 1.,},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief volume_liquid_surface_note -- Return notation for volume of liquid on surface
    \param units
    \return notation for volume of liquid on surface
  */
  const char *
  convert_units::volume_liquid_surface_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3");
      case UNITS_METRIC:
        return  ("sm3");
      case UNITS_FIELD:
        return  ("stb");
      case UNITS_LAB:
        return  ("scc");
      case UNITS_INTERNAL:
        return  ("scf");
      default:
        return 0;
      }
  }

  /*!
    \brief volume_liquid_reservoir_note -- Return notation for volume of liquid in reservoir
    \param units
    \return notation for volume of liquid in reservoir
  */
  const char *
  convert_units::volume_liquid_reservoir_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("rm3");
      case UNITS_METRIC:
        return  ("rm3");
      case UNITS_FIELD:
        return  ("rtb");
      case UNITS_LAB:
        return  ("rcc");
      case UNITS_INTERNAL:
        return  ("rcf");
      default:
        return 0;
      }
  }

  /*!
    \brief volume_gas_surface_note -- Return notation for volume of gas on surface
    \param units
    \return notation for volume of gas on surface
  */
  const char *
  convert_units::volume_gas_surface_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3");
      case UNITS_METRIC:
        return  ("sm3");
      case UNITS_FIELD:
        return  ("Mscf");
      case UNITS_LAB:
        return  ("scc");
      case UNITS_INTERNAL:
        return  ("scf");
      default:
        return 0;
      }
  }

  /*!
    \brief density_mult -- Return multiplier for density conversion
    \return multiplier for density conversion
  */
  double
  convert_units::density_mult (void) const
    {
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL, */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 0.06242796057614461, 0.,
                              0.06242796057614461,},
        { /* UNITS_FIELD:  */ 0., 16.01846337396014, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 16.01846337396014, 1., 0., 1.,},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief pressure_mult -- Return multiplier for pressure conversion
    \return multiplier for pressure conversion
  */
  double
  convert_units::pressure_mult (void) const
    {
#if 1                           // use atm in Metric
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 14.69594877551345, 0., 14.69594877551345,},
        { /* UNITS_FIELD:  */ 0., 0.06804596390987773, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 0.06804596390987773, 1., 0., 1.,},
      };
#else
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 14.50377377302092, 0., 14.50377377302092,},
        { /* UNITS_FIELD:  */ 0., 0.06894757293168361, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 0.06894757293168361, 1., 0., 1.,},
      };
#endif
      return mult[units_in][units_out];
    }

  /*!
    \brief pressure_note -- Return notation for pressure
    \param units
    \return notation for pressure
  */
  const char *
  convert_units::pressure_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("pa");
      case UNITS_METRIC:
        return  ("atm");
      case UNITS_FIELD:
        return  ("psi");
      case UNITS_LAB:
        return  ("psi");
      case UNITS_INTERNAL:
        return  ("psi");
      default:
        return 0;
      }
  }

  /*!
    \brief gas_liquid_mult -- Return multiplier for gas-liquid conversion
    \return multiplier for gas-liquid conversion
  */
  double
  convert_units::gas_liquid_mult (void) const
    {
      /* Metric: m^3/m^3, field: cf/bbl = 1e3/5.614583 Mcf/cf, internal: cf/cf */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 0.005614583333333333, 0., 1.,},
        { /* UNITS_FIELD:  */ 0., 178.1076066790353, 1., 0., 178.1076066790353,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 1., 0.005614583333333333, 0., 1.,},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief gas_liquid_note -- Return notation for gas liquid relation
    \param units
    \return notation for gas liquid relation
  */
  const char *
  convert_units::gas_liquid_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3/sm3");
      case UNITS_METRIC:
        return  ("sm3/sm3");
      case UNITS_FIELD:
        return  ("Mscf/stb");
      case UNITS_LAB:
        return  ("scc/scc");
      case UNITS_INTERNAL:
        return  ("scf/stb");
      default:
        return 0;
      }
  }

  /*!
  \brief liquid_gas_note -- Return notation for liquid gas relation
  \param units
  \return notation for liquid gas relation
  */
  const char *
  convert_units::liquid_gas_note (int units)
    {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3/sm3");
      case UNITS_METRIC:
        return  ("sm3/sm3");
      case UNITS_FIELD:
        return  ("stb/Mscf");
      case UNITS_LAB:
        return  ("scc/scc");
      case UNITS_INTERNAL:
        return  ("stb/scf");
      default:
        return 0;
      }
    }

  /*!
    \brief gas_rate_mult -- Return multiplier for gas surface rate conversion
    \return multiplier for gas surface rate conversion
  */
  double
  convert_units::gas_rate_mult (void) const
    {
      /* Metric: m^3, field: Mcf, internal: cf */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC  , UNITS_FIELD , UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.},
        { /* UNITS_METRIC: */ 0., 1., 0.03531466672148859, 0., 35.31466672148859},
        { /* UNITS_FIELD:  */ 0., 28.31684659200000, 1., 0., 1.e3},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.},
        { /* UNITS_INTERNAL */ 0., 0.028316846592000, 1.e-3, 0., 1.},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief rate_gas_surface_note -- Return notation for rate of gas on surface
    \param units
    \return notation for rate of gas on surface
  */
  const char *
  convert_units::rate_gas_surface_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3/day");
      case UNITS_METRIC:
        return  ("sm3/day");
      case UNITS_FIELD:
        return  ("Mscf/day");
      case UNITS_LAB:
        return  ("scc/hr");
      case UNITS_INTERNAL:
        return  ("scf/day");
      default:
        return 0;
      }
  }

  /*!
    \brief liquid_rate_mult -- Return multiplier for liquid surface rate conversion
    \return multiplier for liquid surface rate conversion
  */
  double
  convert_units::liquid_rate_mult (void) const
    {
      /* Metric: m^3/day = 35.31466 cf/day = 6.289811 bbl/day,
         field: bbl/day = 5.614583 cf/day,
         internal: cf/day */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.},
        { /* UNITS_METRIC: */ 0., 1., 6.289810770432105, 0., 35.31466672148859},
        { /* UNITS_FIELD:  */ 0., 0.158987294928, 1., 0., 5.614583333333333},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.},
        { /* UNITS_INTERNAL */ 0., 0.028316846592, 0.1781076066790353, 0., 1.},
      };
      return mult[units_in][units_out];
    }
  /*!
    \brief gas_liquid_rate_mult -- Return multiplier for gas to liquid rate conversion
    \return multiplier for gas to liquid rate conversion
  */
  double
  convert_units::gas_liquid_rate_mult (void) const
  {
    /* Metric: m^3/m^3, field: bbl/Mcf, internal: cf/cf */
    static double mult[UNITS_TOTAL] = {
      /* In*/
      /* UNITS_SI:     */ 1.,
      /* UNITS_METRIC: */ 1.,
      /* UNITS_FIELD:  */ 0.000178107617253142,
      /* UNITS_LAB:    */ 1.,
      /* UNITS_INTERNAL */ 1.
    };
    return mult[units_in];
  }

  /*!
    \brief rate_liquid_surface_note -- Return notation for rate of liquid on surface
    \param units
    \return notation for rate of liquid on surface
  */
  const char *
  convert_units::rate_liquid_surface_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("sm3/day");
      case UNITS_METRIC:
        return  ("sm3/day");
      case UNITS_FIELD:
        return  ("stb/day");
      case UNITS_LAB:
        return  ("scc/hr");
      case UNITS_INTERNAL:
        return  ("scf/day");
      default:
        return 0;
      }
  }

  /*!
    \brief rate_liquid_reservoir_note -- Return notation for rate of liquid in reservoir
    \param units
    \return notation for rate of liquid in reservoir
  */
  const char *
  convert_units::rate_liquid_reservoir_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("rm3/day");
      case UNITS_METRIC:
        return  ("rm3/day");
      case UNITS_FIELD:
        return  ("rtb/day");
      case UNITS_LAB:
        return  ("rcc/hr");
      case UNITS_INTERNAL:
        return  ("rcf/day");
      default:
        return 0;
      }
  }

  /*!
    \brief density_note -- Return notation for density
    \param units
    \return notation for rate of liquid in reservoir
  */
  const char *
  convert_units::density_note (int units)
  {
    switch (units)
      {
      case UNITS_SI:
        return  ("kg/m3");
      case UNITS_METRIC:
        return  ("kg/m3");
      case UNITS_FIELD:
        return  ("lbl/ft3");
      case UNITS_LAB:
        return  ("kg/m3");
      case UNITS_INTERNAL:
        return  ("rcf/day");
      default:
        return 0;
      }
  }

  /*!
    \brief gas_formation_volume_factor_mult -- Return multiplier for gas formation volume factor conversion
    \return multiplier for gas formation volume factor conversion
  */
  double
  convert_units::gas_formation_volume_factor_mult (void) const
    {
      /* Metric: m^3/m^3, field: bbl/Mcf = 5.614583/1000 cf/cf, internal: cf/cf */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 178.1076066790353, 0., 1.,},
        { /* UNITS_FIELD:  */ 0., 0.005614583333333333, 1., 0.,
                              0.005614583333333333,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 1., 178.1076066790353, 0., 1.,},
      };
      return mult[units_in][units_out];
    }

// Return multiplier for parameter
  double
  convert_units::mult (
    mult_func func_number  //!< enumeration for mul convertion function
  ) const
    {
      double res;
      switch (func_number)
        {
        case length_mult_func:
          res = length_mult ();
          break;
        case density_mult_func:
          res = density_mult ();
          break;
        case pressure_mult_func:
          res = pressure_mult ();
          break;
        case gas_liquid_mult_func:
          res = gas_liquid_mult ();
          break;
        case gas_rate_mult_func:
          res = gas_rate_mult ();
          break;
        case liquid_rate_mult_func:
          res = liquid_rate_mult ();
          break;
        case gas_formation_volume_factor_mult_func:
          res = gas_formation_volume_factor_mult ();
          break;
        case default_mult_func:
        default:
          res = 1.;
        }
      return res;
    }

  /*!
    \brief absolute_temperature_mult -- Return multiplier for absolute temperature conversion
    \return multiplier for absolute temperature conversion
  */
  double
  convert_units::absolute_temperature_mult (void) const
    {
      /*SI: K,
        Metric: K,
        Field: R,
        Lab: K,
        K = R / 1.8
      */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 0.555555555555556, 0., 0.555555555555556,},
        { /* UNITS_FIELD:  */ 0., 1.8, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 1.8, 1., 0., 1.,},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief difference_temperature_mult -- Return multiplier for difference temperature conversion
    \return multiplier for difference temperature conversion
  */
  double
  convert_units::difference_temperature_mult (void) const
    {
      /*SI: C,
        Metric: C,
        Field: F,
        Lab: C,
        F = 32 + 1.8 * C
      */
      static double mult[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 1., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 1., 0.555555555555556, 0., 0.555555555555556,},
        { /* UNITS_FIELD:  */ 0., 1.8, 1., 0., 1.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 1., 0.,},
        { /* UNITS_INTERNAL */ 0., 1.8, 1., 0., 1.,},
      };
      return mult[units_in][units_out];
    }

  /*!
    \brief absolute_temperature_summand -- Return summand for absolute temperature conversion
    \return summand for absolute temperature conversion
  */
  double
  convert_units::absolute_temperature_summand (void) const
    {
      /*SI: K,
        Metric: K,
        Field: R,
        Lab: K,
        K = R / 1.8
      */
      static double summand[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 0., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 0., 0., 0., 0.,},
        { /* UNITS_FIELD:  */ 0., 0., 0., 0., 0.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 0., 0.,},
        { /* UNITS_INTERNAL */ 0., 0., 0., 0., 0.,},
      };
      return summand[units_in][units_out];
    }

  /*!
    \brief difference_temperature_summand -- Return summand for difference temperature conversion
    \return summand for difference temperature conversion
  */
  double
  convert_units::difference_temperature_summand (void) const
    {
      /*SI: C,
        Metric: C,
        Field: F,
        Lab: C,
        F = 32 + 1.8 * C
      */
      static double summand[UNITS_TOTAL][UNITS_TOTAL] =
      {
        /* In / Out         UNITS_SI    , UNITS_METRIC, UNITS_FIELD, UNITS_LAB   , UNITS_INTERNAL */
        { /* UNITS_SI:     */ 0., 0., 0., 0., 0.,},
        { /* UNITS_METRIC: */ 0., 0., -17.777777777777778, 0., -17.777777777777778,},
        { /* UNITS_FIELD:  */ 0., 32., 0., 0., 0.,},
        { /* UNITS_LAB:    */ 0., 0., 0., 0., 0.,},
        { /* UNITS_INTERNAL */ 0., 32., 0., 0., 0.,},
      };
      return summand[units_in][units_out];
    }

  /*!
    \brief absolute_difference_temperature_offset -- Return offset for absolute difference temperature conversion
    \return offset for absolute difference temperature conversion
  */
  double
  convert_units::absolute_difference_temperature_offset (int units) const
    {
      /*SI: K-C,
        Metric: K-C,
        Field: R-F,
        Lab: K-C,
        K = C + 273.15
        R = F + 459.67
      */
      static double offset[UNITS_TOTAL] =
      {
        /* 0 - means currently unsupported */
        /* UNITS_SI:      */ 273.15,
        /* UNITS_METRIC:  */ 273.15,
        /* UNITS_FIELD:   */ 459.67,
        /* UNITS_LAB:     */ 273.15,
        /* UNITS_INTERNAL */ 459.67,
      };
      return offset[units];
    }

} // namespace blue_sky
