/*!
  \file convert_units.cpp
  \brief Contains class convert_units

  Class #convert_units methods.
*/
#include "bs_bos_core_data_storage_stdafx.h"

#include "convert_units.h"
#include "localization.h"

namespace blue_sky
  {

  /*!
    \brief set_units -- set units for output constants
    \param units
    \return non-zero if unsupported units
  */
  int
  physical_constants::set_units (int units)
  {
    /* Conversion factors */
    /* 1kg = 2.2049579 lb */
    /* 1 kg / m^3 = 0.06243745 lb / ft^3 */
    /* 1 m = 3.280839 ft */
    /* 1 Atm = 14.696 psi */
    /* 1 Atm = 1.01325 bar */
    /* 1 mD  = 1.062315e-14 ft^2 */
    /* 1 mD  = 9.869233e-16 m^2 */
    /* 1 cP  = 1.450377e-7 psi*s */
    /* 1 Day = 86400 s */
    /* 1 cP  = 1.678677e-12 psi*day */
    /* 1 cP  = 1.142268e-13 atm*day */
    /* 1 Bbl = 5.614583 ft^3 */
    /* 1 m^3 = 35.31466 ft^3 */
    /* 1 m^3 = 6.289811 Bbl  */

    if (units == 0 || units == 3)
      return -1;

    // Darcy_constant --  value of Darcy constant used
    // to convert permeability*length*pressure/viscousity to rate
    static double darcy_constant_array[UNITS_TOTAL] =
    {
      /* 0 - means currently unsupported */
      /* UNITS_SI:      */ 1.,
      /* UNITS_METRIC:  */ 0.00864003,  /* = (m^2/md)*(cp /(atm*day)) */
      /* UNITS_FIELD:   */ 0.00112712,  /* = (bbl/ft^3)*(ft^2/md)*(cp /(psi*day)) */
      /* UNITS_LAB:     */ 3.6,
      /* UNITS_INTERNAL */ 0.006328309, /* = (ft^2/md)*(cp /(psi*day))            */
    };
    darcy_constant = darcy_constant_array[units];

    // gravity_constant --  value of gravity constant
    // is used to convert density*length to pressure

    static double gravity_constant_array[UNITS_TOTAL] =
    {
      /* 0 - means currently unsupported */
      /* UNITS_SI:     */  9.81,      /*  m/s^2        */
      /* UNITS_METRIC: */  0.000096816838, /* (atm * m^2)/kg  */
      /* UNITS_FIELD:  */  0.006944,   /* (psi * ft^2)/lb */
      /* UNITS_LAB:    */  0.000096816838,
      /* UNITS_INTERNAL */ 0.00694,   /* (psi * ft^2)/lb */
    };
    gravity_constant = gravity_constant_array[units];


    //atrmospheric_pressure --  value of atmospheric pressure

    static double atmospheric_pressure_array[UNITS_TOTAL] =
    {
      /* 0 - means currently unsupported */
      /* UNITS_SI:     */  101325.35,   /*  Pa/atm */
      /* UNITS_METRIC: */  1.0,
      /* UNITS_FIELD:  */  14.695948,   /* psi/atm */
      /* UNITS_LAB:    */  1.0,
      /* UNITS_INTERNAL */ 14.695948,   /* psi/atm */
    };
    atmospheric_pressure = atmospheric_pressure_array[units];

    static double jfunc_constant_array[UNITS_TOTAL] =
    {
      /* 0 - means currently unsupported */
      /* UNITS_SI:     */  3.2120816e+4,
      /* UNITS_METRIC: */  0.314153,
      /* UNITS_FIELD:  */  4.61678,
      /* UNITS_LAB:    */  0.314153,
      /* UNITS_INTERNAL */ 4.61678,
    };
    jfunc_constant = jfunc_constant_array[units];

    // Default injection limit for BHP is 700 atm
    default_injection_bhp_limit = atmospheric_pressure * 700;

    // Minimal possible average reservoir pressure in model 10 atm
    default_minimal_average_pressure = atmospheric_pressure * 10;

    // Default production limit for BHP is 1 atm
    default_production_bhp_limit = atmospheric_pressure * 1.000001;

    static double cubic_meter_per_day_array[UNITS_TOTAL] =
    {
      /* 0 - means currently unsupported */
      /* UNITS_SI:     */  1.1574074E-05,  /*   (m^3/s)/ (m^3/day) */
      /* UNITS_METRIC: */  1.0,
      /* UNITS_FIELD:  */  6.289811,     /*  (bbl/day)/ (m^3/day) */
      /* UNITS_LAB:    */  0,
      /* UNITS_INTERNAL */ 35.31466,     /*  (scf/day)/ (m^3/day) */
    };

    // Default limit for LIQUID RATE is 1000000 m^3 /day
    default_liquid_rate_limit = cubic_meter_per_day_array[units] * 1e6;

    // Default gas rate is 50000000 m^3/day
    default_gas_rate_limit = cubic_meter_per_day_array[units] * 50e6;

#if COMPOSITIONAL
    // gas_constant --  value of gas constant used
    static double gas_constant_array[UNITS_TOTAL] =
    {
      /* UNITS_SI:      */ 8.31,
      /* UNITS_METRIC:  */ 0.083143,
      /* UNITS_FIELD:   */ 10.732,
      /* UNITS_LAB:     */ 820.55776,
      /* UNITS_INTERNAL */ 0,
    };
    gas_constant = gas_constant_array[units];
#endif

    return 0;
  }

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
        return GET_TEXT ("sm3");
      case UNITS_METRIC:
        return GET_TEXT ("sm3");
      case UNITS_FIELD:
        return GET_TEXT ("stb");
      case UNITS_LAB:
        return GET_TEXT ("scc");
      case UNITS_INTERNAL:
        return GET_TEXT ("scf");
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
        return GET_TEXT ("rm3");
      case UNITS_METRIC:
        return GET_TEXT ("rm3");
      case UNITS_FIELD:
        return GET_TEXT ("rtb");
      case UNITS_LAB:
        return GET_TEXT ("rcc");
      case UNITS_INTERNAL:
        return GET_TEXT ("rcf");
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
        return GET_TEXT ("sm3");
      case UNITS_METRIC:
        return GET_TEXT ("sm3");
      case UNITS_FIELD:
        return GET_TEXT ("Mscf");
      case UNITS_LAB:
        return GET_TEXT ("scc");
      case UNITS_INTERNAL:
        return GET_TEXT ("scf");
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
        return GET_TEXT ("pa");
      case UNITS_METRIC:
        return GET_TEXT ("atm");
      case UNITS_FIELD:
        return GET_TEXT ("psi");
      case UNITS_LAB:
        return GET_TEXT ("psi");
      case UNITS_INTERNAL:
        return GET_TEXT ("psi");
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
        return GET_TEXT ("sm3/sm3");
      case UNITS_METRIC:
        return GET_TEXT ("sm3/sm3");
      case UNITS_FIELD:
        return GET_TEXT ("Mscf/stb");
      case UNITS_LAB:
        return GET_TEXT ("scc/scc");
      case UNITS_INTERNAL:
        return GET_TEXT ("scf/stb");
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
        return GET_TEXT ("sm3/sm3");
      case UNITS_METRIC:
        return GET_TEXT ("sm3/sm3");
      case UNITS_FIELD:
        return GET_TEXT ("stb/Mscf");
      case UNITS_LAB:
        return GET_TEXT ("scc/scc");
      case UNITS_INTERNAL:
        return GET_TEXT ("stb/scf");
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
        return GET_TEXT ("sm3/day");
      case UNITS_METRIC:
        return GET_TEXT ("sm3/day");
      case UNITS_FIELD:
        return GET_TEXT ("Mscf/day");
      case UNITS_LAB:
        return GET_TEXT ("scc/hr");
      case UNITS_INTERNAL:
        return GET_TEXT ("scf/day");
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
        return GET_TEXT ("sm3/day");
      case UNITS_METRIC:
        return GET_TEXT ("sm3/day");
      case UNITS_FIELD:
        return GET_TEXT ("stb/day");
      case UNITS_LAB:
        return GET_TEXT ("scc/hr");
      case UNITS_INTERNAL:
        return GET_TEXT ("scf/day");
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
        return GET_TEXT ("rm3/day");
      case UNITS_METRIC:
        return GET_TEXT ("rm3/day");
      case UNITS_FIELD:
        return GET_TEXT ("rtb/day");
      case UNITS_LAB:
        return GET_TEXT ("rcc/hr");
      case UNITS_INTERNAL:
        return GET_TEXT ("rcf/day");
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
        return GET_TEXT ("kg/m3");
      case UNITS_METRIC:
        return GET_TEXT ("kg/m3");
      case UNITS_FIELD:
        return GET_TEXT ("lbl/ft3");
      case UNITS_LAB:
        return GET_TEXT ("kg/m3");
      case UNITS_INTERNAL:
        return GET_TEXT ("rcf/day");
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
