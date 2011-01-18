/*!
  \file convert_units.h
  \brief Class #convert_units header file
*/
#ifndef CONVERT_UNITS_H_
#define CONVERT_UNITS_H_

namespace blue_sky
  {

#define UNITS_SI        0
#define UNITS_METRIC    1
#define UNITS_FIELD     2
#define UNITS_LAB       3
  //! Used in project units: = field, but cubic foots instead of barrels for rates
#define UNITS_INTERNAL  4

  //! Number of supported systems
#define UNITS_TOTAL     5
  //! Default input system
#define UNITS_INPUT_DEFAULT     UNITS_METRIC
  //! Default output system
#define UNITS_OUTPUT_DEFAULT    UNITS_METRIC
  //! Default internal project units
#define UNITS_INTERNAL_DEFAULT  UNITS_METRIC
  //#define UNITS_INTERNAL_DEFAULT  UNITS_INTERNAL

  //! \brief enumeration for mul convertion function
  enum mult_func
  {
    length_mult_func                       = 0,     //!< number of function for length convertion
    density_mult_func,                              //!< number of function for density convertion
    pressure_mult_func,                             //!< number of function for pressure convertion
    gas_liquid_mult_func,                           //!< number of function for gas liquid convertion
    gas_rate_mult_func,                             //!< number of function for gas rate convertion
    liquid_rate_mult_func,                          //!< number of function for liquid rate convertion
    gas_formation_volume_factor_mult_func,          //!< number of function for gas fvf convertion
    default_mult_func                               //!< no multiplier - mult 1
  };
  
  /*!
    \class physical_constants
    \ingroup KeywordLanguage
    \brief Class used to output physical constants for specified unit system
  */

  class BS_API_PLUGIN physical_constants
    {
      public:
      
      physical_constants (int units = UNITS_INTERNAL_DEFAULT)
      {
        set_units(units);
      }
  
      // Set  units. Return non-zero if unsupported units.
      int set_units (int units)
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

      double darcy_constant;                          //!< Return value of Darcy constant
      double gravity_constant;                        //!< Return value of gravity constant
      double atmospheric_pressure;                    //!< Return value of atmospheric pressure
      double jfunc_constant;                          //!< return value of jfunc multiplier constant
      double default_injection_bhp_limit;             //!< Default injection limit for BHP is 1000 atm
      double default_production_bhp_limit;            //!< Default production limit for BHP is 1 atm
      double default_minimal_average_pressure;        //!< Minimal allowable average reservoir pressure in model -- 10 atm
      double default_liquid_rate_limit;               //!< Default limit for LIQUID RATE is 1000000 m^3 /day
      double default_gas_rate_limit;                  //!< Default gas rate is 50000000 m^3/day
      double gas_constant;                            //!< Return value of gas constant
    };

  /*!
    \brief set_units -- set units for output constants
    \param units
    \return non-zero if unsupported units
  */

  /*!
    \class convert_units
    \ingroup KeywordLanguage
    \brief Class used for convertion different types of units
           Class convert_units is a part of UfaSolver keyword input language
  */

  class BS_API_PLUGIN convert_units
    {
    public:
      convert_units (int units_in_i = UNITS_INPUT_DEFAULT,
                     int units_out_i = UNITS_OUTPUT_DEFAULT);

      physical_constants output_constants;
      physical_constants input_constants;
      // Set input units. Return non-zero if unsupported units.
      int set_input_units (int units_in_i);
      int get_input_units () const  //! Return input units.
        {
          return units_in;
        };
      // Set output units. Return non-zero if unsupported units.
      int set_output_units (int output_in_i);
      int get_output_units () const //! Return output units.
        {
          return units_out;
        };
      // Return multiplier for parameter
      double mult (mult_func func_number) const;
      // Return multiplier for length conversion
      double length_mult (void) const;
      // Return multiplier for density conversion
      double density_mult (void) const;
      // Return multiplier for pressure conversion
      double pressure_mult (void) const;
      // Return multiplier for gas-liquid conversion
      double gas_liquid_mult (void) const;
      // Return multiplier for gas to liquid rate conversion
      double gas_liquid_rate_mult (void) const;
      // Return multiplier for gas surface rate conversion
      double gas_rate_mult (void) const;
      // Return multiplier for liquid surface rate conversion
      double liquid_rate_mult (void) const;
      // Return multiplier for gas formation volume factor conversion
      double gas_formation_volume_factor_mult (void) const;
      // Return multiplier for absolute temperature conversion
      double absolute_temperature_mult (void) const;
      // Return multiplier for difference temperature conversion
      double difference_temperature_mult (void) const;
      // Return summand for absolute temperature conversion
      double absolute_temperature_summand (void) const;
      // Return summand for difference temperature conversion
      double difference_temperature_summand (void) const;
      // Return offset for absolute difference temperature conversion
      double absolute_difference_temperature_offset (int units) const;

      // Return notation for pressure
      static const char *pressure_note (int units);
      // Return notation for volume of liquid on surface
      static const char *volume_liquid_surface_note (int units);
      // Return notation for volume of liquid in reservoir
      static const char *volume_liquid_reservoir_note (int units);
      // Return notation for volume of gas on surface
      static const char *volume_gas_surface_note (int units);
      // Return notation for rate of liquid on surface
      static const char *rate_liquid_surface_note (int units);
      // Return notation for rate of liquid in reservoir
      static const char *rate_liquid_reservoir_note (int units);
      // Return notation for rate of gas on surface
      static const char *rate_gas_surface_note (int units);
      // Return notation for gas liquid relation
      static const char *gas_liquid_note (int units);
      // Return notation for liquid gas relation
      static const char *liquid_gas_note (int units);
      // Return notation for density
      static const char *density_note (int units);

    private:
      int units_in;                 //!< Input units
      int units_out;                //!< Output units
    };

} // namespace blue_sky

#endif // CONVERT_UNITS_H_
