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
      int set_units (int units);

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
