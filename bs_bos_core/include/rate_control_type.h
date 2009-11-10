/**
 * \file rate_control_type.h
 * \brief
 * \author Sergey Miryanov
 * \date 24.11.2008
 * */
#ifndef BS_RATE_CONTROL_TYPE_H_
#define BS_RATE_CONTROL_TYPE_H_

namespace blue_sky
  {
  namespace wells
    {

    enum rate_control_type
    {
      null_control,                         //!< null control
      bhp_control,                          //!< control by bhp
      rate_control,													//!< control by liquid injection rate
      liquid_rate_control,                  //!< control by liquid rate
      oil_rate_control,											//!< control by oil rate
      water_rate_control,										//!< control by water rate
      gas_rate_control,											//!< control by gas rate
    };

    rate_control_type
    rate_control_cast (int i);              //!< cast int value to rate_control_type; if value is invalid it returns null_control
    rate_control_type
    rate_control_cast (const std::string &str);

    bool
    is_bhp_control (rate_control_type type);
    bool
    is_rate_control (rate_control_type type);

    inline rate_control_type
    rate_control_cast (int i)
    {
      if (i <= null_control && i > gas_rate_control)
        {
          BS_ASSERT (false && "UNSUPPORTED VALUE") (i);
          return null_control;
        }

      return static_cast <rate_control_type> (i);
    }
    inline rate_control_type
    rate_control_cast (const std::string &str)
    {
      if (str == "BHP")
        return bhp_control;
      else if (str == "RATE")
        return rate_control;
      else if (str == "LRAT")
        return liquid_rate_control;
      else if (str == "ORAT")
        return oil_rate_control;
      else if (str == "WRAT")
        return water_rate_control;
      else if (str == "GRAT")
        return gas_rate_control;
      else if (str == "RESV")
        {
          BS_ASSERT (false && "Unsupported rate_control type") (str);
          return null_control;
        }
      else if (str == "")
        return null_control;
      else
        {
          BS_ASSERT (false && "Unsuported rate_control type") (str);
          throw bs_exception ("rate_control_cast", "Unsupported rate_control type");
        }
    }

    inline bool
    is_bhp_control (rate_control_type type)
    {
      return type == bhp_control;
    }
    inline bool
    is_rate_control (rate_control_type type)
    {
      return type  != bhp_control && type != null_control;
    }


  } // namespace wells
} // namespce blue_sky


#endif  // #ifndef BS_RATE_CONTROL_TYPE_

