/**
 *       \file  rate_control_type.h
 *      \brief  Types of well control
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_RATE_CONTROL_TYPE_H_
#define BS_RATE_CONTROL_TYPE_H_

namespace blue_sky
  {
  namespace wells
    {

    /**
     * \enum  rate_control_type
     * \brief Descripes types of well control
     * */
    enum rate_control_type
    {
      null_control,                         //!< null control
      bhp_control,                          //!< control by bhp
      rate_control,                         //!< control by liquid injection rate
      liquid_rate_control,                  //!< control by liquid rate
      oil_rate_control,                     //!< control by oil rate
      water_rate_control,                   //!< control by water rate
      gas_rate_control,                     //!< control by gas rate
    };

    /**
     * \brief  Converts int value to rate_control_type
     * \param  i Value to convert
     * \return null_control is value is invalid otherwise
     *         element of rate_control_type
     * */
    rate_control_type
    rate_control_cast (int i);

    /**
     * \brief  Converts string value to rate_control_type
     * \param  str Value to convert
     * \return Throws exception if value is invalid,
     *         null_control if value is 'RESV' or
     *         empty string otherwise element of 
     *         rate_control_type
     * */
    rate_control_type
    rate_control_cast (const std::string &str);

    /**
     * \brief  Checks is rate_control_type is bhp control
     * \param  type rate_control_type to check
     * \return True if type is bhp_control
     * */
    bool
    is_bhp_control (rate_control_type type);

    /**
     * \brief  Checks is rate_control_type is rate control
     * \param  type rate_control_type to check
     * \return True is type is not bhp_control and is not 
     *         null_control
     * */
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

