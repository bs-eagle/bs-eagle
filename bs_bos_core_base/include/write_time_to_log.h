/**
 * \file write_time_to_log.h
 * \brief Helper to write time diff into log file
 * \author Miryanov Sergey
 * \date 15.04.2008
 */

#ifndef BS_WRITE_TIME_TO_LOG_H_
#define BS_WRITE_TIME_TO_LOG_H_

#include "bos_report.h"
#include <stdio.h>

namespace blue_sky
  {

  struct write_time_to_log
    {
      clock_t t;
      const std::string str;
      const std::string postfix;

      write_time_to_log (const std::string &str, const std::string &postfix)
      : str (str)
      , postfix (postfix)
      {
        //BOSOUT (section::app_info, level::medium) << "---" << str << "_" << postfix << bs_end;
        t = clock ();
      }

      ~write_time_to_log ()
      {
        clock_t t2 = clock ();
        double diff = double (t2 - t) / double (CLOCKS_PER_SEC);
        if (postfix != "")
          BOSOUT (section::app_info, level::medium) << str << "_" << postfix << " --- " << diff << "s" << bs_end;
        else
          BOSOUT (section::app_info, level::medium) << str << " --- " << diff << "s" << bs_end;
      }
    };

}

#endif // #ifndef BS_WRITE_TIME_TO_LOG_H_
