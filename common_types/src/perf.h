/** 
 * @file perf.h
 * @brief implementation of perf storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef PERF_YSG17OI0

#define PERF_YSG17OI0


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "perf_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN perf : public perf_iface
    {
    
    public: 

      typedef BS_SP (prop_iface)                      sp_prop_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~perf ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop ()
        {
          return sp_prop;
        }

    public:
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN
      
    protected:
      void init_prop ();

      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      sp_prop_t sp_prop;        //!< ptoperties pointer

      BLUE_SKY_TYPE_DECL (perf);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: perf_YSG17OI0 */
