/** 
 * @file frac.h
 * @brief implementation of frac storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef FRAC_YSG17OI0

#define FRAC_YSG17OI0


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "frac_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN frac : public frac_iface
    {
    
    public: 

      typedef BS_SP (prop_iface)                      sp_prop_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~frac ()
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

      BLUE_SKY_TYPE_DECL (frac);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: FRAC_YSG17OI0 */
