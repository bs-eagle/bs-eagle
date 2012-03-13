/** 
 * @file upsc.h
 * @brief implementation of hdm upscaling
 * @author Alina Yapparova
 * @version 
 * @date 2012-20-02
 */
#ifndef UPSC_H

#define UPSC_H


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "upsc_iface.h"



namespace bp = boost::python;

namespace blue_sky
{
  
  class BS_API_PLUGIN upsc : public upsc_iface
    {
    
    public: 

      typedef BS_SP (prop_iface)                      sp_prop_t;

       
      typedef std::multimap <t_double, t_long> layer_mmap_t;
      typedef std::multimap <t_double, t_long>::iterator layer_mmap_it_t;
      typedef std::map <t_long, std::multimap<t_double,t_long>::iterator> layer_map_t;
      typedef std::map <t_long, std::multimap<t_double,t_long>::iterator>::iterator layer_map_it_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~upsc ()
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
      
      virtual bp::tuple king_method (t_long nx, t_long ny, t_long nz, t_long nz_upsc,
                            spv_float vol, spv_float ntg, spv_float poro, spv_float perm, 
                            spv_float swat);
      virtual spv_float upscale_grid ( t_long Nx, t_long Ny, t_long Nz, spv_float zcorn, spv_uint layers );
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
      t_float calc_sum_dW ( t_long k1, t_long k2, t_long Nx, t_long Ny,  
                                                spv_float vol, spv_float ntg, 
                                                spv_float poro, spv_float perm );

      t_int upscale_cubes ( t_long k1, t_long k2, t_long Nx, t_long Ny, 
                      spv_float vol, spv_float ntg, spv_float poro, spv_float swat );
    

      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      sp_prop_t sp_prop;        //!< ptoperties pointer

      BLUE_SKY_TYPE_DECL (upsc);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: UPSC_H */
