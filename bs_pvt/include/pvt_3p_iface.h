#ifndef __PVT_3P_IFACE_H
#define __PVT_3P_IFACE_H

/*
#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"
#include "pvt_dummy_iface.h"
*/

namespace blue_sky
  {
/*  
    class BS_API_PLUGIN pvt_base;
    class BS_API_PLUGIN pvt_oil;
    class BS_API_PLUGIN pvt_dead_oil;
    class BS_API_PLUGIN pvt_gas;
    class BS_API_PLUGIN pvt_water;
*/
  /**
   * \brief pvt_base
   */
   
  class pvt_base;
  class pvt_oil;
  class pvt_dead_oil;
  class pvt_gas;
  class pvt_water;
  
   
  class BS_API_PLUGIN pvt_3p_iface : public objbase
    {
      public:
        typedef pvt_base                                  pvt_base_t;               //!< type of base pvt class
        typedef pvt_oil                                   pvt_oil_t;                //!< pvt_oil type
        typedef pvt_dead_oil                              pvt_dead_oil_t;           //!< pvt_dead_oil type
        typedef pvt_gas                                   pvt_gas_t;                //!< pvt_gas type
        typedef pvt_water                                 pvt_water_t;              //!< pvt_water type
      
        typedef smart_ptr <pvt_base_t, true>              sp_pvt_t;                 //!< smart_ptr to pvt_base type
        typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_oil;               //!< smart_ptr to pvt_dead_oil type
        typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_dead_oil;          //!< smart_ptr to pvt_dead_oil type
        typedef smart_ptr <pvt_gas_t, true>               sp_pvt_gas;               //!< smart_ptr to pvt_gas type
        typedef smart_ptr <pvt_water_t, true>             sp_pvt_water;             //!< smart_ptr to pvt_water type
      
        typedef std::vector< sp_pvt_t >                   sp_pvt_array_t;           //!< type for array of pvt_base objects
        typedef std::vector< sp_pvt_oil >                 sp_pvt_oil_array_t;       //!< type for array of pvt_dead_oil objects
        typedef std::vector< sp_pvt_dead_oil >            sp_pvt_dead_oil_array_t;  //!< type for array of pvt_dead_oil objects
        typedef std::vector< sp_pvt_gas >                 sp_pvt_gas_array_t;       //!< type for array of pvt_gas objects
        typedef std::vector< sp_pvt_water >               sp_pvt_water_array_t;     //!< type for array of pvt_water objects

        //typedef smart_ptr <pvt_dummy_iface, true>         sp_pvt_dummy_iface;
        /**
         * \brief destructor
         */
        virtual ~pvt_3p_iface () {}

        virtual BS_SP (pvt_dead_oil) 
        get_pvt_oil (const t_long index_pvt_region) const = 0;
        
        virtual BS_SP (pvt_gas)
        get_pvt_gas (const t_long index_pvt_region) const = 0;
        
        virtual BS_SP (pvt_water)
        get_pvt_water (const t_long index_pvt_region) const = 0;
        
        virtual sp_pvt_dead_oil_array_t &
        get_pvt_oil_array () = 0;
        
        virtual sp_pvt_gas_array_t &
        get_pvt_gas_array () = 0;
        
        virtual sp_pvt_water_array_t &
        get_pvt_water_array () = 0;

		    /*
		    virtual void
		    init_from_pvt(const sp_pvt_dummy_iface &pvt,
		                  bool is_oil, bool is_gas, bool is_water,
					            t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals) = 0;
			  */

        virtual std::list <BS_SP (table_iface)>
        get_tables (t_long index_pvt_region) const = 0;

        virtual void
        init_pvt_arrays (const t_long n_pvt_regions_, 
                         bool is_oil, bool is_gas, bool is_water) = 0;
        
        virtual BS_SP (table_iface)
				get_table (t_long pvt_fluid_type, t_long index_pvt_region) const = 0;	 

        //! get density data to fill in
        virtual spv_float 
        get_density () = 0;

        //! set density to pvt internal data  
        virtual void 
        set_density_to_pvt_internal () = 0;  
				           
				//! build pvt internal tables 
				virtual void 
				build_pvt_internal (t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals) = 0;           
    };

} // namespace blue_sky


#endif // __PVT_3P_IFACE_H