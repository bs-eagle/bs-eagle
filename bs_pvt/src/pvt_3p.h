#ifndef PVT_3P_H
#define PVT_3P_H

#include "pvt_3p_iface.h"
#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"
#include "pvt_dummy_iface.h"

namespace blue_sky
  {
    class BS_API_PLUGIN idata;

    class BS_API_PLUGIN pvt_3p : public pvt_3p_iface 
    {
      public: 
  
        typedef idata                                     idata_t;                  //!< idata type
        typedef smart_ptr<idata_t, true>                  sp_idata_t;               //!< smart_ptr to idata type

        typedef pvt_3p_iface::pvt_dead_oil_t              pvt_dead_oil_t;           //!< pvt_dead_oil type
        typedef pvt_3p_iface::pvt_gas_t                   pvt_gas_t;                //!< pvt_gas type
        typedef pvt_3p_iface::pvt_water_t                 pvt_water_t;              //!< pvt_water type
      
        typedef pvt_3p_iface::sp_pvt_dead_oil_array_t     sp_pvt_dead_oil_array_t;  //!< type for array of pvt_dead_oil objects
        typedef pvt_3p_iface::sp_pvt_gas_array_t          sp_pvt_gas_array_t;       //!< type for array of pvt_gas objects
        typedef pvt_3p_iface::sp_pvt_water_array_t        sp_pvt_water_array_t;     //!< type for array of pvt_water objects
    
        typedef BS_SP (pvt_dummy_iface)                   sp_pvt_dummy_iface_t;
        typedef std::vector <sp_pvt_dummy_iface_t>        sp_pvt_dummy_iface_array_t;

        typedef std::vector<BS_SP(table_iface)>           sp_pvt_list_t;
        //typedef std::list<BS_SP(table_iface)>::iterator   sp_pvt_list_iter_t;
        
        virtual 
        ~pvt_3p ()   { } 
        
        BS_SP (pvt_dead_oil) 
        get_pvt_oil (const t_long index_pvt_region) const;
        
        BS_SP (pvt_gas)
        get_pvt_gas (const t_long index_pvt_region) const;
        
        BS_SP (pvt_water)
        get_pvt_water (const t_long index_pvt_region) const;
        
        sp_pvt_dead_oil_array_t &
        get_pvt_oil_array () 
          { return pvt_oil_array; }
        
        sp_pvt_gas_array_t &
        get_pvt_gas_array () 
          { return pvt_gas_array; }
        
        sp_pvt_water_array_t &
        get_pvt_water_array () 
          { return pvt_water_array; }  
        
        /*void
        init_pvt_arrays (const t_long n_pvt_regions_, const sp_idata_t idata_,
                         bool is_oil, bool is_gas, bool is_water, 
                         t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals);*/
        
        /*void
        init_pvt_arrays (const t_long n_pvt_regions_, 
                         const sp_pvt_dummy_iface_array_t &pvt_oil_data,
                         const sp_pvt_dummy_iface_array_t &pvt_gas_data,
                         const sp_pvt_dummy_iface_array_t &pvt_water_data,
                         bool is_oil, bool is_gas, bool is_water, 
                         t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals);*/
                         
        void
        init_pvt_arrays (const t_long n_pvt_regions_, 
                         bool is_oil, bool is_gas, bool is_water);
                         
        void
        init_from_pvt(const sp_pvt_dummy_iface &pvt,
                      bool is_oil, bool is_gas, bool is_water,
										  t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals);
 
 
        void
        init_pvt_arrays (const t_long n_pvt_regions_, 
                         const sp_pvt_dummy_iface_array_t &pvt_data,
                         bool is_oil, bool is_gas, bool is_water, 
                         t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals);
                         
        virtual std::list <BS_SP( table_iface)>
        get_table (t_long index_pvt_region) const;
        
        virtual BS_SP (table_iface)
				get_table (t_long index_pvt_region, t_long pvt_fluid_type) const;	            

        virtual spv_float 
        get_density ()
          { return density; }
          
        //! set density to pvt internal data  
        virtual void 
        set_density_to_pvt_internal ();  

				//! build pvt internal tables 
				virtual void 
				build_pvt_internal (t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals);           
                         
      protected: 
        t_long                       n_pvt_regions;                  //!< number of pvt regions 
        
        sp_pvt_dead_oil_array_t      pvt_oil_array;                  //!< array of pvt_oil objects, length == n_pvt_regions
        sp_pvt_water_array_t         pvt_water_array;                //!< array of pvt_water objects, length == n_pvt_regions
        sp_pvt_gas_array_t           pvt_gas_array;                  //!< array of pvt_gas objects, length == n_pvt_regions
         
        spv_float                    density;                        //!< array of densities, length = n_pvt_regions * FI_PHASE_TOT
        
      public:

        BLUE_SKY_TYPE_DECL (pvt_3p);

    };
    
} // namespace blue_sky
#endif // PVT_3P_H