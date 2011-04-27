/** 
 * @file fluids_iface.h
 * @brief Interface for Fluids storage (PVT, SCAL)
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-29
 */
#ifndef FLUIDS_IFACE_YCX3MW5G

#define FLUIDS_IFACE_YCX3MW5G


#include "bs_assert.h"
#include "bs_tree.h"
#include "smart_ptr.h"
#include "bs_array.h"
#include "conf.h"
#include "pool_iface.h"
#include "prop_iface.h"

#include <string>

namespace blue_sky
{
  /** 
   * @brief PVT interface
   */
  class fluids_iface: public objbase
  {
    public:
      typedef h5_pool_iface                     pool_t;
      typedef prop_iface                        prop_t;
      typedef smart_ptr<pool_t, true>           sp_pool_t;
      typedef smart_ptr<prop_t, true>           sp_prop_t;

      /** 
       * @brief destructor
       */
       virtual ~fluids_iface () {}

       /** 
        * @brief Initialize by given pool and property storage 
        * 
        * @param sp_pool        -- <INPUT> pool
        * @param sp_prop        -- <INPUT> prop
        * 
        * @return 0 if ok, < 0 if error occur
        */
       virtual int init (sp_pool_t sp_pool, sp_prop_t sp_prop) = 0;

       /** 
        * @brief return number of phases in the model, 
        *        if class non initialized allready or phases non specified
        *        return -1
        */
       virtual t_long get_n_phases () const = 0;

       /** 
        * @brief change the number of PVT regions
        * 
        * @param n_regs -- <INPUT> new number of PVT regions
        */
       virtual void set_n_pvt_regs (const t_long n_regs) = 0;

       /** 
        * @brief change the number of SCAL regions
        * 
        * @param n_regs -- <INPUT> new number of SCAL regions
        */
       virtual void set_n_scal_regs (const t_long n_regs) = 0;

       /** 
        * @brief change the number of EQUIL regions
        * 
        * @param n_regs -- <INPUT> new number of EQUIL regions
        */
       virtual void set_n_equil_regs (const t_long n_regs) = 0;

       /** 
        * @brief make checking of internal state and prepare properties for 
        *        calculation process
        *        Call this method before visualizing any fluid property
        * 
        * @return 0 if success, < 0 if error occur
        */
       virtual int update () = 0;


       /** 
        * @brief return number of PVT regions
        */
       virtual t_long get_n_pvt_regs () const = 0;

       /** 
        * @brief return number of SCAL regions
        */
       virtual t_long get_n_scal_regs () const = 0;

       /** 
        * @brief return number of EQUIL regions
        */
       virtual t_long get_n_equil_regs () const = 0;

       /** 
        * @brief set phases 
        * 
        * @param water  -- <INPUT> 0 -- with out water, 1 -- with water
        * @param oil    -- <INPUT> 0 -- with out oil, 1 -- oil with out desolvet gas,
        *                          2 -- oil with desolved gas
        * @param gas    -- <INPUT> 0 -- with out gas, 1 -- with dry gas only,
        *                          2 -- gas and vaporised oil
        * 
        * @return 0 if ok, < 0 if error occur
        */
       virtual int set_phases (const t_long water, 
                               const t_long oil, 
                               const t_long gas) = 0;



#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

  };
};

#endif /* end of include guard: FLUIDS_IFACE_YCX3MW5G */

