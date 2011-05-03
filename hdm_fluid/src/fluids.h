/**
 * @file fluids.h
 * @brief fluid information storage
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-30
 */
#ifndef FLUIDS_SA24Z34Q

#define FLUIDS_SA24Z34Q

#include "fluids_iface.h"

namespace blue_sky
{
  /**
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN fluids: public fluids_iface
    {
    public:
      typedef h5_pool_iface                     pool_t;
      typedef prop_iface                        prop_t;
      typedef smart_ptr<pool_t, true>           sp_pool_t;
      typedef smart_ptr<prop_t, true>           sp_prop_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (fluids);
    public:

      //! destructor
      virtual ~fluids ()
        {};

       /**
        * @brief
        */
       virtual int init (sp_pool_t sp_pool, sp_prop_t sp_prop)
         {
           // TODO: ADD
           return 0;
         }


       /**
        * @brief return number of phases in the model,
        *        if class non initialized allready or phases non specified
        *        return -1
        */
       virtual t_long get_n_phases () const
         {
           return n_phases;
         }

       /**
        * @brief change the number of PVT regions
        *
        * @param n_regs -- <INPUT> new number of PVT regions
        */
       virtual void set_n_pvt_regs (const t_long n_regs)
         {
           // TODO: ADD
         }

       /**
        * @brief change the number of SCAL regions
        *
        * @param n_regs -- <INPUT> new number of SCAL regions
        */
       virtual void set_n_scal_regs (const t_long n_regs)
         {
           // TODO: ADD
         }

       /**
        * @brief change the number of EQUIL regions
        *
        * @param n_regs -- <INPUT> new number of EQUIL regions
        */
       virtual void set_n_equil_regs (const t_long n_regs)
         {
           // TODO: ADD
         }

       /**
        * @brief make checking of internal state and prepare properties for
        *        calculation process
        *        Call this method before visualizing any fluid property
        *
        * @return 0 if success, < 0 if error occur
        */
       virtual int update ()
         {
           // TODO: ADD
         }


       /**
        * @brief return number of PVT regions
        */
       virtual t_long get_n_pvt_regs () const
         {
           // TODO: ADD
           return 0;
         }

       /**
        * @brief return number of SCAL regions
        */
       virtual t_long get_n_scal_regs () const
         {
           // TODO: ADD
           return 0;
         }

       /**
        * @brief return number of EQUIL regions
        */
       virtual t_long get_n_equil_regs () const
         {
           // TODO: ADD
           return 0;
         }

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
                               const t_long gas)
         {
           // TODO: ADD
           return 0;
         }


#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "FLUIDS\n";

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN

    //______________________________________
    //  VARIABLES
    //______________________________________
    protected:
      t_long n_phases;
      t_long phases;


    };

}//namespace blue_sky
#endif /* end of include guard: FLUIDS_SA24Z34Q */
