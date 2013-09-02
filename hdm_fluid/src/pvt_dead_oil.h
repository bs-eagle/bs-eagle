/**
 * @file pvt_dead_oil.h
 * @brief realization of Dead Oil PVT properties
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-25
 */
#ifndef PVT_DEAD_OIL_V3X97TI

#define PVT_DEAD_OIL_V3X97TI

#include "pvt_oil_iface.h"

#include <vector>
#include <string>

namespace blue_sky
{
  /**
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN pvt_dead_oil: public pvt_oil_iface
    {
    public:

      typedef smart_ptr<prop_iface, true>       sp_prop_t;
      typedef smart_ptr<table_iface, true>      sp_table_t;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (pvt_dead_oil);
    public:

      //! destructor
      virtual ~pvt_dead_oil ()
        {};

      //-----------------------------------------
      //  pvt_iface METHODS
      //-----------------------------------------

       /**
        * @brief check class consistency
        *
        * @return 0 if ok, < 0 error found
        */
       virtual int check () const;

       /**
        * @brief update calculated part of PVT
        *
        * @return 0 if ok, <0 if error occur
        */
       virtual int update ();

       /**
        * @brief return surface density of the phase
        *
        * @param reg    -- <INPUT> region number
        */
       virtual t_double get_surface_density () const
         {
           return surface_density;
         }

       /**
        * @brief return properties for a given region
        *
        * @param reg    -- <INPUT> given PVT region number
        */
       virtual sp_prop_t get_prop ()
         {
           return sp_in_prop;
         }

       /**
        * @brief check and set new property
        *
        * @param new_prop -- <INPUT> given property
        */
       virtual void set_prop (sp_prop_t new_prop)
         {
           sp_in_prop = new_prop;
         }

       /**
        * @brief return table to the input data
        *
        */
       virtual sp_table_t get_table ()
         {
           return sp_in_table;
         }

       /**
        * @brief check and set given table
        *
        * @param new_table -- <INPUT> given table
        */
       virtual void set_table (sp_table_t new_table)
         {
           sp_in_table = new_table;
         }



#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "PVT\n";

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN

    protected:
      void init_prop ();
    //______________________________________
    //  VARIABLES
    //______________________________________
    protected:
      // INPUT DATA
      //! Input property storage
      sp_prop_t         sp_in_prop;
      //! Input PVT table
      sp_table_t        sp_in_table;

      // CALULATED DATA
      //! PVT table used in calculation process
      sp_table_t        sp_calc_table;
      //! surface density
      t_double          surface_density;
      //! molar density
      t_double          molar_density;

    private:


    };

}//namespace blue_sky

#endif /* end of include guard: PVT_DEAD_OIL_V3X97TI */
