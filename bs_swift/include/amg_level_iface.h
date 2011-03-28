/** 
 * @file amg_level_iface.h
 * @brief interface for geting all data from AMG level
 * @date 2009-12-05
 */
#ifndef __AMG_LEVEL_IFACE_H
#define __AMG_LEVEL_IFACE_H

#include "bs_assert.h"
#include "bs_tree.h"

#include "bcsr_amg_matrix_iface.h"

namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  class amg_level_iface: public bs_node
    {
    public:
      typedef amg_level_iface                                   this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;

      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      //! sp for bcsr matrix
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------

      /** 
       * @return reference to the main matrix of the next level
       */
      const bcsr_t &get_next_level_matrix () const = 0;  

      /** 
       * @return reference to the P(prolongation) matrix 
       */
      const bcsr_t &get_prolongation_matrix () const = 0;  


      /** 
       * @return reference to the right hand side vector for the next level 
       */
      const fp_vector_type_t &get_next_level_rhs () const = 0;

      /** 
       * @return reference to the solution vector for the next level 
       */
      const fp_vector_type_t &get_next_level_solution () const = 0;


    public:


#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    public:
      //! destructor
      virtual ~amg_level_iface ()
        {}

    };

}//namespace blue_sky

#endif //__AMG_LEVEL_IFACE_H


