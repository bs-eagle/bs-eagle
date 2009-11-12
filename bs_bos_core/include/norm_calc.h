/**
 *       \file  norm_calc.h
 *      \brief  Calculates norms
 *     \author  Borschuk Oleg
 *       \date  31.01.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef NORM_CALC_H
#define NORM_CALC_H

namespace blue_sky
  {
  namespace norms
    {

    //! enumerate norms
    enum norms_consts
    {
      C_CPV = 0,                      //!< C-norm normalized by cell pore volume
      C_CPV_GAS,                      //!< Gas C-norm normalized by cell pore volume
      C_CPV_WATER,                    //!< Water C-norm normalized by cell pore volume
      C_CPV_OIL,                      //!< Oil C-norm normalized by cell pore volume
      //==========
      C_ACPV,                         //!< C-norm normalized by average cell pore volume
      C_ACPV_GAS,                     //!< Gas C-norm normalized by average cell pore volume
      C_ACPV_WATER,                   //!< Water C-norm normalized by average cell pore volume
      C_ACPV_OIL,                     //!< Oil C-norm normalized by average cell pore volume
      //==========
      L2_CPV,                         //!< L_2-norm normalized by cell pore volume
      L2_CPV_GAS,                     //!< Gas L_2-norm normalized by cell pore volume
      L2_CPV_WATER,                   //!< Water L_2-norm normalized by cell pore volume
      L2_CPV_OIL,                     //!< Oil L_2-norm normalized by cell pore volume
      //==========
      L2_ACPV,                        //!< L_2-norm normalized by average cell pore volume
      L2_ACPV_GAS,                    //!< Gas L_2-norm normalized by average cell pore volume
      L2_ACPV_WATER,                  //!< Water L_2-norm normalized by average cell pore volume
      L2_ACPV_OIL,                    //!< Oil L_2-norm normalized by average cell pore volume
      //==========
      MB_ERR,                         //!< math balans error
      MB_ERR_GAS,                     //!< math balans error for gas
      MB_ERR_WATER,                   //!< math balans error for water
      MB_ERR_OIL,                     //!< math balans error for oil
      //===========
      S_RHS,                          //!< max abs value of sec_rhs

      //========== next constant should be the last one ====
      NORMS_COUNTER
    };

  } // namespace norms

  /*!
   * \brief class for storing different norms
   */
  template <class strategy_t>
  class norms_storage
  {
  public:
    typedef typename strategy_t::item_array_t         item_array_t;
    typedef typename strategy_t::rhs_item_array_t     rhs_item_array_t;
    typedef seq_vector <const char *>                 name_array_t;

  public:
    //! default constructor
    norms_storage ();
    //! default destructor
    ~norms_storage ();
    //! clear all norms
    void clear ();

    // print norms
    //void print (const sp_mesh_grdecl msh);

    norms_storage &operator= (const norms_storage &rhs);

    // --------------------------------------
    // VARIABLES
    // --------------------------------------
  public:
    name_array_t  name;
    item_array_t  val;
    boost::array <int, norms::NORMS_COUNTER> idx;
    boost::array <int, norms::NORMS_COUNTER> p_flag;
  };

} // ns bs

#endif // NORM_CALC_H
