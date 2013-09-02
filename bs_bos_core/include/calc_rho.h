/**
 *       \file  calc_rho.h
 *      \brief  Density (rho) calculation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  26.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be renamed to calc_density.h
 * */
#ifndef BS_CALC_RHO_H_
#define BS_CALC_RHO_H_

namespace blue_sky
  {
  class calc_model;
  class well;

  /**
   * \class calc_rho_base
   * \brief Base class for density (rho) calculation
   * \todo  Should be removed to calc_density_iface
   * */
  class calc_rho_base : public objbase
    {
    public:
      typedef t_long                          index_t;
      typedef t_double                        item_t;
      typedef spv_double                      item_array_t;

      typedef calc_model                      calc_model_t;
      typedef well                            well_t;
      typedef rs_mesh_iface                   mesh_iface_t;

      typedef smart_ptr <calc_model_t, true>  sp_calc_model_t;
      typedef smart_ptr <well_t, true>        sp_well_t;
      typedef smart_ptr <mesh_iface_t, true>  sp_mesh_iface_t;

    public:
      /**
       * \brief  calc_rho_base dtor
       * */
      virtual ~calc_rho_base () {}

      /**
       * \brief  For each well perforation calculates density value
       * \param  well well to calculate perforation density value
       * \param  calc_model
       * \param  mesh
       * */
      virtual void 
      calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const = 0;
    };

  /**
   * \class calc_total_average_rho
   * \brief Calculates total average density (rho)
   * */
  class calc_total_average_rho : public calc_rho_base
    {
    public:
      typedef t_long                                index_t;
      typedef t_double                              item_t;
      typedef spv_double                            item_array_t;

      typedef calc_model                            calc_model_t;
      typedef well                                  well_t;
      typedef rs_mesh_iface                         mesh_iface_t;

      typedef smart_ptr <calc_model_t, true>        sp_calc_model_t;
      typedef smart_ptr <well_t, true>              sp_well_t;
      typedef smart_ptr <mesh_iface_t, true>        sp_mesh_iface_t;

    public:

      /**
       * \brief  For each well perforation calculates density value as a total average density
       * \param  well well to calculate perforation density value
       * \param  calc_model
       * \param  mesh
       * */
      void 
      calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (calc_total_average_rho);
    };

  /**
   * \brief  Registers calc_rho and calc_total_averga_rho types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return true if types registered successfully
   * */
  bool
  calc_rho_register_types (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky

#endif  // #ifndef BS_CALC_RHO_H_

