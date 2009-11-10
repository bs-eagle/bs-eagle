#ifndef ROCK_GRID_H
#define ROCK_GRID_H

#include "convert_units.h"

namespace blue_sky
{
  class BS_API_PLUGIN fi_params;

  class BS_API_PLUGIN idata;

  template <class strategy_t>
  class BS_API_PLUGIN rock_grid : public objbase
  {
  public:
    // typedefs
    typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
    typedef typename strategy_t::item_t         item_t;         ///< short name to array item type
    typedef typename strategy_t::item_array_t   item_array_t;
    typedef typename strategy_t::index_t        index_t;        ///< short name to matrix's index type
    typedef typename strategy_t::index_array_t  index_array_t;

    typedef idata                               idata_t;

    typedef smart_ptr< fi_params, true>         sp_fi_params;
    typedef smart_ptr< idata_t, true>           sp_idata;

    //! dtor
    virtual ~rock_grid();

    //! init rock_grid data
    void init (const sp_idata &input_data, index_t n_els, index_t n_pvt_regs);

    //! calculate geometric plane transmissibility
    int init_planes_trans (index_t cells_count, const item_array_t &mesh_volumes, const sp_fi_params &ts_params, physical_constants &internal_consts);

    //! init data
    int init_data (index_t cells_count, const index_array_t &index_map, const sp_idata &input_data);

    //! blue-sky class declaration
    BLUE_SKY_TYPE_DECL(rock_grid)

    //protected:
    item_array_t  net_to_gros,       //!< (n_elements) coefficient of Net to Gross Thickness Ratios
    porosity_p_ref,    //!< (n_elements) porosity by pressure p_ref
    permeability,      //!< (PLANE_ORIENTATION_TOTAL * n_elements) permeability for x, y, z directions
    comp_const,        //!< (n_pvt_regions) compress constant
    comp_ref_pressure, //!< (n_pvt_regions) compress ref pressure

    // phisical properties calculed once
    //---------------------------------------------------------------------------------------------------------
    planes_trans,      //!< (n_planes) transmissibility for all planes
    trans_mult,        //!< (3 * n_elements) transmissibility multipliers
    volume,            //!< (n_elements) volume of all elements
    multpv;            //!< (n_elements) poro volume multipliers

  protected:
    index_t n_elements;         //!< number of elements
    index_t n_pvt_regions;      //!< number of pvt regions
    index_t n_planes;           //!< number of planes
    index_t n_boundary_planes;  //!< number of boundary planes
  };
}

#endif // ROCK_GRID_H
