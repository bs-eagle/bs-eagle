#ifndef ROCK_GRID_H
#define ROCK_GRID_H

#include "convert_units.h"
#include "conf.h"

namespace blue_sky
{
  class BS_API_PLUGIN fi_params;

   class BS_API_PLUGIN idata;

  
  class BS_API_PLUGIN rock_grid : public objbase
  {
  public:
    // typedefs
    //typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
    
    typedef idata                                idata_t;

    typedef smart_ptr< fi_params, true>         sp_fi_params;
    typedef smart_ptr< idata_t, true>           sp_idata;

    //! dtor
    virtual ~rock_grid();

    //! init rock_grid data
    void init (const sp_idata &input_data, t_long n_els, t_int n_pvt_regs);

    //! calculate geometric plane transmissibility
    int init_planes_trans (t_long cells_count, const spv_float mesh_volumes, const sp_fi_params &ts_params, physical_constants &internal_consts);

    //! init data
    int init_data (t_long cells_count, const spv_long index_map, const sp_idata &input_data);

    
    //protected:
    stdv_float  net_to_gros,       //!< (n_elements) coefficient of Net to Gross Thickness Ratios
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
    t_long n_elements;         //!< number of elements
    t_int n_pvt_regions;      //!< number of pvt regions
    t_long n_planes;           //!< number of planes
    t_long n_boundary_planes;  //!< number of boundary planes

    //! blue-sky class declaration
    BLUE_SKY_TYPE_DECL(rock_grid)
  };
}

#endif // ROCK_GRID_H
