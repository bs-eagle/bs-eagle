/** 
 * @file upsc.h
 * @brief implementation of hdm upscaling
 * @author Alina Yapparova
 * @version 
 * @date 2012-20-02
 */
#ifndef UPSC_H
#define UPSC_H


#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "upsc_iface.h"
#include "dens_matrix_iface.h"
#include "lsolver_iface.h"

namespace bp = boost::python;

namespace blue_sky
{
  
  class BS_API_PLUGIN upsc : public upsc_iface
    {
    
    public: 

      typedef BS_SP (prop_iface)                      sp_prop_t;
      typedef BS_SP (dens_matrix_iface)               sp_dens_mtx_t;  
      typedef BS_SP (lsolver_iface)                   sp_blu_solver_t;

       
      typedef std::multimap <t_double, t_long> layer_mmap_t;
      typedef std::multimap <t_double, t_long>::iterator layer_mmap_it_t;
      typedef std::map <t_long, std::multimap<t_double,t_long>::iterator> layer_map_t;
      typedef std::map <t_long, std::multimap<t_double,t_long>::iterator>::iterator layer_map_it_t;
      
      typedef grd_ecl::fpoint3d                  fpoint3d_t;
      typedef boost::array <fpoint3d_t, N_PLANE_CORNERS>    plane_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~upsc ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop ()
        {
          return sp_prop;
        }
      
      virtual bp::tuple king_method (t_long nx, t_long ny, t_long nz, t_long nz_upsc,
                                     spv_float vol, spv_float ntg, spv_float poro, spv_float perm); 
      
      virtual spv_float upscale_grid_zcolumn ( t_long Nx, t_long Ny, t_long Nz, spv_float zcorn, spv_uint layers );
      
      //virtual t_int upscale_perm (t_long Nx, t_long Ny, t_long Nz, spv_uint layers, spv_float perm, BS_SP(rs_mesh_iface) mesh);
      
      virtual spv_float upscale_permz_zcolumn (t_long Nx, t_long Ny, t_long Nz, 
                                              spv_uint layers, spv_float permz_, 
                                              spv_uint actnum_, BS_SP(rs_mesh_iface) sp_mesh_iface);

      virtual spv_float upscale_saturation_cube (t_long Nx, t_long Ny, t_long Nz, t_long Nz_upsc,
                                                 spv_uint layers_, spv_float vol_, spv_float ntg_, spv_float poro_, spv_float sat_);

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;
      

#endif //BSPY_EXPORTING_PLUGIN
      
    protected:
      void init_prop ();
      t_float calc_sum_dW ( t_long k1, t_long k2, t_long Nx, t_long Ny,  
                                                spv_float vol, spv_float ntg, 
                                                spv_float poro, spv_float permx );

      t_int upscale_cubes ( t_long k1, t_long k2, t_long Nx, t_long Ny, spv_float vol, spv_float ntg, spv_float poro, spv_float permx );

      t_double solve_pressure_eq (t_long Ny, t_long Nz, t_long i, t_long j, t_long k1, t_long k2, BS_SP(rs_mesh_iface) sp_mesh_iface);
      
      //t_double upsc_permx_zcolumn (t_long Nx, t_long Ny, t_long i, t_long j, t_long k1, t_long k2, spv_float permx_, BS_SP(rs_mesh_iface) sp_mesh_iface);
      
  
      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      sp_prop_t sp_prop;        //!< ptoperties pointer

      BLUE_SKY_TYPE_DECL (upsc);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: UPSC_H */
