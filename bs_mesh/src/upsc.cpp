/** 
 * @file upsc.cpp
 * @brief implementation of hdm upscaling
 * @author Alina Yapparova
 * @version 
 * @date 2012-20-02
 */

#include "bs_mesh_stdafx.h"
#include <iomanip>
#include <iostream>

#include <list>
#include <map>

#include "bs_kernel.h"
#include "upsc.h"

#include "conf.h"

using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

  upsc::upsc (bs_type_ctor_param) 
    {
    
      sp_prop = BS_KERNEL.create_object ("prop");
      if (!sp_prop)
        {
          bs_throw_exception ("Type (prop) not refractered");
        }
      
    }
  upsc::upsc (const upsc& rhs) 
        : bs_refcounter ()
    {
      *this = rhs;
    }
  void
  upsc::init_prop ()
    {

    }
spv_float upsc::upscale_grid ( t_long Nx, t_long Ny, t_long Nz, spv_float zcorn_, spv_uint layers_ )
{   
    t_int k1, k2, new_k;
    t_int new_Nz, new_zcorn_size, layer_size;
    v_uint& layers = *layers_;
    v_uint::iterator lit, lit_next;
    v_float& zcorn = *zcorn_;
    
    new_Nz = layers.size(); 
    new_zcorn_size = 8*Nx*Ny*new_Nz; 
    layer_size = 4*Nx*Ny; 
    
    spv_float new_zcorn = BS_KERNEL.create_object(v_float::bs_type());
    new_zcorn->resize(new_zcorn_size);
    
    new_k = 0;

    // FIXME if layers.size()==1 

    for (lit=layers.begin();lit!=layers.end();lit++)
        {
            k1 = (*lit);
            lit_next = lit + 1;
            //lit_next++;
            if (lit_next!=layers.end())
                k2 = (*lit_next);
            else
                k2 = Nz;
            printf("\n k1 = %d k2 = %d", k1, k2);
            std::copy ( &zcorn[2*k1*layer_size], &zcorn[(2*k1+1)*layer_size], &(*new_zcorn)[(new_k++)*layer_size] );
            std::copy ( &zcorn[(2*k2-1)*layer_size], &zcorn[2*k2*layer_size], &(*new_zcorn)[(new_k++)*layer_size] );
        }

   return new_zcorn; 
}

t_double upsc::calc_sum_dW ( t_long k1, t_long k2, t_long Nx, t_long Ny, 
                                          spv_float vol_, spv_float ntg_, spv_float poro_, spv_float perm_ )
{
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& perm = *perm_;
    
    t_long i, j, n, z1, z2, layer_size;
    t_long index[2];
    t_double p[2];
    t_double dW, sum_dW = 0, sum_vol;
    
    layer_size = Nx*Ny;
    z1 = k1*layer_size;
    z2 = k2*layer_size;
    

    for (j = 0; j < Ny; j++)
      {
        for (i = 0; i < Nx; i++)
          {
            index[0] = i + j * Nx + z1;
            index[1] = i + j * Nx + z2;
            // FIXME
            for (n = 0; n < 2; n++)
                {
                    if (poro[index[n]  != 0])
                        p[n] = ntg[index[n]] * perm[index[n]] / poro[index[n]];
                    else
                        p[n] = 0;
                }
            // FIXME
            sum_vol = (vol[index[0]] + vol[index[1]]);
            if (sum_vol != 0)
                dW = (p[0] - p[1]) * (p[0] - p[1]) * vol[index[0]] * vol[index[1]] / sum_vol;
            else
                dW = 0;
            sum_dW += dW;
          }
      }
    return sum_dW;
}

int upsc::upscale_cubes ( t_long k1, t_long k2, t_long Nx, t_long Ny, 
                          spv_float vol_, spv_float ntg_, spv_float poro_, spv_float swat_ )
{
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& swat = *swat_;
    
    t_long i, j, index[2], z1, z2, layer_size;
    t_double vol_sum, ntg_vol_sum, poro_ntg_vol_sum, swat_poro_ntg_vol_sum; 

    layer_size = Nx*Ny;
    z1 = k1*layer_size;
    z2 = k2*layer_size;

    for (j = 0;j < Ny; j++)
      {
        for (i = 0;i < Nx; i++)
          {
            index[0] = i + j * Nx + z1;
            index[1] = i + j * Nx + z2;

            vol_sum = vol[index[0]] + vol[index[1]];
            ntg_vol_sum = ntg[index[0]]*vol[index[0]] +  ntg[index[1]]*vol[index[1]];
            poro_ntg_vol_sum = poro[index[0]]*ntg[index[0]]*vol[index[0]] + poro[index[1]]*ntg[index[1]]*vol[index[1]];
            swat_poro_ntg_vol_sum =swat[index[0]]*poro[index[0]]*ntg[index[0]]*vol[index[0]] + swat[index[1]]*poro[index[1]]*ntg[index[1]]*vol[index[1]]; 

            // FIXME
            vol[index[0]] = vol_sum;
            if (vol_sum != 0)
                ntg[index[0]] = ntg_vol_sum/vol_sum;
            else
                ntg[index[0]] = 0;
            if (ntg_vol_sum != 0)
                poro[index[0]] = poro_ntg_vol_sum/ntg_vol_sum;
            else
                poro[index[0]] = 0;
            if (swat_poro_ntg_vol_sum != 0)
                swat[index[0]] = swat_poro_ntg_vol_sum/poro_ntg_vol_sum;
            else
                swat[index[0]] = 0;

          }
      }
    return 0;
}

  bp::tuple
  upsc::king_method (t_long Nx, t_long Ny, t_long Nz, t_long Nz_upsc,
                            spv_float vol_, spv_float ntg_, spv_float poro_, 
                            spv_float perm_, spv_float swat_)
  {
    t_int i, k, n, k0, k1, k2, k3;
    t_int layer_size, new_cube_size;
    t_float sum_dW;
    
    //v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& swat = *swat_;

    std::list <t_int> layers;
    std::list <t_int>::iterator lit, tmp_it;

    layer_mmap_t sum_dW2layer;
    layer_mmap_it_t mit;
    layer_map_t layer2mmap_it;
    layer_map_it_t lmit;
    
    //spv_float new_vol = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_ntg = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_poro = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_swat = BS_KERNEL.create_object(v_float::bs_type());
    spv_uint layers_v = BS_KERNEL.create_object(v_uint::bs_type());

    layer_size = Nx*Ny;

    for (i = 0;i < Nz; i++)
        layers.push_back(i);

    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nz*Ny;

    // calculate sum_dW for every pair of adjacent layers
    // and store it in sumdw2layer multimap (sum_dW -> k) sorted in increasing order
    // and in layer2mmapit map (k -> mmap_it) 
    for (k = 0; k < Nz-1; k++)
    {
        sum_dW = calc_sum_dW ( k, k+1,  Nx, Ny, vol_, ntg_, poro_, perm_ );
        mit = sum_dW2layer.insert(std::make_pair(sum_dW, k));
        layer2mmap_it[k] = mit;
    }

    k0 = 0;
    k3 = Nz-1;

    for (n=Nz;n>Nz_upsc;n--)
    {
        // [k0, k1, k2, k3] - indices of 4 consequentive layers
        // [k1, k2]  - indices of 2 layers to be united
        k1 = (*sum_dW2layer.begin()).second;
        lit = find(layers.begin(), layers.end(), k1);
        tmp_it = lit;
        tmp_it ++;
        k2 = (*tmp_it);

        // upscaling of cubes
        upscale_cubes (k1, k2, Nx, Ny, vol_, ntg_, poro_, swat_ );
        
        // erase sum_dW for (k1; k2) pair of layers
        mit = layer2mmap_it[k1];
        sum_dW2layer.erase(mit);
        layer2mmap_it.erase(k1);

        
        // if k2 isn't the last layer
        if (k2 != layers.back())
            {
                // erase sum_dW for (k2; k3) pair of layers
                mit = layer2mmap_it[k2];
                sum_dW2layer.erase(mit);
                layer2mmap_it.erase(k2);
                
                tmp_it ++;
                k3 = (*tmp_it);
                tmp_it --;

                // add sum_dW for (k1; k3) pair of layers
                sum_dW = calc_sum_dW ( k1, k3, Nx, Ny, vol_, ntg_, poro_, perm_ );
                mit = sum_dW2layer.insert ( std::make_pair(sum_dW, k1) );
                layer2mmap_it[k1] = mit;
            }
        
        // delete k2 layer
        layers.erase(tmp_it);

        // if k1 isn't the first layer
        if (k1 != 0)
            {
                tmp_it = lit;
                tmp_it --;
                k0 = (*tmp_it);
                
                // erase sum_dW for (k0; k1) pair of layers
                mit = layer2mmap_it[k0];
                sum_dW2layer.erase(mit);
                layer2mmap_it.erase(k0);

                // add sum_dW for (k0; k1) pair of layers
                sum_dW = calc_sum_dW ( k0, k1, Nx, Ny, vol_, ntg_, poro_, perm_ );
                mit = sum_dW2layer.insert ( std::make_pair(sum_dW, k0) );
                layer2mmap_it[k0] = mit;
            }

    }
    
    // TODO: check that Nz_upsc == layers.size()
    new_cube_size = Nz_upsc*layer_size;

    //new_vol->resize(new_cube_size);
    new_ntg->resize(new_cube_size);
    new_poro->resize(new_cube_size);
    new_swat->resize(new_cube_size);
    layers_v->resize(Nz_upsc);

    i = 0;
    for (lit=layers.begin();lit!=layers.end();lit++)
        {
            k = (*lit);
            //std::copy ( &vol[k*layer_size], &vol[(k+1)*layer_size], &(*new_vol)[i*layer_size] );
            std::copy ( &ntg[k*layer_size], &ntg[(k+1)*layer_size], &(*new_ntg)[i*layer_size] );
            std::copy ( &poro[k*layer_size], &poro[(k+1)*layer_size], &(*new_poro)[i*layer_size] );
            std::copy ( &swat[k*layer_size], &swat[(k+1)*layer_size], &(*new_swat)[i*layer_size] );
            // FIXME
            i++;
        }

    std::copy ( layers.begin(), layers.end(), &(*layers_v)[0] );
    
    return bp::make_tuple (layers_v, /*new_vol,*/ new_ntg, new_poro, new_swat) ;
  }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  upsc::py_str () const
    {
      std::stringstream s;
      s << sp_prop->py_str () << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (upsc);
  BLUE_SKY_TYPE_STD_COPY (upsc);

  BLUE_SKY_TYPE_IMPL(upsc,  upsc_iface, "upsc", "upscaling", "hdm upscaling");

}  // blue_sky namespace
