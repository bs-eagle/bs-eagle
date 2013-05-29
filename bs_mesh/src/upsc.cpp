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
#include "bs_misc.h"

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

spv_float upsc::upscale_grid_zcolumn ( t_long Nx, t_long Ny, t_long Nz, spv_float zcorn_, spv_uint layers_ )
{   
    t_int k1, k2, new_k;
    t_long new_Nz, new_zcorn_size, layer_size;
    v_uint& layers = *layers_;
    v_uint::iterator lit, lit_next;
    v_float& zcorn = *zcorn_;
    
    new_Nz = layers.size(); 
    new_zcorn_size = 8*Nx*Ny*new_Nz; 
    layer_size = 4*Nx*Ny; 
    
    spv_float new_zcorn = BS_KERNEL.create_object(v_float::bs_type());
    new_zcorn->resize(new_zcorn_size);
    
    new_k = 0;

    for (lit=layers.begin();lit!=layers.end();lit++)
        {
            k1 = (*lit);
            lit_next = lit + 1;
            if (lit_next!=layers.end())
                k2 = (*lit_next);
            else
                k2 = Nz;
            std::copy ( &zcorn[2*k1*layer_size], &zcorn[(2*k1+1)*layer_size], &(*new_zcorn)[(new_k++)*layer_size] );
            std::copy ( &zcorn[(2*k2-1)*layer_size], &zcorn[2*k2*layer_size], &(*new_zcorn)[(new_k++)*layer_size] );
        }

   return new_zcorn; 
}

bp::tuple upsc::upscale_grid ( t_long Nx, t_long Ny, t_long Nz, t_long ux, t_long uy, 
                               spv_float coord_, spv_float zcorn_, spv_uint layers_ )
{   
    t_int i, j, n, k1, k2, k, ind;
    t_long new_Nx, new_Ny, new_Nz, new_coord_size, new_zcorn_size, layer_size;
    v_uint& layers = *layers_;
    v_uint::iterator lit, lit_next;
    v_float& coord = *coord_;
    v_float& zcorn = *zcorn_;
    
    new_Nx = ceil(double(Nx)/double(ux));
    new_Ny = ceil(double(Ny)/double(uy));
    new_Nz = layers.size(); 
    new_coord_size = 6*(new_Nx+1)*(new_Ny+1);
    new_zcorn_size = 8*new_Nx*new_Ny*new_Nz; 
    layer_size = 4*Nx*Ny;
    
    spv_float new_coord = BS_KERNEL.create_object(v_float::bs_type());
    new_coord->resize(new_coord_size);
    spv_float new_zcorn = BS_KERNEL.create_object(v_float::bs_type());
    new_zcorn->resize(new_zcorn_size);

    n = 0;
    for (j=0; j<new_Ny; j++)
        {
            for (i=0; i<new_Nx; i++)
                {
                    std::copy ( &coord[6*(uy*j*(Nx+1)+ux*i)], &coord[6*(uy*j*(Nx+1)+ux*i+1)], &(*new_coord)[n] );
                    n += 6;
                }
            std::copy ( &coord[6*(uy*j*(Nx+1)+Nx)], &coord[6*(uy*j*(Nx+1)+Nx+1)], &(*new_coord)[n] );
            n += 6;
        }
    for (i=0; i<new_Nx; i++)
        {
            std::copy ( &coord[6*(Ny*(Nx+1)+ux*i)], &coord[6*(Ny*(Nx+1)+ux*i+1)], &(*new_coord)[n] );
            n += 6;
        }
    std::copy ( &coord[6*(Ny*(Nx+1)+Nx)], &coord[6*(Ny*(Nx+1)+Nx+1)], &(*new_coord)[n] );

    k = 0;
    for (lit=layers.begin();lit!=layers.end();lit++)
        {
            k1 = (*lit);
            lit_next = lit + 1;
            if (lit_next!=layers.end())
                k2 = (*lit_next);
            else
                k2 = Nz;

            // tops
            for (j=0; j<new_Ny-1; j++)
                {
                    for (i=0; i<new_Nx-1; i++)
                        {
                            ind = 2*k1*layer_size + 4*uy*j*Nx;
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                        }
                    ind = 2*k1*layer_size + 4*uy*j*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
                    
                    for (i=0; i<new_Nx-1; i++)
                        {
                            ind = 2*k1*layer_size + (2*uy*(j+1)-1)*2*Nx;
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                        }
                    ind = 2*k1*layer_size + (2*uy*(j+1)-1)*2*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
                }
            
            for (i=0; i<new_Nx-1; i++)
                {
                    ind = 2*k1*layer_size + 4*uy*(new_Ny-1)*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                }
            ind = 2*k1*layer_size + 4*uy*(new_Ny-1)*Nx;
            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
            (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
            
            for (i=0; i<new_Nx-1; i++)
                {
                    ind = 2*k1*layer_size + (2*Ny-1)*2*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                }
            ind = 2*k1*layer_size + (2*Ny-1)*2*Nx;
            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx - 1)];
            (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
            
            // bottoms
            for (j=0; j<new_Ny-1; j++)
                {
                    for (i=0; i<new_Nx-1; i++)
                        {
                            ind = (2*k2-1)*layer_size + 4*uy*j*Nx;
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                        }
                    ind = (2*k2-1)*layer_size + 4*uy*j*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
                    
                    for (i=0; i<new_Nx-1; i++)
                        {
                            ind = (2*k2-1)*layer_size + (2*uy*(j+1)-1)*2*Nx;
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                        }
                    ind = (2*k2-1)*layer_size + (2*uy*(j+1)-1)*2*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
                }
            
            for (i=0; i<new_Nx-1; i++)
                {
                    ind = (2*k2-1)*layer_size + 4*uy*(new_Ny-1)*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                }
            ind = (2*k2-1)*layer_size + 4*uy*(new_Ny-1)*Nx;
            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx-1)];
            (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
            
            for (i=0; i<new_Nx-1; i++)
                {
                    ind = (2*k2-1)*layer_size + (2*Ny-1)*2*Nx;
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*i];
                    (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(i+1) - 1];
                }
            ind = (2*k2-1)*layer_size + (2*Ny-1)*2*Nx;
            (*new_zcorn)[(k++)] = zcorn[ind + 2*ux*(new_Nx - 1)];
            (*new_zcorn)[(k++)] = zcorn[ind + 2*Nx - 1];
        }
   return bp::make_tuple (new_coord, new_zcorn); 
}

t_double upsc::calc_sum_dW ( t_long k1, t_long k2, t_long Nx, t_long Ny, 
                                          spv_float vol_, spv_float ntg_, spv_float poro_, spv_float permx_ )
{
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& permx = *permx_;
    
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
            
            for (n = 0; n < 2; n++)
                {
                    if (poro[index[n]  != 0])
                        p[n] = ntg[index[n]] * permx[index[n]] / poro[index[n]];
                    else
                        p[n] = 0;
                }
            
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

spv_float upsc::upscale_saturation_cube (t_long Nx, t_long Ny, t_long Nz, t_long Nz_upsc,
                                         spv_uint layers_, spv_float vol_, spv_float ntg_, spv_float poro_, spv_float sat_)
{
    t_int i, j, k, n, k1, k2, ind, z;
    t_long layer_size, new_cube_size;
    t_double sum_sat_poro_ntg_vol, sum_poro_ntg_vol;

    v_uint& layers = *layers_;
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& sat = *sat_;
    
    spv_float new_sat = BS_KERNEL.create_object(v_float::bs_type());
    
    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nx*Ny;
    
    layer_size = Nx * Ny;

    for (n = 0; n < Nz_upsc; ++n)
        {
            k1 = layers[n];

            if (n == Nz_upsc-1)
                k2 = Nz;
            else
                k2 = layers[n+1];

            if (k1 != k2)
            for (j = 0; j < Ny; ++j)
                for (i = 0; i < Nx; ++i)
                    {
                        sum_sat_poro_ntg_vol = 0;
                        sum_poro_ntg_vol = 0;

                        for (k = k1; k < k2; ++k)
                            {
                                z = k*layer_size;
                                ind = i + j * Nx + z;
                                sum_sat_poro_ntg_vol += sat[ind]*poro[ind]*ntg[ind]*vol[ind];
                                sum_poro_ntg_vol += poro[ind]*ntg[ind]*vol[ind];
                            }
                        ind = i + j * Nx + k1 * layer_size;
                        if (sum_sat_poro_ntg_vol != 0)
                            sat[ind] = sum_sat_poro_ntg_vol/sum_poro_ntg_vol;
                        else
                            sat[ind] = 0;
                    }

        }
    
    new_cube_size = Nz_upsc*layer_size;
    new_sat->resize(new_cube_size);

    i = 0;
    for ( k = 0; k < Nz_upsc; ++k)
        {
            std::copy ( &sat[k*layer_size], &sat[(k+1)*layer_size], &(*new_sat)[i*layer_size] );
            i++;
        }

    return new_sat;
}

spv_float upsc::upscale_sat_cube_xy (t_long Nx, t_long Ny, t_long ux, t_long uy, t_long Nz, 
                                     spv_uint layers_, spv_float vol_, spv_float ntg_, spv_float poro_, spv_float sat_)
{
    t_int i, j, ii, jj, k, new_k, n, k1, k2, index, ind, new_ind;
    t_long new_Nx, new_Ny, new_Nz;
    t_long layer_size, new_layer_size, new_cube_size;
    t_double sum_sat_poro_ntg_vol, sum_poro_ntg_vol;

    v_uint& layers = *layers_;
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& sat = *sat_;
    
    t_int rx, ry;
    rx = Nx % ux;
    ry = Ny % uy;
    new_Nx = ceil(double(Nx)/double(ux));
    new_Ny = ceil(double(Ny)/double(uy));
    new_Nz = layers.size(); 
    
    spv_float new_sat = BS_KERNEL.create_object(v_float::bs_type());
    new_cube_size = new_Nx*new_Ny*new_Nz;
    new_sat->resize(new_cube_size);
    
    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nx*Ny;
    
    layer_size = Nx * Ny;
    new_layer_size = new_Nx * new_Ny;

    new_k = 0;
    for (n = 0; n < new_Nz; ++n)
      {
        k1 = layers[n];

        if (n == new_Nz-1)
            k2 = Nz;
        else
            k2 = layers[n+1];

        for (j = 0; j < new_Ny-1; ++j)
          {
            for (i = 0; i < new_Nx-1; ++i)
              {
                index  = ux*i + uy*j*Nx + k1*layer_size;
                new_ind = i + j*new_Nx + new_k*new_layer_size;
                
                sum_sat_poro_ntg_vol = 0;
                sum_poro_ntg_vol = 0;
        
                for (k = 0; k < (k2-k1); ++k)
                  {
                    for (jj = 0; jj < uy; jj++)
                      {
                        for (ii = 0; ii < ux; ii++)
                          {
                            ind = index + ii + jj*Nx + k*layer_size;
                            sum_sat_poro_ntg_vol += sat[ind]*poro[ind]*ntg[ind]*vol[ind];
                            sum_poro_ntg_vol += poro[ind]*ntg[ind]*vol[ind];
                          }
                      }
                  }

                if (sum_sat_poro_ntg_vol != 0)
                    (*new_sat)[new_ind] = sum_sat_poro_ntg_vol/sum_poro_ntg_vol;
                else
                    (*new_sat)[new_ind] = 0;
              }

            index  = ux*(new_Nx-1) + uy*j*Nx + k1*layer_size;
            new_ind = new_Nx-1 + j*new_Nx + new_k*new_layer_size;
            
            sum_sat_poro_ntg_vol = 0;
            sum_poro_ntg_vol = 0;
    
            for (k = 0; k < (k2-k1); ++k)
              {
                for (jj = 0; jj < uy; jj++)
                  {
                    for (ii = 0; ii < rx; ii++)
                      {
                        ind = index + ii + jj*Nx + k*layer_size;
                        sum_sat_poro_ntg_vol += sat[ind]*poro[ind]*ntg[ind]*vol[ind];
                        sum_poro_ntg_vol += poro[ind]*ntg[ind]*vol[ind];
                      }
                  }
              }
            
            if (sum_sat_poro_ntg_vol != 0)
                (*new_sat)[new_ind] = sum_sat_poro_ntg_vol/sum_poro_ntg_vol;
            else
                (*new_sat)[new_ind] = 0;
          }
        
        for (i = 0; i < new_Nx-1; ++i)
          {
            index  = ux*i + uy*(new_Ny-1)*Nx + k1*layer_size;
            new_ind = i + (new_Ny-1)*new_Nx + new_k*new_layer_size;
            
            sum_sat_poro_ntg_vol = 0;
            sum_poro_ntg_vol = 0;
    
            for (k = 0; k < (k2-k1); ++k)
              {
                for (jj = 0; jj < ry; jj++)
                  {
                    for (ii = 0; ii < ux; ii++)
                      {
                        ind = index + ii + jj*Nx + k*layer_size;
                        sum_sat_poro_ntg_vol += sat[ind]*poro[ind]*ntg[ind]*vol[ind];
                        sum_poro_ntg_vol += poro[ind]*ntg[ind]*vol[ind];
                      }
                  }
              }
            if (sum_sat_poro_ntg_vol != 0)
                (*new_sat)[new_ind] = sum_sat_poro_ntg_vol/sum_poro_ntg_vol;
            else
                (*new_sat)[new_ind] = 0;
          }

        index  = ux*(new_Nx-1) + uy*(new_Ny-1)*Nx + k1*layer_size;
        new_ind = new_Nx-1 + (new_Ny-1)*new_Nx + new_k*new_layer_size;
        
        sum_sat_poro_ntg_vol = 0;
        sum_poro_ntg_vol = 0;

        for (k = 0; k < (k2-k1); ++k)
          {
            for (jj = 0; jj < ry; jj++)
              {
                for (ii = 0; ii < rx; ii++)
                  {
                    ind = index + ii + jj*Nx + k*layer_size;
                    sum_sat_poro_ntg_vol += sat[ind]*poro[ind]*ntg[ind]*vol[ind];
                    sum_poro_ntg_vol += poro[ind]*ntg[ind]*vol[ind];
                  }
              }
          }  
        if (sum_sat_poro_ntg_vol != 0)
            (*new_sat)[new_ind] = sum_sat_poro_ntg_vol/sum_poro_ntg_vol;
        else
            (*new_sat)[new_ind] = 0;

        new_k ++;
      }

    return new_sat;
}

int upsc::upscale_cubes ( t_long k1, t_long k2, t_long Nx, t_long Ny, 
                          spv_float vol_, spv_float ntg_, spv_float poro_, spv_float permx_ )
{
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& permx = *permx_;
    
    t_long i, j, index[2], z1, z2, layer_size;
    t_double vol_sum, ntg_vol_sum, poro_ntg_vol_sum, permx_ntg_vol_sum;

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
            permx_ntg_vol_sum = permx[index[0]]*ntg[index[0]]*vol[index[0]] + permx[index[1]]*ntg[index[1]]*vol[index[1]];
            
            vol[index[0]] = vol_sum;
            if (vol_sum != 0)
                ntg[index[0]] = ntg_vol_sum/vol_sum;
            else
                ntg[index[0]] = 0;
            
            if (ntg_vol_sum != 0)
                {
                    poro[index[0]] = poro_ntg_vol_sum/ntg_vol_sum;
                    permx[index[0]] = permx_ntg_vol_sum/ntg_vol_sum;
                }
            else
                {
                    poro[index[0]] = 0;
                    permx[index[0]] = 0;
                }
            
          }
      }

    return 0;
}

bp::tuple upsc::upscale_cubes_xy ( t_long Nx, t_long Ny, t_long Nz, t_long ux, t_long uy,
                                   spv_float vol_, spv_float ntg_, spv_float poro_)
{
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    
    spv_float new_vol = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_ntg = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_poro = BS_KERNEL.create_object(v_float::bs_type());

    t_long new_cube_size, new_Nx, new_Ny;
    t_int rx, ry;
    new_Nx = ceil(double(Nx)/double(ux));
    new_Ny = ceil(double(Ny)/double(uy));
    rx = Nx % ux;
    ry = Ny % uy;
    new_cube_size = new_Nx*new_Ny*Nz;
    
    new_vol->resize(new_cube_size);
    new_ntg->resize(new_cube_size);
    new_poro->resize(new_cube_size);
    
    t_long i, j, ii, jj, k, layer_size, new_layer_size;
    t_long ind, index, new_ind;
    t_double vol_sum, ntg_vol_sum, poro_ntg_vol_sum;
    layer_size = Nx*Ny;
    new_layer_size = new_Nx*new_Ny;

    for (k = 0; k < Nz; k++)
      {
        for (j = 0; j < new_Ny-1; j++)
          {
            for (i = 0; i < new_Nx-1; i++)
              {
                ind  = ux*i + uy*j*Nx + k*layer_size;
                new_ind = i + j*new_Nx + k*new_layer_size;

                vol_sum = 0;
                ntg_vol_sum = 0;
                poro_ntg_vol_sum = 0; 

                for (jj = 0; jj < uy; jj++)
                  {
                    for (ii = 0; ii < ux; ii++)
                      {
                        index = ind + ii;
                        vol_sum += vol[index];
                        ntg_vol_sum += ntg[index]*vol[index];
                        poro_ntg_vol_sum += poro[index]*ntg[index]*vol[index];
                      }
                    ind  += Nx;
                  }
                (*new_vol)[new_ind] = vol_sum;
                
                if (vol_sum != 0)
                    (*new_ntg)[new_ind] = ntg_vol_sum/vol_sum;
                else
                    (*new_ntg)[new_ind] = 0;
                
                if (ntg_vol_sum != 0)
                  {
                    (*new_poro)[new_ind] = poro_ntg_vol_sum/ntg_vol_sum;
                  }
                else
                  {
                    (*new_poro)[new_ind] = 0;
                  }
              }
            
            ind  = ux*(new_Nx-1) + uy*j*Nx + k*layer_size;
            new_ind = new_Nx-1 + j*new_Nx + k*new_layer_size;

            vol_sum = 0;
            ntg_vol_sum = 0;
            poro_ntg_vol_sum = 0; 

            for (jj = 0; jj < uy; jj++)
              {
                for (ii = 0; ii < rx; ii++)
                  {
                    index = ind + ii;
                    vol_sum += vol[index];
                    ntg_vol_sum += ntg[index]*vol[index];
                    poro_ntg_vol_sum += poro[index]*ntg[index]*vol[index];
                  }
                ind  += Nx;
              }
            (*new_vol)[new_ind] = vol_sum;
            
            if (vol_sum != 0)
                (*new_ntg)[new_ind] = ntg_vol_sum/vol_sum;
            else
                (*new_ntg)[new_ind] = 0;
            
            if (ntg_vol_sum != 0)
              {
                (*new_poro)[new_ind] = poro_ntg_vol_sum/ntg_vol_sum;
              }
            else
              {
                (*new_poro)[new_ind] = 0;
              }
          }
        
        for (i = 0; i < new_Nx-1; i++)
          {
            ind  = ux*i + uy*(new_Ny-1)*Nx + k*layer_size;
            new_ind = i + (new_Ny-1)*new_Nx + k*new_layer_size;

            vol_sum = 0;
            ntg_vol_sum = 0;
            poro_ntg_vol_sum = 0; 

            for (jj = 0; jj < ry; jj++)
              {
                for (ii = 0; ii < ux; ii++)
                  {
                    index = ind + ii;
                    vol_sum += vol[index];
                    ntg_vol_sum += ntg[index]*vol[index];
                    poro_ntg_vol_sum += poro[index]*ntg[index]*vol[index];
                  }
                ind  += Nx;
              }
            (*new_vol)[new_ind] = vol_sum;
            
            if (vol_sum != 0)
                (*new_ntg)[new_ind] = ntg_vol_sum/vol_sum;
            else
                (*new_ntg)[new_ind] = 0;
            
            if (ntg_vol_sum != 0)
              {
                (*new_poro)[new_ind] = poro_ntg_vol_sum/ntg_vol_sum;
              }
            else
              {
                (*new_poro)[new_ind] = 0;
              }
          }
        
        ind  = ux*(new_Nx-1) + uy*(new_Ny-1)*Nx + k*layer_size;
        new_ind = new_Nx-1 + (new_Ny-1)*new_Nx + k*new_layer_size;

        vol_sum = 0;
        ntg_vol_sum = 0;
        poro_ntg_vol_sum = 0; 

        for (jj = 0; jj < ry; jj++)
          {
            for (ii = 0; ii < rx; ii++)
              {
                index = ind + ii;
                vol_sum += vol[index];
                ntg_vol_sum += ntg[index]*vol[index];
                poro_ntg_vol_sum += poro[index]*ntg[index]*vol[index];
              }
            ind  += Nx;
          }
        (*new_vol)[new_ind] = vol_sum;
        
        if (vol_sum != 0)
            (*new_ntg)[new_ind] = ntg_vol_sum/vol_sum;
        else
            (*new_ntg)[new_ind] = 0;
        
        if (ntg_vol_sum != 0)
          {
            (*new_poro)[new_ind] = poro_ntg_vol_sum/ntg_vol_sum;
          }
        else
          {
            (*new_poro)[new_ind] = 0;
          }
      }

    return bp::make_tuple(new_ntg, new_poro);
}

t_double upsc::solve_pressure_zcolumn (t_long Ny, t_long Nz, t_long i, t_long j, t_long k1, t_long k2, BS_SP(rs_mesh_iface) sp_mesh_iface)
{
    t_long size, k, n, index, ext_ind[2];
    plane_t plane[2];
    mesh_element3d element[2];
    fpoint3d_t center[2];

    t_double tz;
    t_double p_left = 1, p_right = 0;
    t_double dl, dL, dp, dP, K;

    spv_double p = BS_KERNEL.create_object(v_double::bs_type());
    spv_double rhs = BS_KERNEL.create_object(v_double::bs_type());
    spv_float tran_vals = BS_KERNEL.create_object(v_float::bs_type());
    sp_dens_mtx_t tran = BS_KERNEL.create_object("dens_matrix");
    sp_blu_solver_t solver = BS_KERNEL.create_object("blu_solver");

    size = k2 - k1 +1;
    
    p->resize(size);
    rhs->resize(size);

    tran->init(size, size, 60);
    tran_vals = tran->get_values();
    
    v_float& A = *tran_vals;  
    v_double& b = *rhs;

    smart_ptr<bs_mesh_grdecl> sp_mesh(sp_mesh_iface, bs_static_cast());
    mesh_grdecl mesh = sp_mesh->get_wrapped();
    
    // important: ZYX order
    // direction dir along_dim1->X (along_dim2->Y, along_dim3->Z) 

    // k = k1 //////////////////////////////////////////////////////////////////////////////////////////
    ext_ind[0] = k1 + j * Nz + i * Ny * Nz;
    mesh.calc_element (i, j, k1, element[0]);
    center[0] = element[0].get_center();
    element[0].get_plane (z_axis_minus, plane[0]);
    tz = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim3);

    // FIXME
    b[0] = -tz*p_left;
    A[0] -= tz;          
    
    // k1 < k < k2 /////////////////////////////////////////////////////////////////////////////////////
    for (k = k1; k < k2; ++k)
        {
            ext_ind[0] = k + j * Nz + i * Ny * Nz;
            mesh.calc_element (i, j, k, element[0]);
            center[0] = element[0].get_center();
            
            ext_ind[1] = ext_ind[0] + 1;
            mesh.calc_element (i, j, k + 1, element[1]);
            element[0].get_plane (z_axis_plus, plane[0]);
            
            //element[1].get_plane (z_axis_minus, plane[1]);
            
            center[1] = element[1].get_center();
            // don't pass plane[1] to calc_tran
            // because zcolumn cells are always fully adjacent
            //tz = mesh.calc_tran(ext_ind[0], ext_ind[1], plane[0], center[0], center[1], along_dim3, &plane[1]);
            tz = mesh.calc_tran(ext_ind[0], ext_ind[1], plane[0], center[0], center[1], along_dim3);
            
            n = k - k1; 
            // Z-
            index = n + (n+1)*size;
            A[index] = tz;
            index += 1;
            A[index] -= tz;
            
            // Z+
            index = (n+1) + n*size;
            A[index] = tz;
            index -= 1;
            A[index] -= tz;
        }
    
    // k = k2 //////////////////////////////////////////////////////////////////////////////////////////
    ext_ind[0] = k2 + j * Nz + i * Ny * Nz;
    mesh.calc_element (i, j, k2, element[0]);
    center[0] = element[0].get_center();
    element[0].get_plane (z_axis_plus, plane[0]);
    tz = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim3);
    
    b[size-1] = -tz*p_right;
    index  =  size*size - 1;
    A[index] -= tz;          

    // solve system, find pressure
    solver->setup (tran);
    solver->solve (tran, rhs, p);

    dp = (p_left - (*p)[0]);
    dP = ( (p_left + (*p)[0])/2 - ((*p)[size-1] + p_right)/2 );
    
    mesh.calc_element (i, j, k1, element[0]);
    element[0].get_plane (z_axis_minus, plane[0]);
    center[1] = element[0].get_center ();
    get_plane_center (plane[0],  center[0]);
    
    dl = get_len(center[0], center[1])*2;
        
    mesh.calc_element (i, j, k2, element[1]);
    element[1].get_plane (z_axis_plus, plane[1]);
    get_plane_center (plane[1],  center[1]);
    
    dL = get_len(center[0], center[1]);
    
    if (dl < 10*EPSILON)
        {
            dl = 2*EPSILON;
            dL += EPSILON;
        }

    K = (dp/dl)/(dP/dL);

    return K;
}

t_double upsc::solve_pressure_block (t_int direction, t_long Ny, t_long Nz, t_long i1, t_long i2, t_long j1, t_long j2, t_long k1, t_long k2, spv_float ntg_, BS_SP(rs_mesh_iface) sp_mesh_iface)
{
    const t_double cdarcy = 0.008640;
    t_long size, i, j, k, ind[2], ext_ind[2], index, ny, nz;
    plane_t plane[2];
    mesh_element3d element[2];
    fpoint3d_t center[2];

    t_double tz, ty, tx;
    t_double p_left = 1, p_right = 0;
    t_double dL, dP, S, K, Q;

    spv_double p = BS_KERNEL.create_object(v_double::bs_type());
    spv_double rhs = BS_KERNEL.create_object(v_double::bs_type());
    spv_float tran_vals = BS_KERNEL.create_object(v_float::bs_type());
    sp_dens_mtx_t tran = BS_KERNEL.create_object("dens_matrix");
    sp_blu_solver_t solver = BS_KERNEL.create_object("blu_solver");

    size = (i2-i1)*(j2-j1)*(k2-k1);
    ny = j2 - j1;
    nz = k2 - k1;
    
    p->resize(size);
    rhs->resize(size);

    tran->init(size, size, 60);
    tran_vals = tran->get_values();
    
    v_float& A = *tran_vals;
    v_double& b = *rhs;
    v_float& ntg = *ntg_;

    smart_ptr<bs_mesh_grdecl> sp_mesh(sp_mesh_iface, bs_static_cast());
    mesh_grdecl mesh = sp_mesh->get_wrapped();
    
    // important: ZYX order

    // direction dir along_dim1->X (2->Y, 3->Z) 
 
    for (i = i1; i < i2; ++i)
        for (j = j1; j < j2; ++j)
            for (k = k1; k < k2; ++k)
              {
                ext_ind[0] = k + j * Nz + i * Ny * Nz;
                ind[0] = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
                mesh.calc_element (i, j, k, element[0]);
                center[0] = element[0].get_center();
                
                // TRANZ /////////////////////////////////////////////////////////////////////////////////////////////
                if (k!=k2-1)
                    {
                        ext_ind[1] = ext_ind[0] + 1;
                        ind[1] = ind[0] + 1;
                        mesh.calc_element (i, j, k + 1, element[1]);
                        element[0].get_plane (z_axis_plus, plane[0]);
                        element[1].get_plane (z_axis_minus, plane[1]);
                        center[1] = element[1].get_center();
                        tz = mesh.calc_tran(ext_ind[0], ext_ind[1], plane[0], center[0], center[1], along_dim3, &plane[1]);
                        
                        // Z-
                        index = ind[0] + ind[1]*size;
                        A[index] = tz;
                        index += 1;
                        A[index] -= tz;
                        
                        // Z+
                        index = ind[1] + ind[0]*size;
                        A[index] = tz;
                        index -= 1;
                        A[index] -= tz;
                    }
                
                // TRANY /////////////////////////////////////////////////////////////////////////////////////////////
                if (j != j2-1)
                    {
                        ext_ind[1] = ext_ind[0] + Nz;
                        ind[1] = ind[0] + nz;
                        mesh.calc_element (i, j + 1, k, element[1]);
                        element[0].get_plane (y_axis_plus, plane[0]);
                        element[1].get_plane (y_axis_minus, plane[1]);
                        center[1] = element[1].get_center();
                        ty = mesh.calc_tran(ext_ind[0], ext_ind[1], plane[0], center[0], center[1], along_dim2, &plane[1]);          
                        
                        // Y- 
                        index = ind[0] + ind[1]*size;
                        A[index] = ty;
                        index += nz;
                        A[index] -= ty;
                        
                        // Y+
                        index = ind[1] + ind[0]*size;
                        A[index] = ty;
                        index -= nz;
                        A[index] -= ty;
                    }

                // TRANX /////////////////////////////////////////////////////////////////////////////////////////////
                if (i != i2-1)
                    {
                        ext_ind[1] = ext_ind[0] + Ny * Nz;
                        ind[1] = ind[0] + ny*nz;
                        mesh.calc_element (i + 1, j, k, element[1]);
                        element[0].get_plane (x_axis_plus, plane[0]);
                        element[1].get_plane (x_axis_minus, plane[1]);
                        center[1] = element[1].get_center();
                        tx = mesh.calc_tran(ext_ind[0], ext_ind[1], plane[0], center[0], center[1], along_dim1, &plane[1]);
                        
                        // X-
                        index = ind[0] + ind[1]*size;
                        A[index] = tx;
                        index += ny*nz;
                        A[index] -= tx;
                                                
                        // X+
                        index = ind[1] + ind[0]*size;
                        A[index] = tx;
                        index -= ny*nz;
                        A[index] -= tx;

                    }
               // Boundary Conditions /////////////////////////////////////////////////////////////////////////////////
               if (direction == 1)
                 {
                    if (i == i1)
                        {
                           element[0].get_plane (x_axis_minus, plane[0]);
                           tx = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim1);
                           b[ind[0]] = -tx*p_left;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= tx;          
                        }
                    if (i == i2-1)
                        {  
                           element[0].get_plane (x_axis_plus, plane[0]);
                           tx = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim1);
                           b[ind[0]] = -tx*p_right;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= tx;          
                        }
                 }
               else if (direction == 2)
                 {
                    if (j == j1)
                        {
                           element[0].get_plane (y_axis_minus, plane[0]);
                           ty = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim2);
                           b[ind[0]] = -ty*p_left;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= ty;          
                        }
                    if (j == j2-1)
                        {  
                           element[0].get_plane (y_axis_plus, plane[0]);
                           ty = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim2);
                           b[ind[0]] = -ty*p_right;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= ty;          
                        }
                 }
               else if (direction == 3)
                 {
                    if (k == k1)
                        {
                           element[0].get_plane (z_axis_minus, plane[0]);
                           tz = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim3);
                           b[ind[0]] = -tz*p_left;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= tz;          
                        }
                    if (k == k2-1)
                        {  
                           element[0].get_plane (z_axis_plus, plane[0]);
                           tz = mesh.calc_tran_boundary (ext_ind[0], plane[0], center[0], along_dim3);
                           b[ind[0]] = -tz*p_right;
                           index  =  ind[0] + ind[0]*size;
                           A[index] -= tz;          
                        }
                 }

              }

    solver->setup (tran);
    solver->solve (tran, rhs, p);

    // Q=sum(q) q=T*dp
    // Q=k*S*dP/dL
    // for permx and permy don't forget NTG! 
    
    Q = 0;
    S = 0;
    
    if (direction == 1)
      {
        i = i1;
        for (j = j1; j < j2; ++j)
            for (k = k1; k < k2; ++k)
              {
                mesh.calc_element (i, j, k, element[0]);
                element[0].get_plane (x_axis_minus, plane[0]);
                center[0] = element[0].get_center ();
                get_plane_center (plane[0],  center[1]);
                
                index = k + j * Nz + i * Ny * Nz;
                S += ntg[index]*find_area_of_side (plane[0][0], plane[0][1], plane[0][2], plane[0][3]);
                tx = mesh.calc_tran_boundary (index, plane[0], center[0], along_dim1);

                index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
                dP = (p_left - (*p)[index]);
                Q += tx*dP;
              }
      }
    else if (direction == 2)
      {
        j = j1;
        for (i = i1; i < i2; ++i)
            for (k = k1; k < k2; ++k)
              {
                mesh.calc_element (i, j, k, element[0]);
                element[0].get_plane (y_axis_minus, plane[0]);
                center[0] = element[0].get_center ();
                get_plane_center (plane[0],  center[1]);
                
                index = k + j * Nz + i * Ny * Nz;
                S += ntg[index]*find_area_of_side (plane[0][0], plane[0][1], plane[0][2], plane[0][3]);
                ty = mesh.calc_tran_boundary (index, plane[0], center[0], along_dim2);
                
                index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
                dP = (p_left - (*p)[index]);
                Q += ty*dP;
              }
      }
    else if (direction == 3)
      {
        k = k1;
        for (i = i1; i < i2; ++i)
            for (j = j1; j < j2; ++j)
              {
                mesh.calc_element (i, j, k, element[0]);
                element[0].get_plane (z_axis_minus, plane[0]);
                center[0] = element[0].get_center ();
                get_plane_center (plane[0],  center[1]);
                
                index = k + j * Nz + i * Ny * Nz;
                S += find_area_of_side (plane[0][0], plane[0][1], plane[0][2], plane[0][3]);
                tz = mesh.calc_tran_boundary (index, plane[0], center[0], along_dim3);
                
                index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
                dP = (p_left - (*p)[index]);
                Q += tz*dP;
              }
      }

    i = i1; j = j1; k = k1;
    if (direction == 1)
      {
        mesh.calc_element (i, j, k, element[0]);
        element[0].get_plane (x_axis_minus, plane[0]);
        get_plane_center (plane[0],  center[0]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP = (p_left + (*p)[index])/2;
        
        i = i2-1;
        mesh.calc_element (i, j, k, element[1]);
        element[1].get_plane (x_axis_plus, plane[1]);
        get_plane_center (plane[1],  center[1]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP -= ((*p)[index] + p_right)/2;
      }
    else if (direction == 2)
      {
        mesh.calc_element (i, j, k, element[0]);
        element[0].get_plane (y_axis_minus, plane[0]);
        get_plane_center (plane[0],  center[0]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP = (p_left + (*p)[index])/2;
        
        j = j2-1; 
        mesh.calc_element (i, j, k, element[1]);
        element[1].get_plane (y_axis_plus, plane[1]);
        get_plane_center (plane[1],  center[1]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP -= ((*p)[index] + p_right)/2;
      }
    else if (direction == 3)
      {
        mesh.calc_element (i, j, k, element[0]);
        element[0].get_plane (z_axis_minus, plane[0]);
        get_plane_center (plane[0],  center[0]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP = (p_left + (*p)[index])/2;
        
        k = k2-1;
        mesh.calc_element (i, j, k, element[1]);
        element[1].get_plane (z_axis_plus, plane[1]);
        get_plane_center (plane[1],  center[1]);
        
        index = (k-k1) + (j-j1)*nz + (i-i1)*ny*nz;
        dP -= ((*p)[index] + p_right)/2;
      }

    dL = get_len(center[0], center[1]);
    if (S < EPSILON || dP < EPSILON)
        K = 0;
    else
        K = Q*dL/(cdarcy*S*dP);
    
    return K;

}

spv_float upsc::upscale_permz_zcolumn (t_long Nx, t_long Ny, t_long Nz, spv_uint layers_, spv_float permz_, spv_uint actnum_, BS_SP(rs_mesh_iface) sp_mesh_iface)
  {
    t_int i, j, k, n, Nz_upsc, k1, k2, z1, ind, index;
    t_long layer_size, new_cube_size;

    t_double upsc_factor;
    
    spv_float new_permz = BS_KERNEL.create_object(v_float::bs_type());
    v_float& permz = *permz_;
    v_uint& actnum = *actnum_;
    v_uint& layers = *layers_;
    
    smart_ptr<bs_mesh_grdecl> sp_mesh(sp_mesh_iface, bs_static_cast());
    mesh_grdecl mesh = sp_mesh->get_wrapped();
    
    Nz_upsc = layers.size();
    layer_size = Nx*Ny;
    
    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nx*Ny;
    
    for (n = 0; n < Nz_upsc; ++n)
        {
            k1 = layers[n];

            if (n == Nz_upsc-1)
                k2 = Nz - 1;
            else
                k2 = layers[n+1] - 1;

            if (k1 != k2)
                {
                    z1 = k1*layer_size;
                    for (j = 0; j < Ny; ++j)
                        for (i = 0; i < Nx; ++i)
                            {
                                index = i + j * Nx + z1;

                                // in general case call function that finds isolated bodies
                                for (k = k1; k <= k2; k++)
                                    {
                                        ind = i + j * Nx + k * Nx * Ny;
                                        if (!actnum [ind])
                                            {
                                                permz[index] = 0;
                                                break;
                                            }
                                    }
                                if (permz[index])
                                    {
                                        upsc_factor = solve_pressure_zcolumn (Ny, Nz, i, j, k1, k2, sp_mesh_iface);
                                        permz[index] *= upsc_factor;
                                    }

                            }
                }   
        }  

    new_cube_size = Nz_upsc*layer_size;
    new_permz->resize(new_cube_size);
    
    i = 0;
    for ( n = 0; n < Nz_upsc; n++)
        {
            k = layers[n];
            std::copy ( &permz[k*layer_size], &permz[(k+1)*layer_size], &(*new_permz)[i*layer_size] );
            i++;
        }

    return new_permz;
  }

spv_float upsc::upscale_perm_block (t_int dir, t_long Nx, t_long Ny, t_long ux, t_long uy, t_long Nz, spv_uint layers_, spv_float ntg_, spv_uint actnum_, BS_SP(rs_mesh_iface) sp_mesh_iface)
  {
    t_int i, j, n, i1, i2, j1, j2, k1, k2, new_k, new_ind;
    t_long new_Nx, new_Ny, new_Nz;
    t_long new_layer_size, new_cube_size;

    spv_float new_perm = BS_KERNEL.create_object(v_float::bs_type());
    v_uint& actnum = *actnum_;
    v_uint& layers = *layers_;
    
    smart_ptr<bs_mesh_grdecl> sp_mesh(sp_mesh_iface, bs_static_cast());
    mesh_grdecl mesh = sp_mesh->get_wrapped();
    
    t_int rx, ry;
    rx = Nx % ux;
    ry = Ny % uy;
    new_Nx = ceil(double(Nx)/double(ux));
    new_Ny = ceil(double(Ny)/double(uy));
    new_Nz = layers.size(); 
    new_layer_size = new_Nx * new_Ny;
    new_cube_size = new_Nx*new_Ny*new_Nz;
    new_perm->resize(new_cube_size);
    
    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nx*Ny;
    
    // FIXME: call function that finds isolated bodies

    new_k = 0;
    for (n = 0; n < new_Nz; ++n)
      {
        k1 = layers[n];

        if (n == new_Nz-1)
            k2 = Nz;
        else
            k2 = layers[n+1] - 1;

        for (j = 0; j < new_Ny-1; ++j)
          {
            for (i = 0; i < new_Nx-1; ++i)
              {
                i1 = ux*i;
                i2 = i1 + ux;
                j1 = uy*j;
                j2 = j1 + uy;
                new_ind = i + j*new_Nx + new_k*new_layer_size;
                (*new_perm)[new_ind] = solve_pressure_block (dir, Ny, Nz, i1, i2, j1, j2, k1, k2, ntg_, sp_mesh_iface);
              }
            
            if (rx)
              {
                i1 = ux*(new_Nx-1);
                i2 = i1 + rx;
                j1 = uy*j;
                j2 = j1 + uy;
                new_ind = new_Nx-1 + j*new_Nx + new_k*new_layer_size;
                (*new_perm)[new_ind] = solve_pressure_block (dir, Ny, Nz, i1, i2, j1, j2, k1, k2, ntg_, sp_mesh_iface);
              }
          }
           
        if (ry)
          {
            for (i = 0; i < new_Nx-1; ++i)
              {
                i1 = ux*i;
                i2 = i1 + ux;
                j1 = uy*(new_Ny-1);
                j2 = j1 + ry;
                new_ind = i + (new_Ny-1)*new_Nx + new_k*new_layer_size;
                (*new_perm)[new_ind] = solve_pressure_block (dir, Ny, Nz, i1, i2, j1, j2, k1, k2, ntg_, sp_mesh_iface);

              }

            if (rx)
              {
                i1 = ux*(new_Nx-1);
                i2 = i1 + rx;
                j1 = uy*(new_Ny-1);
                j2 = j1 + ry;
                new_ind = new_Nx-1 + (new_Ny-1)*new_Nx + new_k*new_layer_size;
                (*new_perm)[new_ind] = solve_pressure_block (dir, Ny, Nz, i1, i2, j1, j2, k1, k2, ntg_, sp_mesh_iface);
              }
          }

        new_k ++;
      }
    
    return new_perm;
  }


  bp::tuple
  upsc::king_method (t_long Nx, t_long Ny, t_long Nz, t_long Nz_upsc,
                            spv_float vol_, spv_float ntg_, spv_float poro_, spv_float permx_)
  {
    t_int i, k, n, k0, k1, k2, k3;
    t_long layer_size, new_cube_size;
    t_float sum_dW;
    
    v_float& vol = *vol_;
    v_float& ntg = *ntg_;
    v_float& poro = *poro_;
    v_float& permx = *permx_;

    std::list <t_int> layers;
    std::list <t_int>::iterator lit, tmp_it;

    layer_mmap_t sum_dW2layer;
    layer_mmap_it_t mit;
    layer_map_t layer2mmap_it;
    layer_map_it_t lmit;
    
    spv_float new_vol = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_ntg = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_poro = BS_KERNEL.create_object(v_float::bs_type());
    spv_float new_permx = BS_KERNEL.create_object(v_float::bs_type());
    spv_uint layers_v = BS_KERNEL.create_object(v_uint::bs_type());

    layer_size = Nx*Ny;

    for (i = 0;i < Nz; i++)
        layers.push_back(i);

    // important: XYZ order
    // index <-- (i, j, k)
    // index = i + j*Nx + k*Nx*Ny;

    // calculate sum_dW for every pair of adjacent layers
    // and store it in sumdw2layer multimap (sum_dW -> k) sorted in increasing order
    // and in layer2mmapit map (k -> mmap_it) 
    for (k = 0; k < Nz-1; k++)
    {
        sum_dW = calc_sum_dW ( k, k+1,  Nx, Ny, vol_, ntg_, poro_, permx_ );
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
        upscale_cubes (k1, k2, Nx, Ny, vol_, ntg_, poro_, permx_ );
        
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
                sum_dW = calc_sum_dW ( k1, k3, Nx, Ny, vol_, ntg_, poro_, permx_ );
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
                sum_dW = calc_sum_dW ( k0, k1, Nx, Ny, vol_, ntg_, poro_, permx_ );
                mit = sum_dW2layer.insert ( std::make_pair(sum_dW, k0) );
                layer2mmap_it[k0] = mit;
            }

    }
    
    new_cube_size = Nz_upsc*layer_size;

    new_vol->resize(new_cube_size);
    new_ntg->resize(new_cube_size);
    new_poro->resize(new_cube_size);
    new_permx->resize(new_cube_size);
    layers_v->resize(Nz_upsc);

    i = 0;
    for (lit=layers.begin();lit!=layers.end();lit++)
        {
            k = (*lit);
            std::copy ( &vol[k*layer_size], &vol[(k+1)*layer_size], &(*new_vol)[i*layer_size] );
            std::copy ( &ntg[k*layer_size], &ntg[(k+1)*layer_size], &(*new_ntg)[i*layer_size] );
            std::copy ( &poro[k*layer_size], &poro[(k+1)*layer_size], &(*new_poro)[i*layer_size] );
            std::copy ( &permx[k*layer_size], &permx[(k+1)*layer_size], &(*new_permx)[i*layer_size] );
            i++;
        }

    std::copy ( layers.begin(), layers.end(), &(*layers_v)[0] );
    
    return bp::make_tuple (layers_v, new_vol, new_ntg, new_poro, new_permx) ;
  }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  upsc::py_str () const
    {
      std::stringstream s;
      s << wstr2str (sp_prop->py_str ()) << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (upsc);
  BLUE_SKY_TYPE_STD_COPY (upsc);

  BLUE_SKY_TYPE_IMPL(upsc,  upsc_iface, "upsc", "upscaling", "hdm upscaling");

}  // blue_sky namespace
