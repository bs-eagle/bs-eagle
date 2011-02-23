/** 
 * @file dens_matrix_tools.cpp
 * @brief 
 * @date 2009-11-24
 */
#include <time.h>
#include <stdlib.h>
#include <map>

#include "bs_mtx_stdafx.h"
#include "dens_matrix_tools.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  
  dens_matrix_tools::dens_matrix_tools (bs_type_ctor_param) 
        : dens_matrix_tools_iface ()
    {
    }
  
  dens_matrix_tools ::dens_matrix_tools (const dens_matrix_tools & /*src*/) : bs_refcounter ()
     {
     }

 int
dens_matrix_tools::random_init (sp_dens_t matrix, 
                                         const t_long ni, 
                                         const t_long nj,
                                         const t_long block_size,
                                         const t_double rand_value_dispersion 
                                         ) const
{
  if (matrix->init (ni, nj, block_size))
    return -2;

  t_float *v = &(*(matrix->get_values ()))[0];

  srand ((unsigned)time( NULL ));
  
  for (t_long i = 0; i < ni; ++i)
    {
      for (t_long j = 0; j < nj; ++j)
        {
          v[i * nj + j] = (t_float)rand () / (t_float)RAND_MAX * (t_float)rand_value_dispersion;
        }
    }
  return 0;
}
      
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (dens_matrix_tools);
  BLUE_SKY_TYPE_STD_COPY (dens_matrix_tools);

  BLUE_SKY_TYPE_IMPL (dens_matrix_tools, dens_matrix_tools_iface, "dens_matrix_tools", "Tools for Dens Matrix class", "Tools realization of Dens Matricies");
}  // blue_sky namespace
