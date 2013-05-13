/// @file well_completion_coords.h
/// @brief Completion coords structure definition
/// @author uentity
/// @version 1.0
/// @date 13.05.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WELL_COMPLETION_COORDS_2B3UD1U7
#define WELL_COMPLETION_COORDS_2B3UD1U7

#include "rs_smesh_iface.h"

namespace blue_sky
{
  namespace wells
  {
    struct completion_coords
    {
      auto_value <bool, false> use_CCF_flag; 
      auto_value <t_float, 0> completion_length;  //!< length of current complition
      t_float x1[3];   //!< start point of completion 
      t_float x2[3];   //!< end point of completion 
      t_float compute_length ()
        {
          return (completion_length = sqrt ((x2[0] - x1[0]) * (x2[0] - x1[0]) + (x2[1] - x1[1]) * (x2[1] - x1[1]) + (x2[2] - x1[2]) * (x2[2] - x1[2])));
        }
    };

}} // eof blue_sky::wells

#endif /* end of include guard: WELL_COMPLETION_COORDS_2B3UD1U7 */

