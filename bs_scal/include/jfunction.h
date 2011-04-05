/**
 * \file jfunction.h
 * \brief jfunction (for calculate capillary)
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#ifndef BS_JFUNCTION_H_
#define BS_JFUNCTION_H_

#include "plane_orientation.h"

namespace blue_sky
  {

  enum JFUNC_PERM_TYPE_ENUM
  {
    JFUNC_PERM_XY = 0,
    JFUNC_PERM_X,
    JFUNC_PERM_Y,
    JFUNC_PERM_Z,
  };

  class BS_API_PLUGIN jfunction : public objbase
    {
    public:
      //jfunction (JFUNC_PERM_TYPE_ENUM perm_type, double st_phases, double alpha, double beta, bool is_valid = false)
      //  : perm_type (perm_type)
      //  , is_valid (is_valid)
      //{
      //  init (perm_type);
      //}

      typedef t_double                    item_t;
      typedef jfunction                   this_t;

      void init (JFUNC_PERM_TYPE_ENUM perm_type_)
      {
        perm_type = perm_type_;
        is_valid = true;

        if (perm_type == JFUNC_PERM_XY)
          plane_a = PLANE_ORIENTATION_X, plane_b = PLANE_ORIENTATION_Y;
        else if (perm_type == JFUNC_PERM_X)
          plane_a = plane_b = PLANE_ORIENTATION_X;
        else if (perm_type == JFUNC_PERM_Y)
          plane_a = plane_b = PLANE_ORIENTATION_Y;
        else if (perm_type == JFUNC_PERM_Z)
          plane_a = plane_b = PLANE_ORIENTATION_Z;
        else
          {
            is_valid = false;
            throw bs_exception ("jfunction ", "invalid jfunction type");
          }
      }

      bool valid () const
        {
          return is_valid;
        }

      void
      set_valid (bool b)
      {
        is_valid = b;
      }

      item_t get_perm (const item_t *perm_array) const
        {
          BS_ASSERT (is_valid);
          return (perm_array[plane_a] + perm_array[plane_b]) * 0.5;
        }

      BLUE_SKY_TYPE_DECL_T (jfunction);

public:

      item_t    st_phase;
      item_t    alpha;
      item_t    beta;

private:

      bool      is_valid;
      int       plane_a;
      int       plane_b;

      JFUNC_PERM_TYPE_ENUM  perm_type;

    };


} // namespace blue_sky

#endif  // #ifndef BS_JFUNCTION_H_
