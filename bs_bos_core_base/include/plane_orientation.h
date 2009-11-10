/**
 * \file plane_orientation.h
 * \brief plane orientation
 * \author Sergey Miryanov
 * \date 22.05.2008
 * */
#ifndef BS_PLANE_ORIENTATION_H_
#define BS_PLANE_ORIENTATION_H_


namespace blue_sky
  {

  //! Plane orientations
  enum
  {
    PLANE_ORIENTATION_X = 0,
    PLANE_ORIENTATION_Y,
    PLANE_ORIENTATION_Z,

    PLANE_ORIENTATION_TOTAL
  };

#define FI_PLANE_ORI_X_IND(I)  ((I) * PLANE_ORIENTATION_TOTAL + PLANE_ORIENTATION_X)
#define FI_PLANE_ORI_Y_IND(I)  ((I) * PLANE_ORIENTATION_TOTAL + PLANE_ORIENTATION_Y)
#define FI_PLANE_ORI_Z_IND(I)  ((I) * PLANE_ORIENTATION_TOTAL + PLANE_ORIENTATION_Z)

} // namespace blue_sky


#endif  // #ifndef BS_PLANE_ORIENTATION_H_
