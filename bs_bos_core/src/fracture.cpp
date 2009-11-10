#include "stdafx.h"
#include "fracture.h"

namespace blue_sky
  {
  fracture::fracture()
  {
    init();
  }

  fracture::fracture(const fracture &src)
  {
    *this = src;
  }

  fracture::~fracture()
  {
    init();
  }

  fracture &fracture::operator=(const fracture &src)
  {
    name = src.name;

    dim[0] = src.dim[0];
    dim[1] = src.dim[0];
    dim[2] = src.dim[0];
    dim[3] = src.dim[0];
    frac_flag = src.frac_flag;
    frac_len = src.frac_len;
    frac_angle = src.frac_angle;
    frac_skin = src.frac_skin;
    frac_perm = src.frac_perm;
    frac_width = src.frac_width;

    connection_type = src.connection_type;
    start_time = src.start_time;
    is_activated = src.is_activated;
    ignore_fracture = src.ignore_fracture;

    return *this;
  }

  void fracture::init()
  {
    name = "";

    dim[0] = 0;
    dim[1] = 0;
    dim[2] = 0;
    dim[3] = 0;
    frac_flag = 0;
    frac_len = 0.0;
    frac_angle = 0.0;
    frac_skin = 0.0;
    frac_perm = 0.0;
    frac_width = 0.0;

    connection_type = 0;
    start_time = -1;
    is_activated = 0;
    ignore_fracture = 0;
  }

  std::string &fracture::tname()
  {
    return name;
  }

  int fracture::check_data() const
    {
      if (dim[0] < 0 || dim[1] < 0 || dim[2] < 0 || dim[3] < 0)
        return -2;
      if (frac_len < 0 || frac_perm < 0 || frac_width < 0)
        return -3;
      return 0;
    }
}
