/**
 * \file pvt_dummy_iface.h
 * \brief dummy holder for SCAL data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */

#ifndef BS_PVT_DUMMY_IFACE_H_
#define BS_PVT_DUMMY_IFACE_H_

#include "table_iface.h"

namespace blue_sky {


class pvt_dummy_iface : public objbase
  {
    public:

    virtual
    ~pvt_dummy_iface () {}

    virtual void
    init (int pvt_type) = 0;

    virtual std::list <BS_SP( table_iface)>
    get_table () const = 0;

    virtual boost::python::list
    py_get_table () const = 0;

    virtual int
    get_pvt_type () const = 0;

  };

}

#endif //BS_PVT_DUMMY_IFACE_H_
