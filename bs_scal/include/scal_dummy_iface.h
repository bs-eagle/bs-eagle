/**
 * \file scal_dummy_iface.h
 * \brief dummy holder for SCAL data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */
 
#ifndef BS_SCAL_SCAL_DUMMY_IFACE_H_
#define BS_SCAL_SCAL_DUMMY_IFACE_H_

#include "table_iface.h"

namespace blue_sky {

  
class scal_dummy_iface : public objbase
  {
    public:
    
    virtual 
    ~scal_dummy_iface () {}
    
    
    virtual std::pair <BS_SP( table_iface), BS_SP( table_iface)>
    get_table () const = 0;
    
    virtual boost::python::list 
    py_get_table () const = 0;
  
  };

}

#endif //BS_SCAL_SCAL_DUMMY_IFACE_H_