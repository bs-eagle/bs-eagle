/**
 * \file pvt_dummy.h
 * \brief dummy holder for SCAL data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */
#ifndef BS_PVT_DUMMY_H_
#define BS_PVT_DUMMY_H_

#include "pvt_dummy_iface.h"
#include "table_iface.h"


namespace blue_sky {

  
  class BS_API_PLUGIN pvt_dummy : public pvt_dummy_iface
  {
    public:
    virtual ~pvt_dummy () {}
    
    void init (int pvt_type);
    
    
    BS_SP( table_iface)
    get_table () const;
    
    int
    get_pvt_type () const;
    
    private:
      BS_SP (table_iface) pvt_data;
      int                 pvt_type;
    


    public:
    BLUE_SKY_TYPE_DECL (pvt_dummy);
  };

}






#endif //BS_PVT_DUMMY_H_