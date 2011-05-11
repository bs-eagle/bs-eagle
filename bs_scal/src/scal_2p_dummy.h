/**
 * \file scal_2p_dummy.h
 * \brief dummy holder for SCAL data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */
#ifndef BS_SCAL_SCAL_2P_DUMMY_H_
#define BS_SCAL_SCAL_2P_DUMMY_H_

#include "scal_dummy_iface.h"
#include "table_iface.h"


namespace blue_sky {

  
  class BS_API_PLUGIN scal_2p_dummy : public scal_dummy_iface
  {
    public:
    virtual ~scal_2p_dummy () {}
    
    
    std::pair <BS_SP( table_iface), BS_SP( table_iface)>
    get_table () const;
    
    private:
      BS_SP (table_iface) scal_data;
    


    public:
    BLUE_SKY_TYPE_DECL (scal_2p_dummy);
  };


  class BS_API_PLUGIN scal_3p_dummy : public scal_dummy_iface
  {
    public:
    virtual ~scal_3p_dummy () {}
    
    
    std::pair <BS_SP( table_iface), BS_SP( table_iface)>
    get_table () const;
    
    private:
      BS_SP (table_iface) scal_data1;
      BS_SP (table_iface) scal_data2;


    public:
    BLUE_SKY_TYPE_DECL (scal_3p_dummy);
  };
}






#endif //BS_SCAL_SCAL_2P_DUMMY_H_