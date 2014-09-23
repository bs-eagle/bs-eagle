/**
* \file   py_strategies.cpp
* \brief  Python wrapper for strategies
* \author Miryanov Sergey
* \date 2008-04-04
*/

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_strategies.h"

using namespace boost::python;

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_strategies ()
  {
    using namespace boost::python;

    class_ <strategy_iface> ("strategy_iface");

    class_ <base_strategy_fif, bases <strategy_iface> > ("strategy_fif");
    class_ <base_strategy_did, bases <strategy_iface> > ("strategy_did");
    class_ <base_strategy_dif, bases <strategy_iface> > ("strategy_dif");
    class_ <base_strategy_flf, bases <strategy_iface> > ("strategy_flf");
    class_ <base_strategy_dld, bases <strategy_iface> > ("strategy_dld");
    class_ <base_strategy_dlf, bases <strategy_iface> > ("strategy_dlf");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
