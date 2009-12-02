/**
 * \file main.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 23.03.2009
 * */
#include "bs_matrix_stdafx.h"

#include "shared_vector.h"
#include "matrix_base.h"
#include "bcsr_matrix.h"
#include "read_b_matrix.h"
#include "jacobian_matrix.h"

#include "py_bs_assert.h"
#include "py_matrix_base.h"
#include "py_bcsr_matrix.h"
#include "py_naive_file_reader.h"
#include "py_save_seq_vector.h"

#include "py_jacobian_matrix.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_mx", "1.0.0", "Base matrixes for blue_sky", "Base matrixes for blue_sky", "bs_mx")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type(*bs_init.pd_, matrix_base<shared_vector<float>, shared_vector<int> >::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, matrix_base<shared_vector<double>, shared_vector<int> >::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix<shared_vector<float>, shared_vector<int> >::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix<shared_vector<double>, shared_vector<int> >::bs_type()); BS_ASSERT (res);

    res &= blue_sky::jacobian_matrix_register_type (pd); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  //Python vector types
  class_<shared_vector<float> >("vector_float")
    .def(vector_indexing_suite<shared_vector<float> >())
    ;

  class_<shared_vector<double> >("vector_double")
    .def(vector_indexing_suite<shared_vector<double> >())
    ;

  class_ <shared_vector <int> > ("vector_int")
    .def (vector_indexing_suite <shared_vector <int> > ())
    ;

  python::py_export_naive_file_reader ();
  python::py_export_save_seq_vector ();

  python::py_export_base_matrices ();
  python::py_export_jacobian_matrix ();

  python::py_export_read_b_matrix ();

  
}
#endif

