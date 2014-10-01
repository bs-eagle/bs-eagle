/**
 * \file scal_test_model.h
 * \brief model data holder. for testing only.
 * \author Sergey Miryanov
 * \date 23.05.2008
 * */
#ifndef BS_SCAL_TEST_MODEL_H_
#define BS_SCAL_TEST_MODEL_H_

#include "bs_object_base.h"
#include "conf.h"

namespace blue_sky
  {
  namespace scal
    {

    class BS_API_PLUGIN test_model : public objbase
      {

      public:
        typedef t_double                    item_t;
        typedef test_model                  this_t;

        void init (int n_phases, int cell_count);

        item_t *get_relative_perm (int cell_index);
        item_t *get_s_deriv_relativ_perm (int cell_index);

        item_t *get_cap_pressure (int cell_index);
        item_t *get_s_deriv_cap_pressure (int cell_index);

        item_t *get_sat (int cell_index);

        int get_cell_count () const;

      private:

        typedef std::vector <double>                 vector_t;

        int         n_phases;
        int         cell_count;

        vector_t    relative_perm;
        vector_t    s_deriv_relativ_perm;
        vector_t    cap_pressure;
        vector_t    s_deriv_cap_pressure;
        vector_t    sat;

      public:

        BLUE_SKY_TYPE_DECL_T (test_model);

      };

  } // namespace scal
} // namespace blue_sky

#endif // #ifndef BS_SCAL_TEST_MODEL_H_
