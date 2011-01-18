#include "bs_bos_core_data_storage_stdafx.h"

#include "ar_stack.h"
#include "ar_args.h"
#include "ar_operat.h"

namespace blue_sky
  {
  template <typename ar_smth_t>
  ar_stack<ar_smth_t>::ar_stack ()
  {}

  template <typename ar_smth_t>
  ar_stack<ar_smth_t>::~ar_stack ()
  {
    for (;!st.empty ();)
      st.pop ();
  }

  template <typename ar_smth_t>
  void ar_stack<ar_smth_t>::push (const ar_smth_t &arg)
  {
    st.push(arg);
  }

  template <typename ar_smth_t>
  ar_smth_t ar_stack<ar_smth_t>::pop ()
  {
    if (st.empty())
      {
				BS_ASSERT("stack is empty");
        throw bs_exception("ar_stack pop","stack is empty");
      }
    ar_smth_t arg = st.top();
    st.pop();
    return arg;
  }

	template <typename ar_smth_t>
	bool ar_stack<ar_smth_t>::is_empty () const {
		return st.empty();
	}

  template class ar_stack < ar_args <base_strategy_fif> >;
  template class ar_stack < ar_args <base_strategy_did> >;
	template class ar_stack < ar_args <base_strategy_dif> >;
  template class ar_stack < ar_operat <base_strategy_fif> >;
  template class ar_stack < ar_operat <base_strategy_did> >;
	template class ar_stack < ar_operat <base_strategy_dif> >;
}
