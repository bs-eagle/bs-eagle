#ifndef ARITHMETIC_STACK_H
#define ARITHMETIC_STACK_H

#include "ar_args.h"
#include <stack>

#define ST_LEN 512

namespace blue_sky
  {
  template <typename ar_smth_t>
  class BS_API_PLUGIN ar_stack
    {
      // typedefs
      //typedef ar_args <strategy_t>  arg_t;
      typedef std::stack <ar_smth_t>    stk_t;
      //typedef ar_stack <strategy_t> this_t;

    public:
      ar_stack ();
      ~ar_stack ();

      void push (const ar_smth_t &arg);
      ar_smth_t pop ();
			bool is_empty () const;

      //! data
    protected:
      stk_t st;
    };
}

#endif // ARITHMETIC_STACK_H
