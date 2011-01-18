#include "bs_bos_core_data_storage_stdafx.h"
#include "data_class.h"
#include "arrays.h"
#include "main_def.h"

#include <sstream>

// Definition of code number of tokens
//! index for token argument
#define TOKEN_ARG               0
//! index for token function
#define TOKEN_FUNC              1
//! index for token open bracket
#define TOKEN_OPEN_BRACKET      2
//! index for token close bracket
#define TOKEN_CLOSE_BRACKET     3
//! index for token PLUS
#define TOKEN_PLUS              4
//! index for token MINUS
#define TOKEN_MINUS             5
//! index for token MULT
#define TOKEN_MULT              6
//! index for token DEVD
#define TOKEN_DEVD              7
//! index for token COMMA
#define TOKEN_COMMA             8
//! index for initial token all available
#define TOKEN_ALL               9

//! Size of available token table
#define SIZE_TOKEN_TABLE        10
//! Definition of index in token array
#define PREVIOS_TOKEN           0
//! Definition of index in token array
#define CURRENT_TOKEN           1
//! Depth of stack
#define AR_STACK_DEPTH 100

//  Definitions priority for operators in arithmetic
//! priority of function
#define OP_PREORITY_FUNC  4
//! priority of plus
#define OP_PREORITY_ADD  1
//! priority of multiplication
#define OP_PREORITY_MULT  2
//! priority of bracket
#define OP_PREORITY_BRACKET  5

namespace blue_sky {
using namespace pool;
using namespace std;

template <class strategy_t>
  void idata<strategy_t>::build_argument_list (const sp_reader_t &reader) {
	ostringstream out_s;

	int nb = nx * ny * nz;
	if (!nb)
	{
		out_s << "Error in " << reader->get_prefix ()
			<< ": you should specify keyword 'DIMENS' before keyword 'ARITHMETIC'";
		BS_ASSERT(false) (out_s.str());
		throw bs_exception("idata::build_argument_list",out_s.str().c_str());
	}

	BOSOUT (section::arithmetic, level::debug) << "i_pool size = " << i_pool()->size () << bs_end;
	bs_node::n_iterator p_arr, end;
	int nlen;
	size_triple sz;
	for(p_arr = i_pool()->begin(), end = i_pool()->end(); p_arr != end; ++p_arr) {
		BOSOUT (section::arithmetic, level::debug) << "i_pool contains element " << p_arr->name() << bs_end;
		nlen = i_pool()->get_dimens (p_arr->name(), sz);
		ar_args_t a(p_arr->name(), ar_args_t::INT_T, carray_i(p_arr->name()), sz.nx, sz.ny, sz.nz);
		args.insert ( std::pair< std::string, ar_args_t >(a.get_name (),a) );
	}

	BOSOUT (section::arithmetic, level::debug) << "d_pool size = " << d_pool()->size () << bs_end;
	for(p_arr = d_pool()->begin(), end = d_pool()->end(); p_arr != end; ++p_arr) {
		BOSOUT (section::arithmetic, level::debug) << "d_pool contains element " << p_arr->name() << bs_end;
		nlen = d_pool()->get_dimens (p_arr->name(), sz);
		ar_args_t a(p_arr->name(), ar_args_t::FPOINT_T, carray_fp(p_arr->name()), sz.nx, sz.ny, sz.nz);
		args.insert ( std::pair< std::string, ar_args_t > (a.get_name (),a) );
        index_t nlen1 = 0;
	}
}

	template <class strategy_t>
	void idata<strategy_t>::build_operator_list () {
		ops.insert(ar_operat_t(&ar_operation_abs<strategy_t>, OP_PREORITY_FUNC, 1, "ABS"));
		ops.insert(ar_operat_t(&ar_operation_exp<strategy_t>, OP_PREORITY_FUNC, 1, "EXP"));
		ops.insert(ar_operat_t(&ar_operation_log<strategy_t>, OP_PREORITY_FUNC, 1, "LOG"));
		ops.insert(ar_operat_t(&ar_operation_log10<strategy_t>, OP_PREORITY_FUNC, 1, "LOG10"));
		ops.insert(ar_operat_t(&ar_operation_sqrt<strategy_t>, OP_PREORITY_FUNC, 1, "SQRT"));
		ops.insert(ar_operat_t(&ar_operation_tan<strategy_t>, OP_PREORITY_FUNC, 1, "TAN"));
		ops.insert(ar_operat_t(&ar_operation_sin<strategy_t>, OP_PREORITY_FUNC, 1, "SIN"));
		ops.insert(ar_operat_t(&ar_operation_cos<strategy_t>, OP_PREORITY_FUNC, 1, "COS"));
		ops.insert(ar_operat_t(&ar_operation_plus<strategy_t>, OP_PREORITY_ADD, 2, "+"));
		ops.insert(ar_operat_t(0, OP_PREORITY_BRACKET, 1, "("));
		ops.insert(ar_operat_t(0, OP_PREORITY_BRACKET, 1, ")"));
		ops.insert(ar_operat_t(&ar_operation_minus<strategy_t>, OP_PREORITY_ADD, 2, "-"));
		ops.insert(ar_operat_t(&ar_operation_mult<strategy_t>, OP_PREORITY_MULT, 2, "*"));
		ops.insert(ar_operat_t(&ar_operation_dev<strategy_t>, OP_PREORITY_MULT, 2, "/"));
		ops.insert(ar_operat_t(&ar_operation_rminus<strategy_t>, OP_PREORITY_FUNC, 1, "RM"));
		ops.insert(ar_operat_t(&ar_operation_max<strategy_t>, OP_PREORITY_FUNC, 2, "MAX"));
		ops.insert(ar_operat_t(&ar_operation_min<strategy_t>, OP_PREORITY_FUNC, 2, "MIN"));
		//ops.insert(ar_operat_t(&ar_operation_rand, OP_PREORITY_FUNC, 2, "RAND"));
	}

  template <class strategy_t>
  void idata<strategy_t>::clear_argument_list ()
  {
    args.clear ();
    //iargs.clear ();
  }


	template <class strategy_t>
	void idata<strategy_t>::output_argument_list () {
		//std::ostringstream ostr_s;
		std::ofstream ostr_s ("arithmetic.arrays", std::ios_base::out | std::ios_base::trunc);
		for (typename std::map <std::string, ar_args_t>::iterator itr = args.begin(); itr != args.end(); ++itr) {
			ostr_s.clear();
			ostr_s << "\nname = " << itr->second.get_name ()
						<< "\nni = " << itr->second.ni.data ()
						<< "\nnj = " << itr->second.nj.data ()
						<< "\nnk = " << itr->second.nk.data ()
						<< "\ni1 = " << itr->second.i1.data ()
						<< "\ni2 = " << itr->second.i2.data ()
						<< "\nj1 = " << itr->second.j1.data ()
						<< "\nj2 = " << itr->second.j2.data ()
						<< "\nk1 = " << itr->second.k1.data ()
						<< "\nk2 = " << itr->second.k1.data ()
						<< "\nn_size = " << itr->second.n_size.data ()
						<< "\narr_type = " << itr->second.get_type () << "\n";

			for (int k = 0; k < itr->second.nk; ++k)
				for (int j = 0; j < itr->second.nj; ++j) {
					for (int i = 0; i < itr->second.ni; ++i) {
						int c = i + j * itr->second.ni + k * itr->second.ni * itr->second.nj;
						if (itr->second.get_type () == ar_args_t::FPOINT_T)
							ostr_s << ((item_t*)(itr->second.get_array ()))[c] << " ";
						else
							ostr_s << ((int*)(itr->second.get_array ()))[c] << " ";
					}
					ostr_s << "\n";
				}

			//BOSOUT (section::arithmetic, level::debug) << ostr_s.str () << bs_end;
		}
	}

  /*!
  	\brief delete all blanks from string 'buf'
  	\param buf  -- pointer to string
  	\return if success                                      0
  	\return if bad pointer buf                              -1
  */
  template <class strategy_t>
  int idata<strategy_t>::no_blanks (char *buf) const
    {
      int blanks_counter = 0;
      int pos = 0;
      // check input variable
      if (!buf)
        {
          bs_throw_exception ("Internal error (buf is null)");
        }
      do
        {
          if (buf[pos] == ' ' || buf[pos] == '\t')   // if blank skip it
            {
              ++blanks_counter;
              ++pos;
            }
          else                      // copy symbol
            {
              buf[pos - blanks_counter] = buf[pos];
              ++pos;
            }
        }
      while (buf[pos] != '\0');
      buf[pos - blanks_counter] = buf[pos];
      return YS_SUCCESS;
    }

  /*!
  	\brief read left side argument and the limits to change array idata

  *  format of left hand part:      argument(i1:i2,j1:j2,k1:k2) or
  *                                       argument(:i2,j1:j2,k1:k2) or
  *                                       argument(i1:,j1:j2,k1:k2) or
  *                                       argument(,j1:j2,k1:k2) or
  *                                       argument(:,i1:,j1:j2,k1:k2) or
  *                                       and ...
  *       do not forget free memory from pointer 'name'
  	\param s  -- pointer to string
  	\param right
  	\param name
  	\param i1
  	\param i2
  	\param j1
  	\param j2
  	\param k1
  	\param k2
  	\param keyword
  	\return if success                                      0
  	\return if bad pointer                                  -1
  	\return if can not allocate memory                      -2
  */
  template <typename strategy_t>
  void idata<strategy_t>::read_left (const sp_reader_t &reader,
                                     char *s, char **right, std::string &name, int *i1,
                                     int *i2, int *j1, int *j2, int *k1, int *k2, const char *keyword)
  {
    std::ostringstream out_s;
    //local variables
    int pos = 0;
    int c = 0;                 // c - position in buf
    char buf[256];
    char *end;
    int fld = 1;                  // flag for dimention if 1 -- first dimention
    //                       2 -- second
    //                       3 -- third
    int fl_1 = 0;                 // flag for first var read      0 -- first var have not read
    //                              1 -- first var have read
    int fl_2 = 0;                 // flag for second var read     0 -- second var have not read
    //                              1 -- second var have read
    int fl_dp = 0;                // fl_dp == 1 if read part of string have symbol ':'
    // and zero if have not
    int *d1[3];                   // pointer to values of first limit in dimension
    int *d2[3];                   // pointer to values of second limit in dimension

    // initialize pointers
    d1[0] = i1;
    d1[1] = j1;
    d1[2] = k1;
    d2[0] = i2;
    d2[1] = j2;
    d2[2] = k2;

    // initialize section
    *i1 = *j1 = *k1 = 0;
    *i2 = *j2 = *k2 = 0;

    // Check section
    // Check all input variables
    if (s == 0 || i1 == 0 || i2 == 0 || j1 == 0 || j2 == 0 ||
        k1 == 0 || k2 == 0)
      {
        out_s << "Internal error in " << reader->get_prefix ()
        << ": keyword " << keyword;
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::read_left",out_s.str().c_str());
      }
    // count number character in left arguments
		while (s[pos] != '\0' && s[pos] != '=' && s[pos] != '(') {
			name += s[pos];
      ++pos;
		}
    if (s[pos] == '\0')
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": assignement '=' not found in keyword " << keyword;
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::read_left",out_s.str().c_str());
      }

    // check string for dimension limits
    if (s[pos] == '=')            // if no dimension limits
      {
        *right = s + pos + 1;     // set pointer to right hand side
        return;
      }
    else
      {
        for (;;)
          {
            ++pos;                // goto the next symbol
            switch (s[pos])
              {
              case ')':           // if next symbol == ')'
              {
                if (fl_dp == 0) // if last part have not ':'
                  {
                    if (fl_1 == 0)      // if first limit have not read
                      {
                        // set limits to default
                        *d1[fld - 1] = 0;
                        *d2[fld - 1] = 0;
                      }
                    else
                      {
                        buf[c] = '\0';  // finish reading limit value
                        *d1[fld - 1] = strtol (buf, &end, 10);  // convert and set limit value
                        *d2[fld - 1] = *d1[fld - 1];
                        // check for format error
                        if (buf == end) // format error
                          {
                            out_s << "Error in " << reader->get_prefix ()
                            << ": cannot parse integer argument " << buf
                            << " for keyword " << keyword;
                            BS_ASSERT(false) (out_s.str());
                            throw bs_exception("idata::read_left",out_s.str().c_str());
                          }
                      }
                  }
                else            // last part have ':'
                  {
                    if (fl_2 == 0)      // if second limit have not read set it to default
                      {
                        *d2[fld - 1] = 0;
                      }
                    else        // if second limit have read
                      {
                        buf[c] = '\0';  // finish reading limit value
                        *d2[fld - 1] = strtol (buf, &end, 10);
                        // check end
                        if (buf == end) // format error
                          {
                            out_s << "Error in " << reader->get_prefix ()
                            << ": cannot parse integer argument " << buf
                            << " for keyword " << keyword;
                            BS_ASSERT(false) (out_s.str());
                            throw bs_exception("idata::read_left",out_s.str().c_str());
                          }
                      }
                  }

                if (s[++pos] == '=')    // if next symbol is '=' finish reading else format error
                  {
                    *right = s + pos + 1;
                    return;
                  }
                else
                  {
                    out_s << "Error in " << reader->get_prefix ()
                    << ": assignement '=' not found in keyword "
                    << keyword << " argument " << s;
                    BS_ASSERT(false) (out_s.str());
                    throw bs_exception("idata::read_left",out_s.str().c_str());
                  }
                break;
              }
              case ',':           // end of part
              {
                if (fld > 2)    // Too many ',' in input string. Format error
                  {
                    out_s << "Error in " << reader->get_prefix ()
                    << ": too many commas found in keyword " << keyword
                    << " argument " << s;
                    BS_ASSERT(false) (out_s.str());
                    throw bs_exception("idata::read_left",out_s.str().c_str());
                  }
                if (fl_dp == 0) // if no symbol ':' in last part
                  {
                    if (fl_1 == 0)      // if first limit have not read
                      {
                        *d1[fld - 1] = 0;
                        *d2[fld - 1] = 0;
                      }
                    else        // if first limit have read
                      {
                        buf[c] = '\0';
                        *d1[fld - 1] = strtol (buf, &end, 10);
                        *d2[fld - 1] = *d1[fld - 1];
                        // check end
                        if (buf == end) // format error
                          {
                            out_s << "Error in " << reader->get_prefix ()
                            << ": cannot parse integer argument " << buf
                            << " for keyword " << keyword;
                            BS_ASSERT(false) (out_s.str());
                            throw bs_exception("idata::read_left",out_s.str().c_str());
                          }
                      }
                  }
                else            // if last part have symbol ':'
                  {
                    if (fl_2 == 0)      // if second limit have not read
                      {
                        *d2[fld - 1] = 0;
                      }
                    else
                      {
                        buf[c] = '\0';
                        *d2[fld - 1] = strtol (buf, &end, 10);
                        // check end
                        if (buf == end) // format error
                          {
                            out_s << "Error in " << reader->get_prefix ()
                            << ": cannot parse integer argument " << buf
                            << " for keyword " << keyword;
                            BS_ASSERT(false) (out_s.str());
                            throw bs_exception("idata::read_left",out_s.str().c_str());
                          }
                      }
                  }
                fl_1 = 0;       // set flag values to zero and wait for next part
                fl_2 = 0;
                fl_dp = 0;
                ++fld;
                c = 0;          //??????
                break;
              }
              case ':':           // first and second limit separator
                // found
              {
                if (fl_dp == 0)
                  {
                    if (fl_1 == 0)      // if first limit have not read
                      {
                        *d1[fld - 1] = 0;
                        fl_dp = 1;
                        c = 0;  //??????
                      }
                    else
                      {
                        buf[c] = '\0';
                        *d1[fld - 1] = strtol (buf, &end, 10);
                        // check end
                        if (buf == end) // format error
                          {
                            out_s << "Error in " << reader->get_prefix ()
                            << ": cannot parse integer argument " << buf
                            << " for keyword " << keyword;
                            BS_ASSERT(false) (out_s.str());
                            throw bs_exception("idata::read_left",out_s.str().c_str());
                          }
                        c = 0;
                        fl_dp = 1;
                      }
                  }
                else            // Error in string format
                  {
                    out_s << "Error in " << reader->get_prefix ()
                    << ": too many ':' found in keyword " << keyword
                    << " argument " << s;
                    BS_ASSERT(false) (out_s.str());
                    throw bs_exception("idata::read_left",out_s.str().c_str());
                  }
                break;
              }
              case '\0':
              {
                out_s << "Error in " << reader->get_prefix ()
                << ": unknown end of string " << s;
                BS_ASSERT(false) (out_s.str());
                throw bs_exception("idata::read_left",out_s.str().c_str());
                break;
              }
              default:            // copy symbol to bufer
              {
                buf[c] = s[pos];
                ++c;
                if (fl_dp == 0 && fl_1 == 0)
                  fl_1 = 1;
                else if (fl_dp == 1 && fl_2 == 0)
                  fl_2 = 1;
                break;
              }
              }
          }
      }
  }


  /*!
  	\brief main method to do arithmetic action

  *        call No_Blanks to delete blanks from input string.\n
  *        call Read_Left to read left side of string.\n
  *        set all limits for left side argument.\n
  *        call Calculate to parse rite side string.
  	\param buf  -- pointer to string
  	\param keyword
  	\return if success                      YS_SUCCESS
  	\return if error                         < 0
  */

template <class strategy_t>
  void idata<strategy_t>::algorithm_read_and_done (const sp_reader_t &reader, char *buf, const char *keyword)
{
	// local variables
	char * right = 0;
	std::string name;
	int i1,i2,j1,j2,k1,k2;
	ostringstream out_s;
	//std::listar_args_t l

	// check input variables;
	if (!buf || !keyword)
	{
		out_s << "Internal error in "
			<< reader->get_prefix () << ": keyword "
			<< keyword;
		BS_ASSERT(false) (out_s.str());
		throw bs_exception("idata::algorithm_read_and_done",out_s.str().c_str());
	}
	if (strlen (buf) < 3)
	{
		out_s << "Error in " << reader->get_prefix ()
			<< ": unrecognized argument " << buf
			<< " for keyword " << keyword;
		BS_ASSERT(false) (out_s.str());
		throw bs_exception("idata::algorithm_read_and_done",out_s.str().c_str());
	}
	// delete blanks from input string <buf>
	if (no_blanks (buf) < 0)   // return string with out blanks
	{
		out_s << "Internal error in " << reader->get_prefix ()
			<< ": keyword " << keyword;
		BS_ASSERT(false) (out_s.str());
		throw bs_exception("idata::algorithm_read_and_done",out_s.str().c_str());
	}

	read_left (reader,buf,&right,name,&i1,&i2,&j1,&j2,&k1,&k2,keyword);

	typename std::map <std::string,ar_args_t>::iterator l = args.find (name);
	if (l == args.end ())
	{
		out_s << "Error in " << reader->get_prefix ()
			<< ": unknown token " << name;
		BS_ASSERT(false) (out_s.str());

		BOSOUT (section::arithmetic, level::debug) << "Token " << name << " adding..." << bs_end;

		size_triple sz;
		int nlen;
		if (array_sizes_i().find(name) != array_sizes_i().end()) {
			i_pool()->set_array(name, array_info_i(at(array_sizes_i(), name), at(array_initv_i(), name)));
			nlen = i_pool()->get_dimens (name, sz);
			ar_args_t a(name, ar_args_t::INT_T, i_pool()->carray(name), sz.nx, sz.ny, sz.nz);
			args.insert ( std::pair <std::string,ar_args_t> (a.get_name (),a) );
			l = args.find (name);
		}
		else if (array_sizes_fp().find(name) != array_sizes_fp().end()) {
			i_pool()->set_array(name, array_info_fp(at(array_sizes_fp(), name), at(array_initv_i(), name)));
			nlen = i_pool()->get_dimens (name, sz);
			ar_args_t a(name, ar_args_t::FPOINT_T, d_pool()->carray(name), sz.nx, sz.ny, sz.nz);
			args.insert ( std::pair <std::string,ar_args_t> (a.get_name (),a) );
			l = args.find (name);
		}
		else
			throw bs_exception("idata::algorithm_read_and_done",out_s.str().c_str());
	}
	ar_args_t &larg = l->second;
	// check what array have defined
	if (!larg.get_array ())         // if not allocate memory for array
	{
		// try to allocate memory if ni,nj,nk are greate when zero
		if (larg.allocate () != YS_SUCCESS)
		{
			out_s << "Error in " << reader->get_prefix ()
				<< ": cannot allocate memory";
			BS_ASSERT(false) (out_s.str());
			throw bs_exception("idata::algorithm_read_and_done",out_s.str().c_str());
		}
		// set pointer in idata class to allocated area
		//if (l->set_ptr_data () != YS_SUCCESS)
		//  return -9;
		//}
}
// set limits of calculation
// if some limits equals zero, set it to default (max area)
if (i1)
	larg.i1.data () = i1 - 1;
	else
	larg.i1.data () = 0;
if (i2)
	larg.i2.data () = i2;
	else
	larg.i2.data () = larg.ni;
if (j1)
	larg.j1.data () = j1 - 1;
	else
	larg.j1.data () = 0;

if (j2)
	larg.j2.data () = j2;
	else
	larg.j2.data () = larg.nj;
if (k1)
	larg.k1.data () = k1 - 1;
	else
	larg.k1.data () = 0;

if (k2)
	larg.k2.data () = k2;
	else
	larg.k2.data () = larg.nk;

if (larg.i2 > larg.ni)
	larg.i2.data () = larg.ni;
if (larg.j2 > larg.nj)
	larg.j2.data () = larg.nj;
if (larg.k2 > larg.nk)
	larg.k2.data () = larg.nk;

	calculate (reader, larg, right, keyword);    // parse right side, calculate it and copy result to left side argument
	return;
	}

  /*!
  	\brief try to read float value from buf
  	\param reader
  	\param ar
  	\param buf
  	\param keyword
  */
  template <typename strategy_t>
  int idata<strategy_t>::read_numbers (ar_stack <ar_args_t> &st, char **buf)
  {
    typedef typename strategy_t::fp_type_t item_t;

    char *pbeg = 0,*pend = 0;
    item_t t;

    pbeg = *buf;
    // try to read number from string
    t = strtod (pbeg, &pend);
    if (pbeg != pend)
      {
        // set token code
        ar_tokens [CURRENT_TOKEN] = TOKEN_ARG;
        // convert number to ar_args_t class
        st.push (ar_args_t (ar_args_t::FPOINT_T,&t));
      }
		else
			return 0;

    *buf = pend;
		return 1;
  }

  template <typename strategy_t>
  void idata<strategy_t>::calculate (const sp_reader_t &reader, ar_args_t &res, char *buf, const char *keyword)
  {
		std::ostringstream out_s;
    // variables for working with string
    char *p_strt = 0, *p_end = 0;

    //check all input pointers
    if (!buf)
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": not enough valid arguments for keyword " << keyword;
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::calculate",out_s.str().c_str());
      }

    ar_stack <ar_args_t>   ar;
    ar_stack <ar_operat_t> op;

    int arguments_count[AR_STACK_DEPTH];  // Array to check arguments
    // number
    int flag_brek = 0;            // counter of brakets

    p_end = buf;

    // Init arithmetic token codes
    ar_tokens[PREVIOS_TOKEN] = TOKEN_ALL;
    ar_tokens[CURRENT_TOKEN] = TOKEN_ALL;

    while (*p_end != '\0')        // main loop for parsing right side
      {
        p_strt = p_end;
        // 1) try t read operators
        if (!read_operator (reader, res, ar, op, &p_end, &flag_brek, keyword, buf, arguments_count)) {
					// 2) try to read numbers
					if (!read_numbers (ar, &p_end)) {
						//	 3) try to read arguments or function
						if (!read_arg_func (reader, ar, op, &p_end, keyword, arguments_count, flag_brek)) {

							//	 t	ry to read symbol ','
							if (*p_end == ',')
							{
								// move pointer to the next token
								++p_end;
								// se		t token code
								ar_tokens[CURRENT_TOKEN] = TOKEN_COMMA;
								// decrime		nt arguments counter which need for function
								--arguments_count[flag_brek];
								// check number 		of arguments
								if (arguments_count[flag_brek] < 0)
								{
									out_s << "Error in " << reader->get_prefix ()
												<< ": too many arguments in expression " << buf;
									BS_ASSERT(false) (out_s.str());
									throw bs_exception("idata::calculate",out_s.str().c_str());
								}
							}
							else
							{
								out_s << "Error in " << reader->get_prefix ()
											<< ": invalid expression " << buf;
								BS_ASSERT(false) (out_s.str());
								throw bs_exception("idata::calculate",out_s.str().c_str());
							}
						}
					}
				}

        // check sequence of token (previos and current)
        if (!(test_token (ar_tokens[PREVIOS_TOKEN], ar_tokens[CURRENT_TOKEN])))
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": unknown sequence of tokens in expression " << buf;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::calculate",out_s.str().c_str());
          }

				ar_tokens[PREVIOS_TOKEN] = ar_tokens[CURRENT_TOKEN];
      }

    // after main loop calculate all
		BOSOUT (section::arithmetic, level::debug) << "do_calculate start!" << bs_end;
    if (do_calculate (reader, res, ar, op, keyword))
      {
        bs_throw_exception (boost::format ("Error in %s: do_operator failed on keyword %s")
          % reader->get_prefix () % keyword);
      }
		BOSOUT (section::arithmetic, level::debug) << "do_calculate end!" << bs_end;

    // take result from argument stack
		if (!ar.is_empty ()) {
			ar_args_t ws (ar.pop ());
		
			BOSOUT (section::arithmetic, level::debug) << "some checks begin!" << bs_end;
	    // if result have not been initialized print error message
		  if (ws.def_flag)
      {
        bs_throw_exception (boost::format ("Error in %s: value of token %s is undefined")
          % reader->get_prefix () % ws.get_name ());
      }

	    if (ws.ni * ws.nj * ws.nk > 1)
      {
        if (res.i2 > ws.ni || res.j2 > ws.nj || res.k2 > ws.nk)
          {
            bs_throw_exception (boost::format ("Error in %s: index in rhs argument is out of range")
              % reader->get_prefix () % ws.get_name ());
          }
      }

		  // copy all idata from result to left hand argument
			for (int k = res.k1; k < res.k2; ++k)
				for (int j = res.j1; j < res.j2; ++j)
					for (int i = res.i1; i < res.i2; ++i)
					{
            if (ws.ni * ws.nj * ws.nk > 1)
              {
                A_SET_VALUE (res, i + j * res.ni + k * res.ni * res.nj,
                             A_GET_VALUE (ws, i + j * res.ni + k * res.ni * res.nj));
              }
            else
              {
                A_SET_VALUE (res, i + j * res.ni + k * res.ni * res.nj, A_GET_VALUE (ws, 0));
              }
          }

			// set left side argument as initialized
			res.def_flag = 0;

			// check what stack of arguments is empty
			if (!ar.is_empty ()) {
				ws = ar.pop ();

				out_s << "Error in " << reader->get_prefix ()
					<< ": cannot parse expression " << buf;
				BS_ASSERT(false) (out_s.str());
				throw bs_exception("idata::calculate",out_s.str().c_str());
			}

			// check brackets balans
			if (flag_brek)
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": unbalanced brackets in expression " << buf;
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::calculate",out_s.str().c_str());
      }
		} else {
			out_s << "Error in " << reader->get_prefix ()
				<< ": no ws token!";
			BS_ASSERT(false) (out_s.str());
			throw bs_exception("idata::calculate",out_s.str().c_str());
		}
  }

  /*!
  	\brief
  	\param Previos
  	\param Current
  	\return
  */

  template <typename strategy_t>
  int idata <strategy_t>::test_token (int prev, int cur)
  {
    // if value in table are 0 this sequence of previos and current token are not available
    // if value in table are 1 all is OK
    static const int available_token[SIZE_TOKEN_TABLE][SIZE_TOKEN_TABLE] =
    {
      /*prev tok\ curr tok  ARG     Func    '('     ')'     '+'     '-'     '*'     '/'     ','     All            */
      /* Arg  */ {0, 0, 0, 1, 1, 1, 1, 1, 1, 0},
      /* Func */ {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      /* '('  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* ')'  */ {0, 0, 0, 1, 1, 1, 1, 1, 1, 0},
      /* '+'  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* '-'  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* '*'  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* '/'  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* ','  */ {1, 1, 1, 0, 1, 1, 0, 0, 0, 0},
      /* All  */ {1, 1, 1, 1, 1, 1, 1, 1, 1, 0}
    };
    return available_token[prev][cur];
  }

  template <typename strategy_t>
  int idata <strategy_t>::read_operator (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op,
                                          char **start_ptr, int *flag_brek, const char *keyword,
                                          const char *buf, int *arguments_count)
  {
    std::ostringstream out_s;
    char p_end[2];                // buffer to store operator name
//		ar_operat_t opr;           // temp variable

    // check input variables
    if (!buf || !arguments_count || !keyword)
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": no expression(buf), argument_count or keyword";
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::read_operator",out_s.str().c_str());
      }
    if (!start_ptr || !(*start_ptr))
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": empty expression";
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::read_operator",out_s.str().c_str());
      }
    if (!flag_brek || *flag_brek < 0)
      {
        out_s << "Error in " << reader->get_prefix ()
        << ": no flag_brek";
        BS_ASSERT(false) (out_s.str());
        throw bs_exception("idata::read_operator",out_s.str().c_str());
      }

    p_end[0] = **start_ptr;       //copy operator name
    p_end[1] = '\0';

    if (*p_end != '+' && *p_end != '*' && *p_end != '-' && *p_end != '/' &&
        *p_end != '(' && *p_end != ')')
      {
        // nothing to read (cannot read operator)
        return 0;
      }

    // set current token code
    switch (*p_end)
      {
      case '(':
        ar_tokens[CURRENT_TOKEN] = TOKEN_OPEN_BRACKET;
        break;
      case ')':
        ar_tokens[CURRENT_TOKEN] = TOKEN_CLOSE_BRACKET;
        break;
      case '+':
        ar_tokens[CURRENT_TOKEN] = TOKEN_PLUS;
        break;
      case '-':
        ar_tokens[CURRENT_TOKEN] = TOKEN_MINUS;
        break;
      case '*':
        ar_tokens[CURRENT_TOKEN] = TOKEN_MULT;
        break;
      case '/':
        ar_tokens[CURRENT_TOKEN] = TOKEN_DEVD;
        break;
      default:
        ar_tokens[CURRENT_TOKEN] = TOKEN_ALL;
        break;
      }

    // move pointer to the next token
    ++(*start_ptr);
    // if operator is open bracket add bracket counter
    if (*p_end == '(')
      {
        ++(*flag_brek);
      }
    // if operator is close bracket decrement bracket counter
    if (*p_end == ')')
      {
        // check number of arguments
        if (arguments_count[*flag_brek] > 0)
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": not enough arguments in line " << buf;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
        // decrement bracket counter
        --(*flag_brek);
        // check bracket counter
        if (*flag_brek < 0)
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": unbalanced brackets in expression " << buf;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
      }
    if (*p_end == ')')            // if ) calculate all to (
      {
        if (do_calculate (reader, res, ar, op, keyword))
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": do_calculate fault with expression " << buf;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
      }
    else if (*p_end == '(')       // if current token is '(' add it to operator stack
      {
        typename std::set <ar_operat_t>::iterator l = ops.find (ar_operat_t(p_end));
        if (l == ops.end ())
          {
            out_s << "Internal error in " << reader->get_prefix ()
            << ": keyword " << keyword;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
        // push operator to operator stack
        op.push (*l);
      }
    // check for unary MINUS
    else if ((ar_tokens[PREVIOS_TOKEN] == 2 || ar_tokens[PREVIOS_TOKEN] > 3)
             && *p_end == '-')
      {
        // if current token is unary minus find it in operator list and
        // add to the operator stack
        typename std::set <ar_operat_t>::iterator l = ops.find (ar_operat_t("RM"));
        if (l == ops.end ())
          {
            out_s << "Internal error in " << reader->get_prefix ()
            << ": keyword " << keyword;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
        // push operator to operator stack
        op.push (*l);
      }
    // check for unary PLUS
    else if ((ar_tokens[PREVIOS_TOKEN] == 2 || ar_tokens[PREVIOS_TOKEN] > 3)
             && *p_end == '+')
      {
        // if find it do nothing
      }
    else
      {
        // try to find operator in operator list
        typename std::set <ar_operat_t>::iterator l = ops.find (ar_operat_t(p_end));
        if (l == ops.end ())
          {
            out_s << "Internal error in " << reader->get_prefix ()
            << ": unknown operator " << p_end;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
        ar_operat_t opr (*l);
        // calculate with check preority
        if (prior_calculate (reader, res, ar, op, opr.priority, keyword))
          {
            // if error return error code
            out_s << "Error in " << reader->get_prefix ()
            << ": prior_calculate is fault with keyword " << keyword;
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_operator",out_s.str().c_str());
          }
        // push operator to operator stack
        op.push (opr);
      }
		return 1;
  }

  template <typename strategy_t>
  int idata <strategy_t>::prior_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar,
      ar_stack <ar_operat_t> &op, int priority, const char *keyword)
  {
    if (!op.is_empty ()) {
			std::ostringstream out_s;
			ar_operat_t opr (op.pop ());
			// loop while operator is not Open Bracket and stack not empty and priority greater equal <p>
			while (opr.priority >= priority && strcmp (opr.get_name ().c_str (), "("))
      {
        // call function
        if ((*opr.oper)(res,ar))
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": cannot process keyword " << keyword;
            BS_ASSERT(false) (out_s.str());
//				throw bs_exception("Keyword handlers
//				class",out_s.str().c_str());
            return -34;
          }
        opr = op.pop ();
      }
			op.push (opr);
		}
    return YS_SUCCESS;
  }

  template <typename strategy_t>
  int idata <strategy_t>::do_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, const char *keyword)
  {
    std::ostringstream out_s;
    ar_operat_t opr;
		
		if (!op.is_empty ()) {
			opr = op.pop ();

			// loop while operator not Open Bracket or stack not empty
			while (strcmp (opr.get_name ().c_str (), "("))
			{
        // call function
        if ((*opr.oper) (res, ar))        // check error code
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": cannot process keyword " << keyword;
            BS_ASSERT(false) (out_s.str());
//				throw bs_exception("idata::do_calculate",out_s.str().c_str());
            return -34;
          }
        // try to pop operator or function from stack
				if (!op.is_empty ())
					opr = op.pop ();
				else
					break;
      }
		}
    return YS_SUCCESS;
  }

  template <typename strategy_t>
  int idata <strategy_t>::read_arg_func (const sp_reader_t &reader, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op,
                                          char **start_ptr, const char *keyword, int *arguments_count, int flag_brek)
  {
    std::ostringstream out_s;

    char *p_end = 0;              // temp variable
    char n_buf[CHAR_BUF_LEN];     // buffer to store function or argument name
    int i = 0;

    p_end = *start_ptr;

    // name of argument or function can consist from letters or numbers or symbol '_'
    while (isalnum ((int) (*p_end)) || *p_end == '_')     // read name of argument or function
      {
        n_buf[i] = *p_end;        // copy name to buffer
        ++i;
        ++p_end;
      }
    n_buf[i] = '\0';

    // if nothing to read return 0
    if (strlen (n_buf) < 1)
      return 0;

    // try to find token in arguments list
    typename std::map <std::string,ar_args_t>::iterator sc = args.find (n_buf);
		if (sc != args.end ())
      {
        // if token is argument set token code
        ar_tokens[CURRENT_TOKEN] = TOKEN_ARG;
        // push argument to stack
        ar.push (sc->second);
      }
    else
      {
        // if token not found in argument list try to find it in function list
        typename std::set <ar_operat_t>::iterator opr = ops.find (ar_operat_t(n_buf));
        if (opr == ops.end ())
          {
            out_s << "Error in " << reader->get_prefix ()
            << ": unknown argument " << n_buf
            << " for keyword " << keyword << " is found ";
            BS_ASSERT(false) (out_s.str());
            throw bs_exception("idata::read_arg_func",out_s.str().c_str());
          }
        // if found set token code
        ar_tokens[CURRENT_TOKEN] = TOKEN_FUNC;

        // set number of function arguments to the future test
        if (flag_brek + 1 < AR_STACK_DEPTH)
          arguments_count[flag_brek + 1] = opr->num_of_arg - 1;

        // push token to the stack
        op.push (*opr);
      }
    // move pointer to the next token
    *start_ptr = p_end;
		return 1;
  }

#define DEFINE_METHODS(T) \
  template void idata <T>::build_argument_list (const smart_ptr <FRead, true>&); \
  template void idata <T>::build_operator_list (); \
  template void idata <T>::algorithm_read_and_done (const smart_ptr <FRead, true>&, char*, const char*); \
  template void idata <T>::clear_argument_list (); \
  template void idata <T>::output_argument_list (); \
  template int idata <T>::no_blanks (char *buf) const; \
  template void idata <T>::read_left (const smart_ptr <FRead, true>&,char*,char **right,std::string &name,int*,int*,int*,int*,int*,int*,const char*); \
  template void idata <T>::calculate (const sp_reader_t&, ar_args_t&, char*, const char*); \
  template int idata <T>::read_numbers (ar_stack <ar_args_t>&, char**); \
  template int idata <T>::read_operator (const sp_reader_t&, ar_args_t&, ar_stack <ar_args_t>&, ar_stack <ar_operat_t>&, char**, int*, const char*, const char*, int*); \
  template int idata <T>::do_calculate (const sp_reader_t&, ar_args_t&, ar_stack <ar_args_t>&, ar_stack <ar_operat_t>&, const char*); \
  template int idata <T>::prior_calculate (const sp_reader_t&, ar_args_t&, ar_stack <ar_args_t>&, ar_stack <ar_operat_t>&, int, const char*); \
  template int idata <T>::read_arg_func (const sp_reader_t&, ar_stack <ar_args_t>&, ar_stack <ar_operat_t>&, char**, const char*, int*, int); \
  template int idata <T>::test_token (int prev, int cur);

DEFINE_METHODS(base_strategy_fif);
DEFINE_METHODS(base_strategy_did);
DEFINE_METHODS(base_strategy_dif);

}
