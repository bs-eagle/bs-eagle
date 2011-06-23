/*! \file ar_args.h
    \brief Contains a class for store arguments in arithmetic operation.
*/
#ifndef ARITHMETIC_ARGS_H
#define ARITHMETIC_ARGS_H


#define DEF_SIZE 100

namespace blue_sky
  {
  // amap_strategy_xi using here
  class BS_API_PLUGIN ar_args
    {
    public:
      // typedefs
      typedef strategy_t::item_t         item_t;
      typedef strategy_t::index_t        index_t;
      //typedef typename strategy_t::item_array_t   item_array_old_t;
      //typedef bos_array <item_array_old_t>        item_array_t;
      //typedef smart_ptr <item_array_t>            sp_item_array_t;
      //typedef typename strategy_t::index_array_t  index_array_t;

      //typedef ar_args <strategy_t> this_t;

      typedef enum
      {
        INT_T,
        FPOINT_T
      } itype;

      //! ctors
      ar_args(const std::string &tname,
              //const sp_item_array_t &tarray,
              itype itype_,
              void *tarray,
//						lva_t &tlarray,
              const int tni,        //!< first dimension size
              const int tnj,        //!< second dimension size
              const int tnk,        //!< third dimension size
              const int ti1 = 0,    //!< begin chek array dimmention 1
              const int ti2 = 0,    //!< end of chek array dimmention 1
              const int tj1 = 0,    //!< begin of chek array dimmention 2
              const int tj2 = 0,    //!< end of chek array dimmention 2
              const int tk1 = 0,    //!< begin of chek array dimmention 3
              const int tk2 = 0     //!< end of chek array dimmention 3.
             );

      ar_args(const ar_args &a);
      ar_args(itype itype_, void *tarray/*const sp_item_array_t &tarray*/, const std::string &tname = "");

			ar_args(const std::string &tname);

      ~ar_args();

      int set_ptr_data();

      const std::string &get_name() const;

      //const sp_item_array_t &get_array();
      void *get_array ();
      //item_t &get (unsigned i);

      itype get_type () const;

      //! Assignment operator
      ar_args &operator=(const ar_args &a);

      //! Compare operator
      int operator==(const std::string &tname) const;
			//bool operator==(const std::string &tname) const;
      bool operator<(const ar_args &arg) const;

      /*!
       \brief function -- if memory have not been allocate before, allocate memory for array
       \return        if success                      YS_SUCCESS
       \return        if nb == 0                      -2
       \return        if memory already allocate      -3
       \return        if cannot allocate memory       -4
      */
      int allocate();

      ///////////////////////////////////
      ///////////////Data////////////////
      ///////////////////////////////////

      auto_value<int> ni;                       //!< first dimension size
      auto_value<int> nj;                       //!< second dimension size
      auto_value<int> nk;                       //!< third dimension size
      auto_value<int> i1, i2;                   //!< Calculate from i1 to i2
      auto_value<int> j1, j2;                   //!< Calculate from j1 to j2
      auto_value<int> k1, k2;                   //!< Calculate from k1 to k2
      auto_value<int> def_flag;                 //!< if def_flag = 1 then array undefined
      auto_value<int> n_size;

    protected:
      itype arr_type;
      void *array;
//		sp_item_array_t array;       //!< array of data
//		lva_t larray;     //!< list of previous arrays
      std::string name; //!< name of array
      auto_value<int> flag;         //!< if flag == 1 destructor must call free for array
    };
}

#endif // ARITHMETIC_ARGS_H
