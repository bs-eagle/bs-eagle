#include "bs_bos_core_data_storage_stdafx.h"
#include "ar_args.h"

namespace blue_sky
  {
  template <typename strategy_t>
  ar_args<strategy_t>::ar_args(const std::string &tname,
                   itype itype_,
                   void *tarray,
//					const sp_item_array_t &tarray,
//					lva_t &tlarray,
                   const int tni,        //!< first dimension size
                   const int tnj,        //!< second dimension size
                   const int tnk,        //!< third dimension size
                   const int ti1,    //!< begin chek array dimmention 1
                   const int ti2,    //!< end of chek array dimmention 1
                   const int tj1,    //!< begin of chek array dimmention 2
                   const int tj2,    //!< end of chek array dimmention 2
                   const int tk1,    //!< begin of chek array dimmention 3
                   const int tk2     //!< end of chek array dimmention 3.
                  )
  {
    arr_type = itype_;
    if (tname.length() == 0)
      throw bs_exception("ar_args ctor","Can't create array with 'empty string'-like name");
    name = tname;
    std::transform(name.begin(),name.end(),name.begin(),toupper);

    //int arr_size = tarray->size();
    //def_flag = (arr_size ? 1 : 0);
    def_flag = (tarray ? 1 : 0);

    if (!(n_size = tni * tnj * tnk))
      throw bs_exception("ar_args ctor","Error: you should specify keyword 'DIMENS' before keyword 'ARITHMETIC'");

    //BS_ASSERT(n_size == arr_size);

    array = tarray;
    //array.resize (0);
    //array.assign (tarray.begin (), tarray.end ());
    flag = 0;
    ni = tni;
    nj = tnj;
    nk = tnk;
    i1 = ti1;
    i2 = ti2;
    j1 = tj1;
    j2 = tj2;
    k1 = tk1;
    k2 = tk2;

//			larray = tlarray;

  }

	template <typename strategy_t>
  ar_args<strategy_t>::ar_args(const ar_args &a)
  {
    ni = a.ni;
    nj = a.nj;
    nk = a.nk;
    i1 = a.i1;
    i2 = a.i2;
    j1 = a.j1;
    j2 = a.j2;
    k1 = a.k1;
    k2 = a.k2;
    flag = 1;
    name = a.name;
    def_flag = 0;
    n_size = ni * nj * nk;
    arr_type = a.arr_type;

    //int arr_size = a.array->size();

    //BSOUT << n_size << ", " << array.size()
    //BS_ASSERT(n_size == arr_size);

    if (a.array && n_size)
      {
        //array.resize (0);
        //array.assign (a.array.begin (), a.array.end ());
        //array = a.array;
				if (arr_type == FPOINT_T) {
					array = (void*)(new item_t[n_size.data ()]);
					memcpy (array, a.array, n_size * sizeof (item_t));
				}
				else if (arr_type == INT_T) {
					array = (void*)(new int[n_size.data ()]);
					memcpy (array, a.array, n_size * sizeof (int));
				}
      }
  }

	template <typename strategy_t>
  ar_args<strategy_t>::ar_args(itype itype_, void *tarray, /*const sp_item_array_t &tarray,*/ const std::string &tname)
  {
    i1 = 0;
    i2 = 0;
    j1 = 0;
    j2 = 0;
    k1 = 0;
    k2 = 0;
    ni = 0;
    nj = 0;
    nk = 0;
    name = "";
    flag = 0;
    def_flag = 0;
    arr_type = itype_;

    if (tarray)
    {
        n_size = 1;
        if (arr_type == FPOINT_T)
					array = (void*)(new item_t[1]);
				else if (arr_type == INT_T)
					array = (void*)(new int[1]);
				else {
					BOSERR (section::arithmetic, level::error) << "unknown array type in arithmetic" << bs_end;
				}
				if (array == 0) {
					BOSERR (section::arithmetic, level::error) << "no memory allocated for array in arithmetic." << bs_end;
				}
				else {
					if (arr_type == FPOINT_T)
						((item_t*)array)[0] = *((item_t*)tarray);
					else
						((int*)array)[0] = *((int*)tarray);
					ni = nj = nk = 1;
					flag = 1;
					if (tname.length()) {
						name = tname;
						std::transform(name.begin(),name.end(),name.begin(),toupper);
					}
				}
		}
  }

	template <typename strategy_t>
	ar_args<strategy_t>::ar_args(const std::string &tname) {
		name = tname;
		std::transform(name.begin(),name.end(),name.begin(),toupper);
	}

	template <typename strategy_t>
  ar_args<strategy_t>::~ar_args()
  {
    i1 = 0;
    i2 = 0;
    j1 = 0;
    j2 = 0;
    k1 = 0;
    k2 = 0;
    ni = 0;
    nj = 0;
    nk = 0;
    n_size = 0;
    def_flag = 0;
    name = "";
    if (flag == 1)
      {
        if (arr_type == FPOINT_T)
          delete [] (item_t*)array;
        else if (arr_type == INT_T)
          delete [] (int*)array;
        //else
          //delete [] array;
        //array->clear();
      }
    flag = 0;
  }

	template <typename strategy_t>
  int ar_args<strategy_t>::set_ptr_data()
  {
//		int arr_size = 0; //array->size ();

    //if (!arr_size)
    if (!array)
      throw bs_exception("ar_args ctor","Array is empty, but set_ptr_data try to use it");
    flag = 0;
//		larray.clear();
//		larray.push_back(array);
    return YS_SUCCESS;
  }

	template <typename strategy_t>
  const std::string &ar_args<strategy_t>::get_name() const
    {
      //typename ar_args<strategy_t>::index_t
      //ar_args<strategy_t>::get_name () const {
      return name;
    }

	template <typename strategy_t>
  /*const typename ar_args<strategy_t>::sp_item_array_t &*/
  void *ar_args<strategy_t>::get_array()
  {
    return array;
  }

//	template <typename strategy_t>
//	typename ar_args<strategy_t>::item_t &ar_args<strategy_t>::get (unsigned i) {
//		return array->operator [](i);
//	}

	template <typename strategy_t>
  ar_args<strategy_t> &ar_args<strategy_t>::operator=(const ar_args &a)
  {
    if (this == &a)
      return *this;

		i1 = 0;
		i2 = 0;
		j1 = 0;
		j2 = 0;
		k1 = 0;
		k2 = 0;
		ni = 0;
		nj = 0;
		nk = 0;
		name = "";
		flag = 0;
		def_flag = 0;
		arr_type = a.arr_type;

		if (a.array)
		{
			n_size = 1;
			if (arr_type == FPOINT_T)
				array = (void*)(new item_t[1]);
			else if (arr_type == INT_T)
				array = (void*)(new int[1]);
			else {
				BOSERR (section::arithmetic, level::error) << "unknown array type in arithmetic" << bs_end;
			}
			if (array == 0) {
				BOSERR (section::arithmetic, level::error) << "no memory allocated for array in arithmetic." << bs_end;
			}
			else {
				if (arr_type == FPOINT_T)
					((item_t*)array)[0] = *((item_t*)a.array);
				else
					((int*)array)[0] = *((int*)a.array);
				ni = nj = nk = 1;
				flag = 1;
				if (a.name.length()) {
					name = a.name;
					std::transform(name.begin(),name.end(),name.begin(),toupper);			
				}
			}
		}

    return *this;
  }

	template <typename strategy_t>
  typename ar_args<strategy_t>::itype ar_args<strategy_t>::get_type () const
    {
      return arr_type;
    }

	template <typename strategy_t>
  int ar_args<strategy_t>::operator==(const std::string &tname) const
    {
      return strcmp(name.c_str(),tname.c_str());
    }

	//bool ar_args::operator==(const std::string &tname) const {
	//	return (!strcmp(name.c_str(),tname.c_str())) ? true : false;
	//}

	template <typename strategy_t>
  bool ar_args<strategy_t>::operator<(const ar_args &arg) const
    {
      return (name < arg.name);
    }

	template <typename strategy_t>
  int ar_args<strategy_t>::allocate()
  {
    //int arr_size = array->size ();

    if (array)
      {
        BOSERR (section::arithmetic, level::error) << "Internal error. Memory already allocate!" << bs_end;
        return -1;
      }
    int nb = ni * nj * nk;
    if (!nb)
      {
        BOSERR (section::arithmetic, level::error) << "Error: you should specify keyword 'DIMENS' before keyword 'ARITHMETIC'" << bs_end;
        return -2;
      }
    n_size = nb;
    if (arr_type == FPOINT_T)
      {
        array = (void *)(new item_t[nb]);
        memset (array, 0, sizeof (item_t) * nb);
      }
    else if (arr_type == INT_T)
      {
        array = (void *)(new int[nb]);
        memset (array, 0, sizeof (int) * nb);
      }
    else
      {
        BOSERR (section::arithmetic, level::error) << "Unknown array type in ARITHMETIC" << bs_end;
        return -1;
      }

    if (!array)
      {
        BOSERR (section::arithmetic, level::error) << "Not enought memory" << bs_end;
        return -3;
      }

    //array->resize (n_size);

    def_flag = 1;
    flag = 1;
    return YS_SUCCESS;
  }

//	BLUE_SKY_TYPE_IMPL_T_EXT(1, (mesh_grdecl_mpfa<base_strategy_fif>) , 1 , (mesh_grdecl<base_strategy_fif>), "mesh_grdecl_mpfa<float, int,float>", "GRD_ECL Mesh class for mpfa calculating", "Grid Ecllipse Mesh class  for mpfa calculating", false);
//	BLUE_SKY_TYPE_IMPL_T_EXT(1, (mesh_grdecl_mpfa<base_strategy_did>) , 1 , (mesh_grdecl<base_strategy_did>), "mesh_grdecl_mpfa<double, int,double>", "GRD_ECL Mesh class for mpfa calculating", "Grid Ecllipse Mesh class for mpfa calculating", false);

//	BLUE_SKY_TYPE_STD_CREATE_T_DEF(mesh_grdecl_mpfa, (class));
//	BLUE_SKY_TYPE_STD_COPY_T_DEF(mesh_grdecl_mpfa, (class));

  template class ar_args <base_strategy_fif>;
  template class ar_args <base_strategy_did>;
  template class ar_args <base_strategy_dif>;
}
