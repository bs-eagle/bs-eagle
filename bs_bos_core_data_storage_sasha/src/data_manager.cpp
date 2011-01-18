#include "bs_bos_core_data_storage_stdafx.h"

#include "data_manager.h"
#include "main_def.h"
#include "arrays.h"
#include "localization.h"
//#include "arrays_tables.h"
#include "data_class.h"

/////////////////////////////////////////////////////////////////////////////
// You can add new variable just find next Words
// for symple variable                  -- VAR_V
// for array                                    -- ARRAY_V
// for strings                                  -- STRING_V

//! Default minimal DX value allowed for active cells
#define DEFAULT_MINIMAL_DX              1.e-6
//! Default minimal DY value allowed for active cells
#define DEFAULT_MINIMAL_DY              1.e-6
//! Default minimal DZ value allowed for active cells
#define DEFAULT_MINIMAL_DZ              1.e-10
//! Default minimal PERM value allowed for active cells
#define DEFAULT_MINIMAL_PERM            1e-16
#define DEFAULT_MINIMAL_NTG             1.e-8
#define DEFAULT_MINIMAL_TOPS            1.e-8
//! default blank value
#define DEFAULT_BLANK_VALUE             1.e43
//! check blank value
#define CHECK_BLANK_VALUE(X)  (fabs ((X) - DEFAULT_BLANK_VALUE) < 10)

//! Default minimal PORO value allowed for active cells
#define DEFAULT_MINIMAL_PORO                        1.e-10
#define DEFAULT_MINIMAL_PORO_VOLUME_MULT            1.e-12
/////////////////////////////////////////////////////////////////////////////


// Please don't delete commented old style outputs .. like rep->blah blah blah
namespace blue_sky
  {
#define DM_ASSERT_EXCEPTION \
    BS_ASSERT(false) (out_s.str());\
    throw bs_exception("Data manager class",out_s.str().c_str());

  template <class strategy_t>
  data_manager<strategy_t>::~data_manager ()
  {

  }

  template <class strategy_t>
  data_manager<strategy_t>::data_manager(bs_type_ctor_param /*param*/): lkeeper ("C", LC_ALL)
  {
    this->reader = BS_KERNEL.create_object(FRead::bs_type());
    //this->key_handlers = BS_KERNEL.create_object(keyword_manager <strategy_t>::bs_type());
    this->data = BS_KERNEL.create_object(idata <strategy_t>::bs_type());
  }

  template <class strategy_t>
  data_manager<strategy_t>::data_manager(const data_manager& src)
  : bs_refcounter (src), objbase (src), lkeeper ("C", LC_ALL)
  {
    *this = src;
  }

  //! macros for checking
#define CHECK_FOR_TWO_PARAMS(VAL,ACTNUM,MIN_DEF,COUNTER)                                \
	if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL))                                      \
{                                                                                   \
	if ((ACTNUM))                                                                     \
	{                                                                               \
		int ix, iy, iz;                                                               \
		iz = i / (data->nx * data->ny);                                               \
		iy = (i - iz * (data->nx * data->ny)) / data->nx;                             \
		ix = i - iz * (data->nx * data->ny) - iy * data->nx;                          \
		BOSWARN (section::check_data, level::low)                                     \
		<< "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
		<< "] will be set inactive because of " << #VAL << " = " << VAL << bs_end;  \
		(ACTNUM) = 0;                                                                 \
		++(COUNTER);                                                                  \
		continue;                                                                     \
	}                                                                               \
}

#define CHECK_FOR_THREE_PARAMS(VAL,ACTNUM,MIN_DEF,MAX_DEF,COUNTER)                      \
	if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL) || (VAL) > (MAX_DEF))                 \
{                                                                                   \
	if ((ACTNUM))                                                                     \
	{                                                                               \
		int ix, iy, iz;                                                               \
		iz = i / (nx * ny);                                                           \
		iy = (i - iz * (nx * ny)) / nx;                                               \
		ix = i - iz * (nx * ny) - iy * nx;                                            \
		BOSWARN (section::check_data, level::low)                                     \
		<< "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
		<< "] will be set inactive because of " << #VAL << bs_end;                  \
		(ACTNUM) = 0;                                                                 \
		++(COUNTER);                                                                  \
		continue;                                                                     \
	}                                                                               \
}

#define CHECK_DIRECTION_FOR_XYZ(case,num)   \
	if ((case) && actnum [(num)])             \
{                                       \
	n++;                                  \
	valx += dx[(num)];                    \
	valy += dy[(num)];                    \
	valz += dz[(num)];                    \
}

#define CHECK_DIRECTION_FOR_TOPS(case,num)  \
	if ((case) && actnum[(num)])              \
{                                       \
	n++;                                  \
	valtops += tops[(num)];               \
}

  /*
namespace blue_sky {
using namespace pool;

#define DM_ASSERT_EXCEPTION \
BS_ASSERT(false) (out_s.str());\
throw bs_exception("Data manager class",out_s.str().c_str());

template <class strategy_t>
data_manager<strategy_t>::~data_manager ()
{

}

template <class strategy_t>
data_manager<strategy_t>::data_manager(bs_type_ctor_param param): lkeeper ("C", LC_ALL)
{
	this->reader = BS_KERNEL.create_object(FRead::bs_type());
	//this->key_handlers = BS_KERNEL.create_object(keyword_handlers <strategy_t>::bs_type());
	this->data = BS_KERNEL.create_object(idata <strategy_t>::bs_type());
}

template <class strategy_t>
data_manager<strategy_t>::data_manager(const data_manager& src): lkeeper ("C", LC_ALL)
{
	*this = src;
}

template <class strategy_t>
void data_manager<strategy_t>::update_geometry() const
{
	int i, ix, iy, iz, n;
	std::ostringstream out_s;
	int nb = data->nx * data->ny * data->nz;

	int n_tops, n_xyz;
	double dx_var, dy_var, dz_var, tops_var;
	double valx, valy, valz, valtops;

	bs_array_fp& dx = data->get_non_empty_fp("DX");
	bs_array_fp& dy = data->get_non_empty_fp("DY");
	bs_array_fp& dz = data->get_non_empty_fp("DZ");
	bs_array_fp& tops = data->get_non_empty_fp("TOPS");
	bs_array_i& actnum = data->get_non_empty_i("ACTNUM");

	// In case of restart we got all geometry from mesh restored from
	// the restart file. So we skip this test.
	if (data->restart)
		return;

	if (!nb)
	{
		bs_throw_exception ("update_geometry: all dimensions must be greater than 0");
	}

	dx_var = dy_var = dz_var = tops_var = 0.;
	valx = valy = valz = valtops = 0.;
	n_tops = n_xyz = 0;
	if (data->geom_def_flag == GEOM_FLAG_DX_DY_DZ)
	{
		for (iz = 0, i = 0; iz < data->nz; ++iz)
		{
			for (iy = 0; iy < data->ny; ++iy)
			{
				for (ix = 0; ix < data->nx; ++ix, ++i)
				{
					if (actnum[i])
					{
						dx_var += dx[i];
						dy_var += dy[i];
						dz_var += dz[i];
						n_xyz++;
					}
				}
			}
		}

		if (!n_xyz)
		{
			bs_throw_exception ("update_geometry: No active blocks");
		}

		for (iy = 0, i = 0; iy < data->ny; ++iy)
		{
			for (ix = 0; ix < data->nx; ++ix, ++i)
			{
				if (actnum[i])
				{
					tops_var += tops[i];
					++n_tops;
				}
			}
		}
		if (!n_tops)
		{
			bs_throw_exception ("update_geometry: No active blocks");
		}

		dx_var /= n_xyz;
		dy_var /= n_xyz;
		dz_var /= n_xyz;
		tops_var /= n_tops;
		for (iy = 0, i = 0; iy < data->ny; ++iy)
		{
			for (ix = 0; ix < data->nx; ++ix, ++i)
			{
				if (!actnum[i])
				{
					valtops = 0.;
					n = 0;
					CHECK_DIRECTION_FOR_TOPS (ix > 0, i - 1);
					CHECK_DIRECTION_FOR_TOPS (iy > 0, i - data->nx);
					CHECK_DIRECTION_FOR_TOPS (ix < data->nx - 1 , i + 1);
					CHECK_DIRECTION_FOR_TOPS (iy < data->ny - 1 , i + data->nx);
					if (n)
					{
						tops[i] = valtops / n;
					}
					else
					{
						tops[i] = tops_var;

						BOSWARN (section::check_data, level::low) 
							<< "new tops in node (" << (ix + 1)
							<< ", " << (iy + 1) << ", " << (iz + 1) << ") = " << tops[i]
							<< " (average = " << tops_var << ")" << bs_end;
					}
				}
			}
		}

		for (iz = 0, i = 0; iz < data->nz; ++iz)
		{
			for (iy = 0; iy < data->ny; ++iy)
			{
				for (ix = 0; ix < data->nx; ++ix, ++i)
				{
					if (!actnum[i])
					{
						valx = valy = valz = 0.;
						n = 0;
						CHECK_DIRECTION_FOR_XYZ (ix > 0, i - 1);
						CHECK_DIRECTION_FOR_XYZ (iy > 0, i - data->nx);
						CHECK_DIRECTION_FOR_XYZ (iz > 0, i - data->nx * data->ny );
						CHECK_DIRECTION_FOR_XYZ (ix < data->nx - 1 , i + 1);
						CHECK_DIRECTION_FOR_XYZ (iy < data->ny - 1 , i + data->nx);
						CHECK_DIRECTION_FOR_XYZ (iz < data->nz - 1 , i + data->nx * data->ny);
						if (n)
						{
							dx[i] = valx / n;
							dy[i] = valy / n;
							dz[i] = valz / n;
						}
						else
						{
							dx[i] = dx_var;
							dy[i] = dy_var;
							dz[i] = dz_var;
						}

						BOSOUT (section::check_data, level::low)
							<< "Info: new dx in node (" << (ix + 1) << ", "
							<< (iy + 1) << ", " << (iz + 1) << ") = " << dx[i] << " (aver dx = "
							<< dx_var << ")" 
                              << "\n"
							<< "Info: new dy in node (" << (ix + 1) << ", "
							<< (iy + 1) << ", " << (iz + 1) << ") = " << dy[i] << " (aver dy = "
							<< dy_var << ")" 
                              << "\n"
							<< "Info: new dz in node (" << (ix + 1) << ", "
							<< (iy + 1) << ", " << (iz + 1) << ") = " << dz[i] << " (aver dz = "
							<< dz_var << ")" 
							<< bs_end;
					}
				}
			}
		}
	}
	else if (data->geom_def_flag == GEOM_FLAG_ZCORN_COORD)
	{
#if 0
		for (iz = 0, i = 0; iz <= nz; ++iz)
		{
			for (iy = 0; iy <= ny; ++iy)
			{
				for (ix = 0; ix <= nx; ++ix, ++i)
				{
					if ((*depth)[i] < 0)
					{
						BOSERR (section::check_data, level::error) << "Depth of node (" << (ix + 1) << ", "
							<< (iy + 1) << ", " << (iz + 1) << ") = " << depth[i]
							<< " is out of range" << bs_end;
						return YS_BAD_VALUE;
					}
				}
			}
		}
#endif
	}
	else
	{
		bs_throw_exception ("Unsupported geometry specification");
	}
}
*/
template <class strategy_t>
void data_manager<strategy_t>::check_arrays_for_inactive_blocks () const
{
	std::ostringstream out_s;
	int nb = data->nx * data->ny * data->nz;
	int permx_counter = 0;
	int permy_counter = 0;
	int permz_counter = 0;
	int poro_counter = 0;
	int ntg_counter = 0;

	bs_array_i& actnum = data->get_non_empty_i("ACTNUM");
	bs_array_fp& permx  = data->get_non_empty_fp("PERMX");
	bs_array_fp& permy  = data->get_non_empty_fp("PERMY");
	bs_array_fp& permz  = data->get_non_empty_fp("PERMZ");
	bs_array_fp& poro   = data->get_non_empty_fp("PORO");
	bs_array_fp& ntg    = data->get_non_empty_fp("NTG");

	for (int i = 0; i < nb; ++i)
	{
		if (!actnum[i])
		{
			continue;
		}

		if (actnum[i] && permx[i] < DEFAULT_MINIMAL_PERM && permy[i] < DEFAULT_MINIMAL_PERM && permz[i] < DEFAULT_MINIMAL_PERM)
		{
			int ix, iy, iz;
			iz = i / (data->nx * data->ny);
			iy = (i - iz * (data->nx * data->ny)) / data->nx;
			ix = i - iz * (data->nx * data->ny) - iy * data->nx;

			BOSWARN (section::check_data, level::low)
				<< "blocks " << i << " [" << (ix + 1) << ", " << (iy + 1) << ", "
				<< (iz + 1) << "] will be set inactive because of PERMEABILITY" << bs_end;

			actnum[i] = 0;
			++permx_counter;
			continue;
		}

		CHECK_FOR_TWO_PARAMS (poro[i], actnum[i], DEFAULT_MINIMAL_PORO, poro_counter);
          if (!ntg.empty ())
            {
              CHECK_FOR_TWO_PARAMS (ntg[i], actnum[i], DEFAULT_MINIMAL_NTG, ntg_counter);
            }
	}

	if (permx_counter)
		BOSWARN (section::check_data, level::warning) << permx_counter << " blocks will be set inactive because of PERMX" << bs_end;
	if (permy_counter)
		BOSWARN (section::check_data, level::warning) << permy_counter << " blocks will be set inactive because of PERMY" << bs_end;
	if (permz_counter)
		BOSWARN (section::check_data, level::warning) << permz_counter << " blocks will be set inactive because of PERMZ" << bs_end;
	if (poro_counter)
		BOSWARN (section::check_data, level::warning) << poro_counter << " blocks will be set inactive because of PORO" << bs_end;
	if (ntg_counter)
		BOSWARN (section::check_data, level::warning) << ntg_counter << " blocks will be set inactive because of NTG" << bs_end;
}

//bs stuff
	BLUE_SKY_TYPE_STD_CREATE_T_DEF(data_manager, (class))
BLUE_SKY_TYPE_STD_COPY_T_DEF(data_manager, (class))
	BLUE_SKY_TYPE_IMPL_T_SHORT(data_manager<base_strategy_did>, objbase, "BOS_Core data_manager class")
	BLUE_SKY_TYPE_IMPL_T_SHORT(data_manager<base_strategy_fif>, objbase, "BOS_Core data_manager class")
	BLUE_SKY_TYPE_IMPL_T_SHORT(data_manager<base_strategy_dif>, objbase, "BOS_Core data_manager class")

}//ns bs
