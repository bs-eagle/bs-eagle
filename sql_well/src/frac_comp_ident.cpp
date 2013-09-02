/// @file frac_comp_ident.cpp
/// @brief Implementation of finding fractions/completions intersections with mesh
/// @author uentity
/// @version
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "frac_comp_ident.h"
#include "bcsr_matrix_iface.h"
#include "wpi_algo.h"
#include "wpi_strategy_3d.h"
#include "export_python_wrapper.h"
#include "py_list_converter.h"
#include "sql_well.h"

#include "compdat_traits.h"
#include "fracture_traits.h"

#include <boost/format.hpp>
#include <iosfwd>
#include <iostream>

namespace blue_sky { namespace fci {
// impl details
using namespace std;

typedef wpi::algo< wpi::strategy_3d > wpi_algo;

// hidden details
namespace {

// simple dump
void dump(std::ostream& os, const compdat& cd) {
	os << setw(15) << cd.well_name << ' ' << setw(15) << cd.branch_name << ' ';
	os << setw(10) << std::fixed << std::setprecision(3) << cd.md;
	os << setw(10) << std::fixed << std::setprecision(3) << cd.len;
	for(int i = 0; i < 4; ++i)
		os << setw(4) << cd.cell_pos[i];
	os << setw(2) << cd.dir;
	os << setw(10) << std::fixed << std::setprecision(3) << cd.kh_mult << std::endl;
}

// simple fracture dump
void dump(std::ostream& os, const fracture& frac) {
	os << setw(15) << frac.well_name << ' ' << setw(15) << frac.branch_name << ' ';
	os << setw(3) << frac.frac_status << ' ';
	for(int i = 0; i < 4; ++i)
		os << setw(4) << frac.cell_pos[i];
	os << setw(10) << std::fixed << std::setprecision(3) << frac.frac_half_length_1;
	os << setw(10) << std::fixed << std::setprecision(3) << frac.frac_half_length_2;
	os << setw(10) << std::fixed << std::setprecision(3) << frac.frac_angle;
	os << setw(10) << std::fixed << std::setprecision(3) << frac.frac_half_thin;
	os << setw(10) << std::fixed << std::setprecision(3) << frac.frac_perm << std::endl;
}

template <class x_storage>
void dump(std::ostream& os, const x_storage& cfs) {
	typename x_storage::const_iterator pc;
	typename x_storage::const_iterator end;
	for(pc = cfs.begin(), end = cfs.end(); pc != end; ++pc)
		dump(os, *pc);
}

template< class brick >
struct choose_traits {
	typedef compl_traits type;
};
template< >
struct choose_traits< fracture > {
	typedef fract_traits type;
};


}  // oef hidden namespace

/*-----------------------------------------------------------------
 * compdat
 *----------------------------------------------------------------*/
compdat::compdat(const string& well_name_, const string& branch_name_, const pos_i& cell_pos_, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), dir(' '), kh_mult(0), md(0), len(0), skin (0), diam (0.146), kh (0), status (0)
{
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	cell_id_ = wpi_algo::encode_cell_id(cell_pos_, mesh_size);

  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

compdat::compdat(const string& well_name_, const string& branch_name_, ulong cell_id, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), dir(' '), kh_mult(0), md(0), len(0), skin (0), diam (0.146), kh (0), status (0)
{
	init(cell_id, mesh_size);
}

compdat::compdat(ulong cell_id)
	: dir(' '), kh_mult(0), md(0), len(0), skin (0), diam (0.146), kh (0), status (0), cell_id_(cell_id)
{
  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

void compdat::init(ulong cell_id, const pos_i& mesh_size) {
	// convert cell id -> cell pos
	pos_i cell_pos_;
	wpi_algo::decode_cell_id(cell_id, cell_pos_, mesh_size);
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	// store cell_id
	cell_id_ = cell_id;
  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

/*-----------------------------------------------------------------
 * fracture
 *----------------------------------------------------------------*/
fracture::fracture(const string& well_name_, const string& branch_name_, const pos_i& cell_pos_, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), frac_status(0), frac_half_length_1 (0),
	  frac_half_length_2 (0), frac_angle (0), frac_half_thin (0), frac_perm (0), frac_main_k (-1)
{
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &md_cell_pos[0]);
	cell_id_ = wpi_algo::encode_cell_id(cell_pos_, mesh_size);
}

fracture::fracture(const string& well_name_, const string& branch_name_, ulong cell_id, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_),frac_status(0), frac_half_length_1 (0),
	  frac_half_length_2 (0),  frac_angle (0), frac_half_thin (0), frac_perm (0), frac_main_k (-1)
{
	init(cell_id, mesh_size);
}

fracture::fracture(ulong cell_id)
	: frac_status(0), frac_half_length_1 (0), frac_half_length_2 (0),
	  frac_angle (0), frac_half_thin (0), frac_perm (0), frac_main_k (-1), cell_id_(cell_id)
{
}

void fracture::init(ulong cell_id, const pos_i& mesh_size) {
	// convert cell id -> cell pos
	pos_i cell_pos_;
	wpi_algo::decode_cell_id(cell_id, cell_pos_, mesh_size);
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &md_cell_pos[0]);
	// store cell_id
	cell_id_ = cell_id;
}

/*-----------------------------------------------------------------
 * builder implementation for any brick type
 *----------------------------------------------------------------*/
template< class brick >
class builder< brick >::impl : public frac_comp_builder< typename choose_traits< brick >::type > {};

template< class brick >
builder< brick >::builder() : pimpl_(new impl) {}

template< class brick >
void builder< brick >::init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn);
}

template< class brick >
void builder< brick >::init(t_ulong nx, t_ulong ny, sp_obj trim_backend) {
	pimpl_->init(nx, ny, trim_backend);
}

template< class brick >
void builder< brick >::init(smart_ptr< well_pool_iface, true > src_well) {
	pimpl_->init(src_well);
}

template< class brick >
const typename builder< brick >::storage_t& builder< brick >::build(double date) {
	pimpl_->build(date);
	return storage();
}

template< class brick >
const typename builder< brick >::storage_t& builder< brick >::storage() const {
	return pimpl_->s_;
}

template< class brick >
void builder< brick >::clear() {
	pimpl_->s_.clear();
}

template< class brick >
template< class B >
void builder< brick >::share_cache_with(const builder< B >& rhs) {
	pimpl_->init_cache(rhs->pimpl_->xp_cache_, rhs->pimpl_->cache_limit_);
}

/*-----------------------------------------------------------------
 * compdat_builder implementation
 *----------------------------------------------------------------*/
compdat_builder::compdat_builder() {}

compdat_builder::compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
	builder< compdat >::init(nx, ny, coord, zcorn);
}

compdat_builder::compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
{
	builder< compdat >::init(nx, ny, coord, zcorn);
	builder< compdat >::init(src_well);
}

/*-----------------------------------------------------------------
 * fracture_builder implementation
 *----------------------------------------------------------------*/
fracture_builder::fracture_builder() {}

fracture_builder::fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	builder< fracture >::init(nx, ny, coord, zcorn);
}

fracture_builder::fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
{
	builder< fracture >::init(nx, ny, coord, zcorn);
	builder< fracture >::init(src_well);
}

/*-----------------------------------------------------------------
 * fracture_builder::impl
 *----------------------------------------------------------------*/
class compl_n_frac_builder::impl {
public:
	typedef frac_comp_builder< compl_traits >::xp_cache_t xp_cache_t;
	typedef frac_comp_builder< compl_traits >::spxp_cache_t spxp_cache_t;
	typedef smart_ptr< well_pool_iface, true > sp_srcwell;

	frac_comp_builder< compl_traits > cb_;
	frac_comp_builder< fract_traits > fb_;
	spxp_cache_t xp_cache_;

	// ctors
	impl() {
		init_cache();
	}
	impl(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn) {
		init(nx, ny, coord, zcorn);
		init_cache();
	}
	impl(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn,
		const sp_srcwell& src_well)
	{
		init(nx, ny, coord, zcorn);
		init(src_well);
		init_cache();
	}

	void init(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn) {
		cb_.init(nx, ny, coord, zcorn);
		fb_.init(nx, ny, coord, zcorn);
	}

	void init(t_ulong nx, t_ulong ny, const sp_obj& trim_backend) {
		cb_.init(nx, ny, trim_backend);
		fb_.init(nx, ny, trim_backend);
	}

	void init(const sp_srcwell& src_well) {
		cb_.init(src_well);
		fb_.init(src_well);
	}

	void init_cache(const spxp_cache_t xc = NULL, const ulong cache_limit = 0) {
		if(xc)
			xp_cache_ = xc;
		else
			xp_cache_ = new xp_cache_t;
		cb_.init_cache(xp_cache_, cache_limit);
		fb_.init_cache(xp_cache_, cache_limit);
	}

	const cd_storage& compl_build(double date) {
		cb_.build(date);
		return cb_.s_;
	}

	const frac_storage& frac_build(double date) {
		fb_.build(date);
		return fb_.s_;
	}

	void clear() {
		cb_.s_.clear();
		fb_.s_.clear();
	}

	void share_cache_with(const compl_n_frac_builder& rhs) {
		init_cache(rhs.pimpl_->xp_cache_, rhs.pimpl_->cb_.cache_limit_);
	}
};


/*-----------------------------------------------------------------
 * comp_n_frac_builder implementation
 *----------------------------------------------------------------*/
compl_n_frac_builder::compl_n_frac_builder()
	: pimpl_(new impl())
{}

compl_n_frac_builder::compl_n_frac_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
	: pimpl_(new impl(nx, ny, coord, zcorn))
{}

compl_n_frac_builder::compl_n_frac_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
	: pimpl_(new impl(nx, ny, coord, zcorn, src_well))
{}

void compl_n_frac_builder::init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn);
}

void compl_n_frac_builder::init(t_ulong nx, t_ulong ny, sp_obj trim_backend) {
	pimpl_->init(nx, ny, trim_backend);
}

void compl_n_frac_builder::init(smart_ptr< well_pool_iface, true > src_well) {
	pimpl_->init(src_well);
}

const cd_storage& compl_n_frac_builder::compl_build(double date) {
	return pimpl_->compl_build(date);
}

const frac_storage& compl_n_frac_builder::frac_build(double date) {
	return pimpl_->frac_build(date);
}

const cd_storage& compl_n_frac_builder::storage_compdat ()  const {
	return pimpl_->cb_.s_;
}

const frac_storage& compl_n_frac_builder::storage_fracture () const {
	return pimpl_->fb_.s_;
}

void compl_n_frac_builder::clear() {
	pimpl_->clear();
}

void compl_n_frac_builder::share_cache_with(const compl_n_frac_builder& rhs) {
	pimpl_->share_cache_with(rhs);
}

/*-----------------------------------------------------------------
 * global functions impl
 *----------------------------------------------------------------*/
spv_float completions_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	compdat_builder b(nx, ny, coord, zcorn, src_well);
	b.build(date);
	dump(std::cout, b.storage());
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

spv_float fractures_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	fracture_builder b(nx, ny, coord, zcorn, src_well);
	b.build(date);
	dump(std::cout, b.storage());
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

} // eof blue_sky::fci

/*-----------------------------------------------------------------
 * Python export routines
 *----------------------------------------------------------------*/

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cdb_build_overl, blue_sky::fci::compdat_builder::build, 1, 2)

namespace python {
namespace bp = boost::python;

static void print_compdat(const fci::compdat& cd) {
	fci::dump(std::cout, cd);
}

static void print_fracture(const fci::fracture& frac) {
	fci::dump(std::cout, frac);
}


void py_export_compdat_ident() {
	using namespace fci;

	// export global functions
	bp::def("completions_ident", &fci::completions_ident);
	bp::def("fractures_ident", &fci::fractures_ident);

	// export compdat
	bp::class_< compdat >("compdat", bp::no_init)
		.def_readwrite("well_name", &compdat::well_name)
		.def_readwrite("branch_name", &compdat::branch_name)
		.def_readwrite("dir", &compdat::dir)
		.def("dump", &print_compdat)
	;

	// export cd_storage as list
	typedef bspy_converter< list_traits< cd_storage, 2 > > py_cds_conv;
	py_cds_conv::register_to_py();
	py_cds_conv::register_from_py();

	void (compdat_builder::*init1)(t_ulong, t_ulong, spv_float, spv_float) = &compdat_builder::init;
	void (compdat_builder::*init2)(smart_ptr< well_pool_iface, true >) = &compdat_builder::init;
	void (compdat_builder::*init3)(t_ulong, t_ulong, sp_obj) = &compdat_builder::init;
	bp::class_< compdat_builder >("compdat_builder",
		bp::init< t_ulong, t_ulong, spv_float, spv_float >())
		.def(bp::init< t_ulong, t_ulong, spv_float, spv_float, smart_ptr< well_pool_iface, true> >())
		.def("init", init1)
		.def("init", init2)
		.def("init", init3)
		.def("build", &compdat_builder::build,
			bp::return_value_policy< bp::copy_const_reference >(),
			cdb_build_overl())
		.def("clear", &compdat_builder::clear)
		.def("storage", &compdat_builder::storage,
			bp::return_value_policy< bp::copy_const_reference >())
	;

	// export fracture
	bp::class_< fracture >("fracture", bp::no_init)
		.def_readwrite("well_name", &fracture::well_name)
		.def_readwrite("branch_name", &fracture::branch_name)
		//.def_readwrite("dir", &compdat::dir)
		.def("dump", &print_fracture)
		;

	// export cd_storage as list
	typedef bspy_converter< list_traits< frac_storage, 2 > > py_fracs_conv;
	py_fracs_conv::register_to_py();
	py_fracs_conv::register_from_py();

	// export fracture_builder
	void (fracture_builder::*init4)(t_ulong, t_ulong, spv_float, spv_float) = &fracture_builder::init;
	void (fracture_builder::*init5)(smart_ptr< well_pool_iface, true >) = &fracture_builder::init;
	void (fracture_builder::*init6)(t_ulong, t_ulong, sp_obj) = &fracture_builder::init;
	bp::class_< fracture_builder >("fracture_builder",
			bp::init< t_ulong, t_ulong, spv_float, spv_float >())
		.def(bp::init< t_ulong, t_ulong, spv_float, spv_float, smart_ptr< well_pool_iface, true> >())
		.def("init", init4)
		.def("init", init5)
		.def("init", init6)
		.def("build", &fracture_builder::build,
			bp::return_value_policy< bp::copy_const_reference >(),
			cdb_build_overl())
		.def("clear", &fracture_builder::clear)
		.def("storage", &fracture_builder::storage,
			bp::return_value_policy< bp::copy_const_reference >())
		;

	// export compl_n_frac_builder
	void (compl_n_frac_builder::*init7)(t_ulong, t_ulong, spv_float, spv_float) = &compl_n_frac_builder::init;
	void (compl_n_frac_builder::*init8)(smart_ptr< well_pool_iface, true >) = &compl_n_frac_builder::init;
	void (compl_n_frac_builder::*init9)(t_ulong, t_ulong, sp_obj) = &compl_n_frac_builder::init;
	bp::class_< compl_n_frac_builder >("compl_n_frac_builder",
			bp::init< t_ulong, t_ulong, spv_float, spv_float >())
		.def(bp::init< t_ulong, t_ulong, spv_float, spv_float, smart_ptr< well_pool_iface, true> >())
		.def("init", init7)
		.def("init", init8)
		.def("init", init9)
		.def("compl_build", &compl_n_frac_builder::compl_build,
			bp::return_value_policy< bp::copy_const_reference >(),
			cdb_build_overl())
		.def("frac_build", &compl_n_frac_builder::frac_build,
			bp::return_value_policy< bp::copy_const_reference >(),
			cdb_build_overl())
		.def("clear", &compl_n_frac_builder::clear)
		.def("storage_compdat", &compl_n_frac_builder::storage_compdat,
			bp::return_value_policy< bp::copy_const_reference >())
		.def("storage_fracture", &compl_n_frac_builder::storage_fracture,
			bp::return_value_policy< bp::copy_const_reference >())
		;
}

} /* blue_sky::python */
} /* blue_sky */

