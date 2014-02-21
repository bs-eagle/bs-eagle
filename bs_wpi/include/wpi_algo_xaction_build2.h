/// @file wpi_algo_xaction_build2.h
/// @brief Third approach to finding intersection points
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_BUILD2_9FCQZ4H8
#define WPI_ALGO_XACTION_BUILD2_9FCQZ4H8

#include <CGAL/box_intersection_d.h>
#include "loki/TypeManip.h"
#include "wpi_algo_xaction.h"

#include <boost/pool/object_pool.hpp>

namespace blue_sky { namespace wpi {

template< class strat_t >
class intersect_builder2 : public intersect_base< strat_t > {
public:
	typedef intersect_base< strat_t > base_t;

	// inherit types
	typedef typename base_t::cell_data    cell_data;
	typedef typename base_t::well_data    well_data;
	typedef typename base_t::trimesh      trimesh;
	typedef typename base_t::well_path    well_path;
	typedef typename base_t::mesh_part    mesh_part;
	typedef typename base_t::hit_idx_t    hit_idx_t;

	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Point        Point;
	typedef typename strat_t::Segment      Segment;

	typedef typename base_t::mesh_tools_t mesh_tools_t;
	typedef typename mesh_part::container_t parts_container;
	typedef typename parts_container::iterator part_iterator;
	typedef typename parts_container::const_iterator part_citerator;

	enum { D = strat_t::D };
	typedef typename base_t::template xbbox< D > xbbox_t;

	// base functions
	using base_t::check_intersection;

	// base members
	using base_t::hit_idx_;
	using base_t::m_;
	//using base_t::m_size_;
	using base_t::wp_;

	/*-----------------------------------------------------------------
	* Box description
	*----------------------------------------------------------------*/
	enum fish_id {
		CELL_BOX = 0,
		WELL_BOX = 1,
		MESH_BOX = 2
	};

	template< class Id, class = void >
	struct fish2box {
		// default value
		typedef ulong fish_t;
		typedef cell_data data_t;

		static bool less(const fish_t& lhs, const fish_t& rhs) {
			return lhs < rhs;
		}
	};
	// overload for well path
	template< class unused >
	struct fish2box< Loki::Int2Type< 1 >, unused > {
		typedef ulong fish_t;
		typedef well_data data_t;

		static bool less(const fish_t& lhs, const fish_t& rhs) {
			return lhs < rhs;
		}
	};
	// overload for mesh_part
	template< class unused >
	struct fish2box< Loki::Int2Type< 2 >, unused > {
		typedef const mesh_part* fish_t;
		typedef mesh_part data_t;

		static bool less(const fish_t& lhs, const fish_t& rhs) {
			return *lhs < *rhs;
		}
	};

	// structure to help identify given boxes
	class box_handle {
	public:
		virtual ~box_handle() {};
		virtual int type() const = 0;
		virtual bool less(const box_handle& rhs) const = 0;
		// for storing into ordered containers
		bool operator<(const box_handle& rhs) const {
			return less(rhs);
		}
	};
	// pointer is really stored as box handle
	typedef st_smart_ptr< box_handle > sp_bhandle;

	template< int Id >
	class box_handle_impl : public box_handle {
		typedef fish2box< Loki::Int2Type< Id > > fish2box_t;

	public:
		enum { id = Id };
		typedef typename fish2box_t::fish_t fish_t;
		typedef typename fish2box_t::data_t data_t;

		box_handle_impl(const fish_t& f) : f_(f) {}

		int type() const {
			return int(id);
		}

		fish_t data() const {
			return f_;
		}

		// for storing into ordered containers
		bool less(const box_handle& rhs) const {
			if(rhs.type() != id) return false;
			return fish2box_t::less(f_, static_cast< const box_handle_impl& >(rhs).f_);
		}

	private:
		fish_t f_;
	};
	// handy typedefs
	typedef box_handle_impl< CELL_BOX > cell_box_handle;
	typedef box_handle_impl< WELL_BOX > well_box_handle;
	typedef box_handle_impl< MESH_BOX > mesh_box_handle;

	// box intersections storage
	typedef CGAL::Box_intersection_d::Box_with_handle_d< double, D, box_handle* > Box;
	//typedef CGAL::Box_intersection_d::Box_with_handle_d< double, D, sp_bhandle > Box;

	/*-----------------------------------------------------------------
	 * action taken when boxes intersection detected
	 *----------------------------------------------------------------*/
	struct leafs_builder {
		typedef std::vector< Segment > Segments;

		leafs_builder(base_t& A, const Segments& s, std::vector< ulong >& hit,
			parts_container& leafs, const ulong min_split_threshold = 0
		)
			: A_(A), s_(s), hit_(hit), leafs_(leafs), m_size_(A.m_.size_flat())
		{
			if(min_split_threshold == 0)
				min_spl_thresh_ = 10 * (1 << D);
			else
				min_spl_thresh_ = min_split_threshold;
		}

		//leafs_builder(const leafs_builder& rhs)
		//	: A_(rhs.A_), s_(rhs.s_), hit_(rhs.hit_), leafs_(rhs.leafs_)
		//{}

		void operator()(const Box* mb, const Box* wb) {
			// Box handles are raw pointers
			const mesh_part* pm = static_cast< mesh_box_handle* >(mb->handle())->data();
			ulong wseg_id = static_cast< well_box_handle* >(wb->handle())->data();
			// Box handlex are st_smart_ptrs
			//const mesh_part* pm = static_cast< mesh_box_handle* >(mb->handle().get())->data();
			//ulong wseg_id = static_cast< well_box_handle* >(wb->handle().get())->data();

			// if size of mesh_part <= min split threshold, then every cell in that mesh part
			// will be checked for intersections with given well path segment
			// instead of marking for further split
			if(pm->size() > min_spl_thresh_) {
				leafs_.insert(*pm);
				return;
			}

			// otherwise check for real intersections in every cell of mesh part
			for(ulong i = 0, sz = pm->size(); i < sz; ++i) {
				const ulong cell_id = pm->ss_id(i);
				cell_data cell = A_.m_[cell_id];
				// check if segment start & finish are inside the cell
				if(hit_[wseg_id] >= m_size_) {
					const Point& start = A_.wp_[wseg_id].start();
					if(
						mesh_tools_t::point_inside_bbox(cell.bbox(), start) &&
						cell.contains(start)
					)
						hit_[wseg_id] = cell_id;
				}
				// check segment end for last segment
				if(wseg_id == A_.wp_.size() - 1 && hit_[wseg_id + 1] >= m_size_) {
					const Point& finish = A_.wp_[wseg_id].finish();
					if(
						mesh_tools_t::point_inside_bbox(cell.bbox(), finish) &&
						cell.contains(finish)
					)
						hit_[wseg_id + 1] = cell_id;
				}

				// check for intersections with segment
				// update: assume that in 95% of cases cells are rectangular,
				// and intersection most probably exists
				//if(!s_[wseg_id].is_degenerate() && strat_t::bbox_segment_x(cell.bbox(), s_[wseg_id]))
				if(!s_[wseg_id].is_degenerate() && CGAL::do_intersect(xbbox_t::get(cell), s_[wseg_id]))
				//if(!s_[wseg_id].is_degenerate())
					A_.check_intersection(cell_id, wseg_id, s_[wseg_id]);
			}
		}

		base_t& A_;
		const Segments& s_;
		std::vector< ulong >& hit_;
		parts_container& leafs_;
		const ulong m_size_;
		// store mesh parts for detailed check here
		std::set< mesh_part* > parts2check_;
		// splitting threshold
		ulong min_spl_thresh_;
	};

	/*-----------------------------------------------------------------
	 * main body
	 *----------------------------------------------------------------*/
	// propagate ctor
	intersect_builder2(trimesh& mesh, well_path& wp)
		: base_t(mesh, wp)
	{}

	// branch & bound algorithm for finding cells that really intersect with well
	// works using boxes intersection
	hit_idx_t& build() {
		typedef std::vector< Segment > Segments;

		Segments wseg(wp_.size());
		std::vector< Box > well_boxes(wp_.size());
		std::vector< Box* > well_boxes_p(wp_.size());
		// pool for allocating well box handles
		// they all will be detroyed automatically on exit
		boost::object_pool< well_box_handle > wbh_pool;
		// cache well segments
		// and make boxes for them
		for(ulong i = 0; i < wp_.size(); ++i) {
			wseg[i] = wp_[i].segment();
			well_boxes[i] = Box(wp_[i].bbox(), new(wbh_pool.malloc()) well_box_handle(i));
			//well_boxes[i] = Box(wp_[i].bbox(), new well_box_handle(i));
			well_boxes_p[i] = &well_boxes[i];
		}

		// list of mesh parts
		parts_container parts;
		parts.insert(mesh_part(m_));

		// result - hit index
		hit_idx_.resize(wp_.size() + 1);
		std::fill(hit_idx_.begin(), hit_idx_.end(), m_.size_flat());

		// actual intersector object
		leafs_builder B(*this, wseg, hit_idx_, parts);
		// pool for allocating mesh box handles
		// they all will be detroyed automatically on exit
		boost::object_pool< mesh_box_handle > mbh_pool;

		// let's go
		while(parts.size()) {
			// we need container to hold all mesh parts
			// list can be used because all mesh_parts in parts
			// are different at this point, so results of division will also
			// be unqie
			std::list< mesh_part > leafs;

			// split every part
			// and make boxes around splitted parts
			// try to use single boxes container over all iterations - no profit
			const ulong max_boxes = parts.size() * (1 << D);
			std::vector< Box > mp_boxes;
			std::vector< Box* > mp_boxes_p;
			mp_boxes.reserve(max_boxes);
			mp_boxes_p.reserve(max_boxes);

			//ulong i = 0;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				const parts_container kids = p->divide();
				//ulong i = 0;
				for(part_citerator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
					mp_boxes.push_back(
						Box(pk->bbox(),
						new(mbh_pool.malloc()) mesh_box_handle(&*leafs.insert(leafs.end(),*pk)))
					);
					mp_boxes_p.push_back(&mp_boxes[mp_boxes.size() - 1]);
				}
			}

			// clear parts container before storing boxes to split on next step
			parts.clear();
			// do intersect boxes
			CGAL::box_intersection_d(
				mp_boxes_p.begin(),   mp_boxes_p.end(),
				well_boxes_p.begin(), well_boxes_p.end(),
				B
			);
		}

		return hit_idx_;
	}
};

}} /* blue_sky::wpi */


#endif /* end of include guard: WPI_ALGO_XACTION_BUILD2_9FCQZ4H8 */

