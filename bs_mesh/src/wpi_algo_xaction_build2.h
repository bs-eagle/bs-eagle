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
#include "wpi_algo_xaction.h"

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

	enum { D = strat_t::D };
	typedef typename base_t::template xbbox< D > xbbox_t;

	// base functions
	using base_t::check_intersection;

	// base members
	using base_t::hit_idx_;
	using base_t::m_;
	using base_t::m_size_;
	using base_t::wp_;

	/*-----------------------------------------------------------------
	* Box description
	*----------------------------------------------------------------*/
	enum fish_id {
		CELL_BOX,
		WELL_BOX,
		MESH_BOX
	};

#ifdef _WIN32
	template< fish_id Id >
	struct fish2box {
		// default value
		typedef ulong fish_t;
		typedef cell_data data_t;
	};
	// overload for well path
	template< >
	struct fish2box< WELL_BOX > {
		typedef ulong fish_t;
		typedef well_data data_t;
	};
	// overload for mesh_part
	template< >
	struct fish2box< MESH_BOX > {
		typedef const mesh_part* fish_t;
		typedef mesh_part data_t;
	};
#else
	template< fish_id Id, class = void >
	struct fish2box {
		// default value
		typedef ulong fish_t;
		typedef cell_data data_t;
	};
	// overload for well path
	template< class unused >
	struct fish2box< WELL_BOX, unused > {
		typedef ulong fish_t;
		typedef well_data data_t;
	};
	// overload for mesh_part
	template< class unused >
	struct fish2box< MESH_BOX, unused > {
		typedef const mesh_part* fish_t;
		typedef mesh_part data_t;
	};
#endif
	// structure to help identify given boxes
	class box_handle {
	public:
		virtual int type() const = 0;
	};
	// pointer is really stored as box handle
	typedef st_smart_ptr< box_handle > sp_bhandle;

	template< fish_id Id >
	class box_handle_impl : public box_handle {
		typedef fish2box< Id > fish2box_t;

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

	private:
		fish_t f_;
	};
	// handy typedefs
	typedef box_handle_impl< CELL_BOX > cell_box_handle;
	typedef box_handle_impl< WELL_BOX > well_box_handle;
	typedef box_handle_impl< MESH_BOX > mesh_box_handle;

	// box intersections storage
	typedef CGAL::Box_intersection_d::Box_with_handle_d< double, D, sp_bhandle > Box;

	/*-----------------------------------------------------------------
	 * action taken when boxes intersection detected
	 *----------------------------------------------------------------*/
	struct leafs_builder {
		typedef typename mesh_part::container_t parts_container;
		typedef std::vector< Segment > Segments;
		//typedef xbbox< D > xbbox_t;

		leafs_builder(intersect_builder2& A, const Segments& s, std::vector< ulong >& hit, parts_container& leafs)
			: A_(A), s_(s), hit_(hit), leafs_(leafs)
		{}

		//leafs_builder(const leafs_builder& rhs)
		//	: A_(rhs.A_), s_(rhs.s_), hit_(rhs.hit_), leafs_(rhs.leafs_)
		//{}

		void operator()(const Box* mb, const Box* wb) {
			const mesh_part* pm = static_cast< mesh_box_handle* >(mb->handle().get())->data();
			ulong wseg_id = static_cast< well_box_handle* >(wb->handle().get())->data();

			// if mesh_part contains > 1 cells then just push it for further splitting
			// otherwise check for real intersections
			if(pm->size() == 1) {
				ulong cell_id = pm->ss_id(0);
				cell_data& cell = A_.m_[cell_id];
				// check if segment start & finish are inside the cell
				if(hit_[wseg_id] >= A_.m_.size()) {
					const Point& start = A_.wp_[wseg_id].start();
					if(
						mesh_tools_t::point_inside_bbox(cell.bbox(), start) &&
						cell.contains(start)
						)
						hit_[wseg_id] = cell_id;
				}
				// check segment end for last segment
				if(wseg_id == A_.wp_.size() - 1 && hit_[wseg_id + 1] >= A_.m_.size()) {
					const Point& finish = A_.wp_[wseg_id].finish();
					if(
						mesh_tools_t::point_inside_bbox(cell.bbox(), finish) &&
						cell.contains(finish)
						)
						hit_[wseg_id + 1] = cell_id;
				}

				// check for intersections with segment
				//if(!s_[wseg_id].is_degenerate() && strat_t::bbox_segment_x(cell.bbox(), s_[wseg_id]))
				if(!s_[wseg_id].is_degenerate() && CGAL::do_intersect(xbbox_t::get(cell), s_[wseg_id]))
					A_.check_intersection(cell_id, wseg_id, s_[wseg_id]);
			}
			else
				leafs_.insert(*pm);
		}

		intersect_builder2& A_;
		const Segments& s_;
		std::vector< ulong >& hit_;
		parts_container& leafs_;
	};

	/*-----------------------------------------------------------------
	 * main body
	 *----------------------------------------------------------------*/
	// propagate ctor
	intersect_builder2(trimesh& mesh, well_path& wp, const vertex_pos_i& mesh_size)
		: base_t(mesh, wp, mesh_size)
	{}

	// branch & bound algorithm for finding cells that really intersect with well
	// works using boxes intersection
	hit_idx_t& build() {
		typedef typename mesh_part::container_t parts_container;
		typedef typename parts_container::iterator part_iterator;
		typedef std::vector< Segment > Segments;

		// cache well segments
		// and make boxes for them
		Segments wseg(wp_.size());
		std::vector< Box > well_boxes(wp_.size());
		std::vector< Box* > well_boxes_p(wp_.size());
		//wseg.reserve(wp_.size());
		for(ulong i = 0; i < wp_.size(); ++i) {
			wseg[i] = wp_[i].segment();
			well_boxes[i] = Box(wp_[i].bbox(), new well_box_handle(i));
			well_boxes_p[i] = &well_boxes[i];
		}

		// list of mesh parts
		parts_container parts;
		parts.insert(mesh_part(m_, m_size_));

		// result - hit index
		hit_idx_.resize(wp_.size() + 1);
		std::fill(hit_idx_.begin(), hit_idx_.end(), m_.size());

		// let's go
		leafs_builder B(*this, wseg, hit_idx_, parts);
		//std::vector< Box > mp_boxes;
		while(parts.size()) {
			// we need container to hold all mesh parts
			parts_container leafs;

			// split every part
			// and make boxes around splitted parts
			// try to use single boxes container over all iterations - no profit
			const ulong max_boxes = parts.size() * (1 << D);
			//if(mp_boxes.size() < max_boxes)
			//	mp_boxes.resize(max_boxes);
			std::vector< Box > mp_boxes;
			std::vector< Box* > mp_boxes_p;
			mp_boxes.reserve(max_boxes);
			mp_boxes_p.reserve(max_boxes);

			//ulong i = 0;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				parts_container kids = p->divide();
				//ulong i = 0;
				for(part_iterator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
					mp_boxes.push_back(Box(pk->bbox(), new mesh_box_handle(&*leafs.insert(*pk).first)));
					mp_boxes_p.push_back(&mp_boxes[mp_boxes.size() - 1]);
				}
			}

			// do intersect boxes
			//parts_container new_kids;
			parts.clear();
			CGAL::box_intersection_d(
				mp_boxes_p.begin(), mp_boxes_p.end(),
				well_boxes_p.begin(), well_boxes_p.end(),
				B
			);

			// update parts
			//parts = new_kids;
		}

		return hit_idx_;
	}
};

}} /* blue_sky::wpi */


#endif /* end of include guard: WPI_ALGO_XACTION_BUILD2_9FCQZ4H8 */

