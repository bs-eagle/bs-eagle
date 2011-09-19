/// @file wpi_algo_xaction.h
/// @brief Intersect action class actually track intersections of well & cell boxes
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_PLJRVQ8B
#define WPI_ALGO_XACTION_PLJRVQ8B

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>

#include "wpi_algo_pod.h"
#include "wpi_algo_meshp.h"

#define MD_TOL 0.000001

namespace blue_sky { namespace wpi {

template< class strat_t >
struct wpi_algo_xaction : public wpi_algo_helpers< strat_t > {
	// basic types
	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::cell_pos     cell_pos;

	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;

	// import global consts
	enum { D = strat_t::D, CVN = strat_t::CVN, inner_point_id = strat_t::inner_point_id };

	// import pods
	typedef wpi_algo_pod< strat_t > wpi_pod;
	typedef typename wpi_pod::cell_data cell_data;
	typedef typename wpi_pod::sp_cell_data sp_cell_data;
	typedef typename wpi_pod::trimesh trimesh;
	typedef typename wpi_pod::trim_iterator trim_iterator;
	typedef typename wpi_pod::ctrim_iterator ctrim_iterator;

	typedef typename wpi_pod::well_data well_data;
	typedef typename wpi_pod::well_path well_path;
	typedef typename wpi_pod::wp_iterator wp_iterator;
	typedef typename wpi_pod::cwp_iterator cwp_iterator;

	typedef typename wpi_pod::well_hit_cell well_hit_cell;
	typedef typename wpi_pod::intersect_path intersect_path;

	// mesh_part
	typedef typename wpi_algo_meshp< strat_t >::mesh_part mesh_part;

	/*-----------------------------------------------------------------
	* Box description
	*----------------------------------------------------------------*/
	// structure to help identify given boxes
	class box_handle {
	public:
		enum {
			CELL_BOX,
			WELL_BOX,
			MESH_BOX
		};

		virtual int type() const = 0;

	protected:
		template< class fish_t, class = void >
		struct fish2box_t {
			// default value
			enum { type = CELL_BOX };
			typedef cell_data data_t;
		};
		// overload for well path
		template< class unused >
		struct fish2box_t< wp_iterator, unused > {
			enum { type = WELL_BOX };
			typedef well_data data_t;
		};
		// overload for mesh_part
		template< class unused >
		struct fish2box_t< mesh_part*, unused > {
			enum { type = MESH_BOX };
			typedef mesh_part data_t;
		};
	};
	// pointer is really stored as box handle
	typedef st_smart_ptr< box_handle > sp_bhandle;

	template< class fish >
	class box_handle_impl : public box_handle {
	public:
		typedef fish fish_t;
		typedef typename box_handle::template fish2box_t< fish_t >::data_t data_t;

		box_handle_impl(const fish_t& f) : f_(f) {}

		int type() const {
			return box_handle::template fish2box_t< fish_t >::type;
		}

		fish_t data() const {
			return f_;
		}

	private:
		fish_t f_;
	};
	// handy typedefs
	typedef box_handle_impl< trim_iterator > cell_box_handle;
	typedef box_handle_impl< wp_iterator > well_box_handle;
	typedef box_handle_impl< mesh_part* > mesh_box_handle;

	// box intersections storage
	typedef CGAL::Box_intersection_d::Box_with_handle_d< double, D, sp_bhandle > Box;

	/*-----------------------------------------------------------------
	* callback functor that triggers on intersecting bboxes
	*----------------------------------------------------------------*/
	class intersect_action {
	public:
		typedef typename intersect_path::iterator x_iterator;

		//template< int N >
		struct spatial_sort {
			typedef int dirvec_t[D];
			typedef typename intersect_path::iterator x_iterator;

			spatial_sort(const dirvec_t& dir, const vertex_pos_i& mesh_size)
				: dir_(dir), m_size_(mesh_size)
			{}

			spatial_sort(const spatial_sort& rhs)
				: dir_(rhs.dir_), m_size_(rhs.m_size_)
			{}

			bool operator()(const x_iterator& r1, const x_iterator& r2) const {
				// check if r1 or r2 are node points
				if(r1->is_node) {
					if(!r2->is_node)
						return true;
					else
						// merge node points in same position
						return false;
				}
				else if(r2->is_node)
					return false;

				// cell ids
				vertex_pos_i c1, c2;
				decode_cell_id(r1->cell->first, c1, m_size_);
				decode_cell_id(r2->cell->first, c2, m_size_);

				//bool res = false;
				for(uint i = 0; i < D; ++i) {
					int f = greater(i, c1[i], c2[i]);
					if(f > 0)
						return true;
					else if(f == 0)
						return false;
				}
				return false;
			}

			int greater(uint ndim, ulong v1, ulong v2) const {
				if(v1 == v2) return -1;
				return dir_[ndim] == 0 ? v1 > v2 : v2 > v1;
			}

			const dirvec_t& dir_;
			const vertex_pos_i& m_size_;
		};

		// ctor
		intersect_action(trimesh& mesh, well_path& wp, intersect_path& X, const vertex_pos_i& mesh_size)
			: m_(mesh), wp_(wp), x_(X)
		{
			ca_assign(m_size_, mesh_size);
		}

		static double distance(const Point& p1, const Point& p2) {
			return std::sqrt(Segment(p1, p2).squared_length());
		}

		static double calc_md(const wp_iterator& fish, const Point& target) {
			// walk all segments before the last one;
			double md = fish->second.md();
			// append tail
			md += distance(fish->second.start(), target);
			return md;
		}

		void operator()(const Box& bc, const Box& bw) {
			typedef typename strat_t::xpoints_list xpoints_list;

			//trim_iterator cell_fish = static_cast< cell_box_handle* >(bc.handle().get())->data();
			//wp_iterator well_fish = static_cast< well_box_handle* >(bw.handle().get())->data();
			//cell_data& c = cell_fish->second;
			//well_data& w = well_fish->second;

			cell_box_handle* cell_h = static_cast< cell_box_handle* >(bc.handle().get());
			well_box_handle* well_h = static_cast< well_box_handle* >(bw.handle().get());

			xpoints_list res = strat_t::on_boxes_intersect(cell_h, well_h, bc, bw);

			// add all points to interscetion path
			for(typename xpoints_list::iterator px = res.begin(), end = res.end(); px != end; ++px) {
				x_.insert(well_hit_cell(
					px->first, well_h->data(), cell_h->data(),
					calc_md(well_h->data(), px->first),
					px->second
				));
			}

			//strat_t::on_boxes_intersect< wpi_impl >(bc, bw, x_);
		}

		// run it after all dups killed
		void append_wp_nodes(const std::vector< ulong >& hit_idx) {
			if(!wp_.size()) return;

			// walk through the intersection path and add node points
			// of well geometry to the cell with previous intersection
			x_iterator px = x_.begin();
			//ulong facet_id;

			wp_iterator pw = wp_.begin();
			ulong node_idx = 0;
			for(wp_iterator end = wp_.end(); pw != end; ++pw, ++node_idx) {
				//const well_data& wseg = pw->second;
				// lower_bound
				while(px != x_.end() && px->md < pw->second.md())
					++px;

				px = insert_wp_node(hit_idx[node_idx], pw, px);
			}

			// well path doesn't contain the end-point of trajectory
			// add it manually
			insert_wp_node(hit_idx[node_idx], --pw, x_.end(), true);
		}

		//template< int N >
		void remove_dups2() {
			typedef int dirvec_t[D];
			//typedef intersect_path::iterator x_iterator;

			struct top_surv {
				typedef std::set< x_iterator, spatial_sort > spat_storage_t;
				typedef typename spat_storage_t::iterator spat_iterator;

				top_surv(const dirvec_t& dir, intersect_path& x, const vertex_pos_i& mesh_size)
					: dir_(dir), x_(x), m_size_(mesh_size)
				{}

				x_iterator operator()(x_iterator from, x_iterator to) {
					spat_storage_t sr(spatial_sort(dir_, m_size_));

					// spatially sort iterators
					for(x_iterator px = from; px != to; ++px)
						sr.insert(px);
					// save only frst element
					while(from != to) {
						if(from != *sr.begin())
							x_.erase(from++);
						else ++from;
					}

					return *sr.begin();
				}

				const dirvec_t& dir_;
				intersect_path& x_;
				const vertex_pos_i& m_size_;
			};

			// main processing cycle
			// sanity check
			if(x_.size() < 2) return;

			// position on first intersection
			x_iterator px = x_.begin();
			// walk the nodes and determine direction of trajectory
			double max_md;
			dirvec_t dir;
			for(wp_iterator pw = wp_.begin(), end = wp_.end(); pw != end; ++pw) {
				// identify direction
				const well_data& seg = pw->second;
				// calc direction vector
				Point start = seg.start();
				Point finish = seg.finish();
				for(uint i = 0; i < 2; ++i)
					dir[i] = start[i] <= finish[i] ? 0 : 1;
				// judge
				top_surv judge(dir, x_, m_size_);

				// remove dups lying on current well segment
				max_md = seg.md() + seg.len();
				for(; px != x_.end() && px->md <= max_md; ++px) {
					// skeep well node points if any
					//if(px->facet == 4)
					//	continue;

					// find range of cross points with equal MD
					// upper_bound
					ulong cnt = 0;
					x_iterator pn = px;
					for(; pn != x_.end() && abs(px->md - pn->md) < MD_TOL; ++pn, ++cnt)
					{}
					// if we have nonempty range - leave only 1 element
					if(cnt)
						px = judge(px, pn);
				}
			}
		}

		spv_float export_1d() const {
			spv_float res = BS_KERNEL.create_object(v_float::bs_type());
			res->resize(x_.size() * 7);
			vf_iterator pr = res->begin();

			for(typename intersect_path::const_iterator px = x_.begin(), end = x_.end(); px != end; ++px) {
				// cell num
				*pr++ = px->cell->first;
				// MD
				*pr++ = px->md;
				// intersection point
				for(uint i = 0; i < D; ++i)
					*pr++ = px->where[i];
				// facet id
				*pr++ = px->facet;
				// node flag
				*pr++ = px->is_node;
			}

			return res;
		}

	private:
		x_iterator insert_wp_node(ulong cell_id, wp_iterator pw, x_iterator px, bool end_point = false) {
			// initialization
			const well_data& wseg = pw->second;
			Point where = wseg.start();
			t_float wp_md = wseg.md();
			if(end_point) {
				where = wseg.finish();
				wp_md += wseg.len();
				//wp_md = px->md + distance(px->where, where);
			}

			// check if current or prev intersection match with node
			uint facet_id = inner_point_id;
			if(px != x_.end() && abs(px->md - wp_md) < MD_TOL)
				facet_id = px->facet;
			else if(px != x_.begin()) {
				--px;
				if(abs(px->md - wp_md) < MD_TOL)
					facet_id = px->facet;
			}

			// if node point coinside with existing intersection
			// then just change is_node flag (dont't affect ordering)
			// otherwise insert new intersection point
			if(facet_id == inner_point_id)
				px = x_.insert(well_hit_cell(
					where, pw, m_.find(cell_id), wp_md,
					facet_id, true
				));
			else {
				well_hit_cell& x = const_cast< well_hit_cell& >(*px);
				x.is_node = true;
			}
			return px;
		}

		// mesh
		trimesh& m_;
		// well path
		well_path& wp_;
		// well-mesh intersection path
		intersect_path& x_;
		// mesh size
		vertex_pos_i m_size_;
	};

};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_XACTION_PLJRVQ8B */

