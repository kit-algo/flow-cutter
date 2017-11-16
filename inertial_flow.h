#ifndef INERTIAL_FLOW_H
#define INERTIAL_FLOW_H

#include "array_id_func.h"
#include "id_func.h"
#include "tiny_id_func.h"
#include "back_arc.h"
#include "geo_pos.h"
#include "edmond_karp.h"
#include "dinic.h"
#include <vector>
#include <cassert>

namespace inertial_flow{

	struct Cut{
		BitIDFunc is_on_smaller_side;
		int smaller_side_size;
		int cut_size;
	};

	template<class InvTail, class Head, class Flow, class Source>
	Cut extract_cut_from_maximum_unit_flow(
		const InvTail&inv_tail, const Head&head,
		const Flow&flow,
		const Source&source_list
	){
		const int node_count = head.image_count();
		
		BitIDFunc reachable(node_count);
		int reachable_count = 0;	
		reachable.fill(false);

		ArrayIDFunc<int>stack(node_count);
		int stack_end = 0;

		for(int i=0; i<source_list.preimage_count(); ++i){
			stack[stack_end++] = source_list(i);
			reachable.set(source_list(i), true);
			++reachable_count;
		}


		while(stack_end != 0){
			auto x = stack[--stack_end];
			for(auto xy:inv_tail(x)){
				if(!flow(xy)){
					auto y = head(xy);
					if(!reachable(y)){
						reachable.set(y, true);
						stack[stack_end++] = y;
						++reachable_count;
					}
				}
			}
		}

		int cut_size = 0;
		for(int x=0; x<node_count; ++x)
			for(auto xy:inv_tail(x))
				if(reachable(x) != reachable(head(xy)))
					++cut_size;
		cut_size /= 2;

		if(reachable_count <= node_count/2)
			return Cut{std::move(reachable), reachable_count, cut_size};
		else
			return Cut{~std::move(reachable), node_count-reachable_count, cut_size};
	}

	template<class CompGeoPos>
	void build_source_and_target_list(int node_count, double min_balance, const CompGeoPos&comp_geo_pos, ArrayIDIDFunc&source_list, ArrayIDIDFunc&target_list){
		ArrayIDIDFunc node_order = identity_permutation(node_count);
		std::sort(node_order.begin(), node_order.end(), comp_geo_pos);
		int min_side_size = std::max(static_cast<int>(min_balance*node_count), 1);

		source_list = id_id_func(min_side_size, node_count, [&](int x){return node_order[x];});
		target_list = id_id_func(min_side_size, node_count, [&](int x){return node_order[node_count-x-1];});

	}

	template<class InvTail, class Head, class BackArc, class GetGeoPos>
	Cut compute_inertial_flow_cut(
		const InvTail&inv_tail, const Head&head, const BackArc&back_arc,
		const GetGeoPos&geo_pos, 
		double min_balance
	){
		const int node_count = head.image_count();
		ArrayIDIDFunc 
			source_list[4], 
			target_list[4];

		build_source_and_target_list(
			node_count, min_balance, 
			[&](int l, int r)->bool{return geo_pos(l).lon < geo_pos(r).lon;}, 
			source_list[0], target_list[0]
		);
		build_source_and_target_list(
			node_count, min_balance, 
			[&](int l, int r)->bool{return geo_pos(l).lat < geo_pos(r).lat;}, 
			source_list[1], target_list[1]
		);
		build_source_and_target_list(
			node_count, min_balance, 
			[&](int l, int r)->bool{return geo_pos(l).lon+geo_pos(l).lat < geo_pos(r).lon+geo_pos(r).lat;}, 
			source_list[2], target_list[2]
		);
		build_source_and_target_list(
			node_count, min_balance, 
			[&](int l, int r)->bool{return geo_pos(l).lon-geo_pos(l).lat < geo_pos(r).lon-geo_pos(r).lat;}, 
			source_list[3], target_list[3]
		);

		max_flow::UnitDinicAlgo<InvTail, Head, BackArc, ArrayIDIDFunc, ArrayIDIDFunc> instance [] = { 
			{inv_tail, head, back_arc, source_list[0], target_list[0]},
			{inv_tail, head, back_arc, source_list[1], target_list[1]},
			{inv_tail, head, back_arc, source_list[2], target_list[2]},
			{inv_tail, head, back_arc, source_list[3], target_list[3]}
		};

		for(;;){
			int next_instance = 0;
			for(int i=1; i<4; ++i)
				if(instance[i].get_current_flow_intensity() < instance[next_instance].get_current_flow_intensity())
					next_instance = i;

			if(instance[next_instance].is_finished())
				return extract_cut_from_maximum_unit_flow(inv_tail, head, instance[next_instance].move_saturated_flags(), source_list[next_instance]);
			instance[next_instance].advance();
		}
	}

	template<class Tail, class Head, class GetGeoPos>
	Cut compute_inertial_flow_cut(
		const Tail&tail, const Head&head, 
		const GetGeoPos&geo_pos, 
		double min_balance
	){
		if(std::is_sorted(tail.begin(), tail.end()))
			return compute_inertial_flow_cut(invert_sorted_id_id_func(tail), head, compute_back_arc_permutation(tail, head), geo_pos, min_balance);
		else
			return compute_inertial_flow_cut(invert_id_id_func(tail), head, compute_back_arc_permutation(tail, head), geo_pos, min_balance);
	}

	template<class Tail, class Head, class GetGeoPos>
	std::vector<int> compute_inertial_flow_separator(const Tail&tail, const Head&head, const GetGeoPos&geo_pos, double min_balance){
		const int arc_count = head.preimage_count();	
		const int node_count = head.image_count();

		std::vector<int>sep;
		if(node_count == 1){
			sep = {0};
		} else {
			Cut c = compute_inertial_flow_cut(tail, head, geo_pos, min_balance);
			
			for(int i=0; i<arc_count; ++i)
				if(c.is_on_smaller_side(tail(i)) && !c.is_on_smaller_side(head(i)))
					sep.push_back(head(i));
		}
		return sep; // NVRO
	}

	template<class GetGeoPos>
	struct InertialFlowSeparator{
		InertialFlowSeparator(const GetGeoPos&geo_pos, double min_balance):
			geo_pos(&geo_pos), min_balance(min_balance){}

		template<class Tail, class Head, class InputNodeID, class ArcWeight>
		std::vector<int>operator()(const Tail&tail, const Head&head, const InputNodeID& input_node_id, const ArcWeight&arc_weight)const{
			const int node_count = head.image_count();
			return compute_inertial_flow_separator(
				tail, head, 
				id_func(node_count, [&](int x){return (*geo_pos)(input_node_id(x));}),
				min_balance
			);
		}

		const GetGeoPos*geo_pos;
		double min_balance;
	};

	template<class GetGeoPos>
	InertialFlowSeparator<GetGeoPos>
		ComputeSeparator(const GetGeoPos&geo_pos, double min_balance){
		return {geo_pos, min_balance};
	}
}

#endif

