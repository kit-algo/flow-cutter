#ifndef REFINE_CUT_H
#define REFINE_CUT_H
#include "tiny_id_func.h"
#include "id_multi_func.h"
#include "array_id_func.h"



template<class Tail, class Head>
std::vector<int> cycle_refine_cut(const Tail&tail, const Head&head, std::vector<int>cut){
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	BitIDFunc side(node_count);
	BitIDFunc in_cut(arc_count);
		
	
	auto out_arc = invert_sorted_id_id_func(tail);

	{
		auto back_arc = [&](int xy)->int{
			auto x = tail(xy), y = head(xy);
			for(auto yz:out_arc(y))
				if(head(yz) == x)
					return yz;
			return -1;
		};

		in_cut.fill(false);
		for(auto a:cut){
			in_cut.set(a, true);
			in_cut.set(back_arc(a), true);
		}

		side.fill(0);

		ArrayIDFunc<int>stack(node_count);
		int stack_end = 0;
		auto push = [&](int x){
			if(side(x) == 0){
				side.set(x, 1);	
				stack[stack_end++] = x;
			}
		};

		auto pop = [&]{
			return stack[--stack_end];
		};

		auto is_empty = [&]{
			return stack_end == 0;
		};

		push(0);
		while(!is_empty()){
			auto x = pop();
			for(auto xy:out_arc(x))
				if(!in_cut(xy))
					push(head(xy));
		}
	}

	ArrayIDFunc<int> move_score(node_count);
	move_score.fill(0);

	for(int xy=0; xy<arc_count; ++xy){
		auto x = tail(xy), y = head(xy);
		if(in_cut(xy)){
			++move_score[x];
			++move_score[y];
		}else{
			--move_score[x];
			--move_score[y];
		}
	}

	int side0_cand = -1, side1_cand = -1;

	for(int x=0; x<node_count; ++x){
		if(move_score(x) > 0){
			if(side(x) == 0){
				if(side0_cand==-1)
					side0_cand = x;
			}else{
				if(side1_cand==-1)
					side1_cand = x;
			}
		}
		if(side0_cand != -1 && side1_cand != -1)
			break;
	}

	if(side0_cand != -1 && side1_cand != -1){
		side.set(side0_cand, 1);
		side.set(side1_cand, 0);
	}
	

	cut.clear();
	for(int i=0; i<arc_count; ++i)
		if(side(tail(i)) != side(head(i)))
			cut.push_back(i);

	return std::move(cut);
}

/*
struct RefinedCut{
	std::vector<int>cut;
	int small_side_size;
	int large_side_size;
};

template<class Tail, class Head>
RefinedCut perfectly_balance_cut(const Tail&tail, const Head&head, std::vector<int>cut){
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	BitIDFunc side(node_count);
	side.fill(0);

	auto out_arc = invert_sorted_id_id_func(tail);

	{
		auto back_arc = [&](int xy)->int{
			auto x = tail(xy), y = head(xy);
			for(auto yz:out_arc(y))
				if(head(yz) == x)
					return yz;
			return -1;
		};

		BitIDFunc in_cut(arc_count);
		in_cut.fill(false);
		for(auto a:cut){
			in_cut.set(a, true);
			in_cut.set(back_arc(a), true);
		}

		ArrayIDFunc<int>stack(node_count);
		int stack_end = 0;
		auto push = [&](int x){
			if(side(x) == 0){
				side.set(x, 1);	
				stack[stack_end++] = x;
			}
		};

		auto pop = [&]{
			return stack[--stack_end];
		};

		auto is_empty = [&]{
			return stack_end == 0;
		};

		push(0);
		while(!is_empty()){
			auto x = pop();
			for(auto xy:out_arc(x))
				if(!in_cut(xy))
					push(head(xy));
		}
	}

	bool small_side;
	bool large_side;
	int small_side_size;
	int large_side_size;
	
	{
		int side_node_count[2] = {0,0};
		for(int x=0; x<node_count; ++x)
			++side_node_count[side(x)];
		if(side_node_count[0] <= side_node_count[1]){
			small_side = 0;
			small_side = 1;
			small_side_size = side_node_count[0];
			large_side_size = side_node_count[1];
		} else {
			small_side = 1;
			large_side = 0;
			small_side_size = side_node_count[1];
			large_side_size = side_node_count[0];
		}
	}

	
	for(int threshold = 0; large_side_size - small_side_size > 1; ++threshold){
		for(int x=0; x<node_count && large_side_size - small_side_size > 1; ++x){
			if(side(x) == large_side){
				int large_side_neighbor_count = 0;
				int small_side_neighbor_count = 0;
				for(auto xy:out_arc(x)){
					auto y = head(xy);
					if(side(y) == small_side)
						++small_side_neighbor_count;
					else
						++large_side_neighbor_count;
				}
				if(large_side_neighbor_count <= small_side_neighbor_count + threshold){
					side.set(x, small_side);
					--large_side_size;
					++small_side_size;
				}
			}
		}
	}

	cut.clear();
	for(int i=0; i<arc_count; ++i)
		if(side(tail(i)) == small_side && side(head(i)) == large_side)
			cut.push_back(i);

	return {std::move(cut), small_side_size, large_side_size};
}*/

#endif
