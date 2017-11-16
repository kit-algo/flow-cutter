#ifndef SMALL_TREE_WIDTH_ORDER_H
#define SMALL_TREE_WIDTH_ORDER_H

#include "min_fill_in.h"
#include "contraction_graph.h"
#include "node_flow_cutter.h"
#include "flow_cutter.h"

namespace small_tree_width{

	template<class ComputeOrder>
	ArrayIDIDFunc eliminate_simplicial_nodes(
		ArrayIDIDFunc tail, ArrayIDIDFunc head,
		ArrayIDIDFunc input_node_id,
		const ComputeOrder&compute_order
	){
		assert(is_connected(tail, head));
		assert(std::is_sorted(tail.begin(), tail.end()));

		const int node_count = tail.image_count();
		const int arc_count = tail.preimage_count();
		auto out_arc = invert_sorted_id_id_func(tail);

		int order_pos = 0;
		ArrayIDIDFunc order(node_count, input_node_id.image_count());

		BitIDFunc was_eliminated(node_count);
		was_eliminated.fill(false);

		ArrayIDFunc<int>deg(node_count);
		deg.fill(0);
		for(int xy=0; xy<arc_count; ++xy)
			++deg[head(xy)];

		// Eliminate degree 0 and degree 1 nodes
		{
			ArrayIDIDFunc stack(node_count, node_count);
			int stack_end = 0;
			for(int x=0; x<node_count; ++x){
				if(deg[x] <= 1){
					stack[stack_end++] = x;
				}
			}

			
			while(stack_end != 0){
				int x = stack[--stack_end];
				order[order_pos++] = input_node_id(x);
				was_eliminated.set(x, true);
				for(int xy:out_arc(x)){
					int y = head(xy);
					if(!was_eliminated(y)){
						--deg[y];
						if(deg[y] <= 1){
							stack[stack_end++] = y;
						}
					}
				}
			}
		}


		#ifndef NDEBUG
		for(int x=0; x<node_count; ++x){
			if(!was_eliminated(x))
				assert(deg[x] > 1);
		}
		{
			ArrayIDFunc<int>correct_deg(node_count);
			correct_deg.fill(0);

			for(int xy=0; xy<arc_count; ++xy)
				if(!was_eliminated(tail(xy)))
					++correct_deg[head(xy)];
			for(int x=0; x<node_count; ++x)
				if(!was_eliminated(x))
					assert(deg[x] == correct_deg[x]);

		}
		#endif

		// Eliminate higher degree simplicial nodes
		{

			ArrayIDFunc<std::vector<int>>neighbors(node_count);
			for(int xy=0; xy<arc_count; ++xy){
				int x = tail(xy), y = head(xy);	
				if(!was_eliminated(x) && !was_eliminated(y))
					neighbors[x].push_back(y);
			}
	
			for(int x=0; x<node_count; ++x)
				neighbors[x].push_back(x);

			for(auto&x:neighbors){
				std::sort(x.begin(), x.end());
			}
		

			auto prune_neighbors = [&](int x){
				neighbors[x].erase(
					std::remove_if(
						neighbors[x].begin(),
						neighbors[x].end(),
						was_eliminated
					),
					neighbors[x].end()
				);
			};

			auto is_simplicial = [&](int x){
				prune_neighbors(x);

				for(int y:neighbors[x]){
					if(y == x)
						continue;
					prune_neighbors(y);

					if(
						!std::includes(
							neighbors[y].begin(), neighbors[y].end(),
							neighbors[x].begin(), neighbors[x].end()
						)
					)
						return false;
				}
				return true;
			};

			bool was_simplicial_node_found;
			do{
				was_simplicial_node_found = false;
				for(int x=0; x<node_count; ++x){
					if(!was_eliminated(x)){
						if(is_simplicial(x)){
							order[order_pos++] = input_node_id(x);

							was_eliminated.set(x, true);
							was_simplicial_node_found = true;
						}
					}
				}
			}while(was_simplicial_node_found);
		}

		if(order_pos == node_count)
			return order;

		// Remove eliminated nodes and use other ordering strategy 
		{
			BitIDFunc node_keep_flag = id_func(node_count, [&](int x){return !was_eliminated(x);});

			int new_node_count = count_true(node_keep_flag);

			BitIDFunc arc_keep_flag(tail.preimage_count());
			for(int i=0; i<tail.preimage_count(); ++i)
				arc_keep_flag.set(i, node_keep_flag(tail(i)) && node_keep_flag(head(i)));
			int new_arc_count = count_true(arc_keep_flag);

			tail = keep_if(arc_keep_flag, new_arc_count, move(tail));
			head = keep_if(arc_keep_flag, new_arc_count, move(head));

			auto node_keep_perm = compute_keep_function(node_keep_flag, new_node_count); 
			head = chain(std::move(head), node_keep_perm);
			tail = chain(std::move(tail), node_keep_perm);
			input_node_id = keep_if(node_keep_flag, new_node_count, std::move(input_node_id));

			for(auto x:compute_order(tail, head, input_node_id))
				order[order_pos++] = x;
		}

		return order;
	}




	template<class ComputeOrder>
	ArrayIDIDFunc compute_order_by_decomposing_along_articulation_points(ArrayIDIDFunc tail, ArrayIDIDFunc head, ArrayIDIDFunc input_node_id, const ComputeOrder&compute_order){

		// This is suboptimal

		int node_count = tail.image_count();
		int arc_count = tail.preimage_count();

		ArrayIDIDFunc arc_component;

		{
			auto out_arc = invert_sorted_id_id_func(tail);
			auto back_arc = compute_back_arc_permutation(tail, head);
			arc_component = compute_biconnected_components(out_arc, head, back_arc);
			auto largest_component = max_preimage_over_id_func(compute_histogram(arc_component));

			for(auto&x:arc_component){
				if(x == arc_component.image_count()-1)
					x = largest_component;
				else if(x == largest_component)
					x = arc_component.image_count()-1;
			}
		}

		ArrayIDIDFunc node_max_component(node_count, arc_component.image_count());
		node_max_component.fill(0);
		
		for(int xy=0; xy<arc_count; ++xy){
			max_to(node_max_component[tail(xy)], arc_component(xy));
			max_to(node_max_component[head(xy)], arc_component(xy));
		}


		// Remove all arcs that enter or leave the largest biconnected component

		{
			BitIDFunc keep_flag = id_func(
				arc_count, 
				[&](int a){
					return node_max_component(tail(a)) == node_max_component(head(a));
				}
			);

			arc_count = count_true(keep_flag);
			tail = keep_if(keep_flag, arc_count, std::move(tail));
			head = keep_if(keep_flag, arc_count, std::move(head));
		}

		BitIDFunc node_in_largest_biconnected_component = id_func(node_count, [&](int x){return node_max_component(x) == node_max_component.image_count()-1;});


		auto orderer = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id, ArrayIDFunc<int>a_weight){
			(void)a_weight;
			return compute_order(std::move(a_tail), std::move(a_head), std::move(a_input_node_id));
		};

		auto weight = id_func(tail.preimage_count(), [](int){return 0;});
		return cch_order::reorder_nodes_in_preorder_and_compute_unconnected_graph_order_if_component_is_non_trivial(
			std::move(tail), std::move(head), 
			std::move(input_node_id), std::move(weight),
			orderer, std::move(node_in_largest_biconnected_component)
		);
	}

	template<class ComputeOrder>
	ArrayIDIDFunc compute_order_by_applying_one_round_of_reduction_rules(
		ArrayIDIDFunc tail, ArrayIDIDFunc head,
		ArrayIDIDFunc input_node_id,
		const ComputeOrder&compute_order
	){

		auto orderer = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id, ArrayIDFunc<int>a_weight){
			(void)a_weight;
			return compute_order(std::move(a_tail), std::move(a_head), std::move(a_input_node_id));
		};

		return eliminate_simplicial_nodes(
			std::move(tail), std::move(head), std::move(input_node_id),
			[&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
				auto weight = id_func(a_tail.preimage_count(), [](int){return 0;});
				return cch_order::compute_graph_order_with_degree_two_chain_at_the_begin(
					std::move(a_tail), std::move(a_head), std::move(a_input_node_id), std::move(weight), 
					orderer 
				);
			}
		);
	}

	template<class ComputeOrder>
	ArrayIDIDFunc compute_order_by_applying_reduction_rules(
		ArrayIDIDFunc tail, ArrayIDIDFunc head,
		ArrayIDIDFunc input_node_id,
		const ComputeOrder&compute_order
	){

		auto orderer5 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_one_round_of_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), compute_order);
		};

		auto orderer4 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_one_round_of_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer5);
		};

		auto orderer3 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_one_round_of_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer4);
		};

		auto orderer2 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_one_round_of_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer3);
		};

		auto orderer1 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_one_round_of_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer2);
		};

		return compute_order_by_applying_one_round_of_reduction_rules(std::move(tail), std::move(head), std::move(input_node_id), orderer1);
	}



	template<class ComputeOrder>
	ArrayIDIDFunc compute_order_by_upholding_separator(
		ArrayIDIDFunc tail, ArrayIDIDFunc head, ArrayIDIDFunc input_node_id,
		std::vector<int>separator,
		const ComputeOrder&compute_order
	){
		const int node_count = tail.image_count();
		const int arc_count = tail.preimage_count();

		BitIDFunc in_separator(node_count);
		in_separator.fill(false);
		for(auto x:separator)
			in_separator.set(x, true);

		BitIDFunc keep_arc_flag = id_func(
			arc_count, 
			[&](int a){
				return in_separator(tail(a)) == in_separator(head(a));
			}
		);

		int new_arc_count = count_true(keep_arc_flag);
		tail = keep_if(keep_arc_flag, new_arc_count, std::move(tail));
		head = keep_if(keep_arc_flag, new_arc_count, std::move(head));

		auto orderer = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id, ArrayIDFunc<int>a_weight){
			(void)a_weight;

			static int recurse_level = 0;
			++recurse_level;
			ArrayIDIDFunc order;
			if(recurse_level < 3)
				order = compute_order(std::move(a_tail), std::move(a_head), std::move(a_input_node_id));
			else
				order = chain(compute_greedy_min_shortcut_order(a_tail, a_head), a_input_node_id);
			--recurse_level;
			return order;
		};
		auto weight = id_func(tail.preimage_count(), [](int){return 0;});
		return cch_order::reorder_nodes_in_preorder_and_compute_unconnected_graph_order_if_component_is_non_trivial(
			std::move(tail), std::move(head), 
			std::move(input_node_id), std::move(weight), 
			orderer, std::move(in_separator)
		);
	}


	inline
	int compute_tree_width(ArrayIDIDFunc tail, ArrayIDIDFunc head, ArrayIDIDFunc order){

		auto inv_order = inverse_permutation(order);
		tail = chain(tail, inv_order);
		head = chain(head, inv_order);


		int max_up_deg = 0;
		int current_tail = -1;
		int current_tail_up_deg = 0;
		compute_chordal_supergraph(
			tail, head, 
			[&](int x, int y){
				if(current_tail != x){
					current_tail = x;
					max_to(max_up_deg, current_tail_up_deg);
					current_tail_up_deg = 0;
				}
				++current_tail_up_deg;
			}
		);
		return max_up_deg;
	}

	class ComputeSeparatorSet{
	public:
		explicit ComputeSeparatorSet(flow_cutter::Config config):config(config){}

		template<class Tail, class Head>
		std::vector<std::vector<int>> operator()(const Tail&tail, const Head&head)const{

			const int node_count = tail.image_count();
			const int arc_count = tail.preimage_count();

			auto out_arc = invert_sorted_id_id_func(tail);
			auto back_arc = compute_back_arc_permutation(tail, head);

			auto arc_weight = id_func(arc_count, [](int){ return 0; });

			auto expanded_graph = flow_cutter::expanded_graph::make_graph(
				make_const_ref_id_id_func(tail), 
				make_const_ref_id_id_func(head), 
				make_const_ref_id_id_func(back_arc), 
				make_const_ref_id_id_func(arc_weight), 
				make_const_ref_id_func(out_arc)
			);		

			auto cutter = flow_cutter::make_simple_cutter(expanded_graph, config);
			auto pairs = flow_cutter::select_random_source_target_pairs(node_count, config.cutter_count, config.random_seed);
			cutter.init(flow_cutter::expanded_graph::expand_source_target_pair_list(pairs), config.random_seed);

			std::vector<flow_cutter::expanded_graph::Separator>separator_set;

			for(;;){
				auto s = flow_cutter::expanded_graph::extract_original_separator(tail, head, cutter);
				if(s.small_side_size != 0)
					separator_set.push_back(std::move(s));
				if(!cutter.advance())
					break;
			}

			std::vector<std::vector<int>>ret;
			if(!separator_set.empty()){
				for(int i=separator_set.size()-1; i>0; --i){
					separator_set[i].small_side_size -= separator_set[i-1].small_side_size;		
				}
				separator_set.front().small_side_size = 0;

				std::sort(
					separator_set.begin(), separator_set.end(), 
					[](const flow_cutter::expanded_graph::Separator&l, const flow_cutter::expanded_graph::Separator&r){
						return l.small_side_size * r.sep.size() > r.small_side_size * l.sep.size();
					}
				);

			
				for(int i=0; i<std::min(config.branch_factor, (int)separator_set.size()); ++i)
					ret.push_back(std::move(separator_set[i].sep));


			}

			/*for(auto&x:separator_set)
				std::sort(x.sep.begin(), x.sep.end());

			while(ret.size() < config.branch_factor && !separator_set.empty()){


				int best = -1;
				int best_size = std::numeric_limits<int>::max();
				int best_small_side_size = 0;

				for(int i=0; i<separator_set.size(); ++i){
					if((long long)separator_set[i].small_side_size * (long long)best_size >= (long long)best_small_side_size * (long long)separator_set[i].sep.size()){
						best = i;
						best_size = separator_set[i].sep.size();
						best_small_side_size = separator_set[i].small_side_size;
					}
				}

				assert(best != -1);

				auto c = std::move(separator_set[best].sep);

				if(best != separator_set.size()-1)
					separator_set[best] = std::move(separator_set.back());
				separator_set.pop_back();

				auto compute_set_intersection_size = [&](const std::vector<int>&l, const std::vector<int>&r){
					int n = 0;
					auto l_iter = l.begin(), l_end = l.end();
					auto r_iter = r.begin(), r_end = r.end();
					while(l_iter != l_end && r_iter != r_end){
						if(*l_iter == *r_iter){
							++l_iter;
							++r_iter;
							++n;
						}else if(*l_iter < *r_iter){
							++l_iter;
						}else if(*l_iter > *r_iter){
							++r_iter;
						}
					}
					return n;
				};

				auto is_too_similar = [&](const flow_cutter::expanded_graph::Separator&s){
					return 3*compute_set_intersection_size(c, s.sep) > (int)c.size()/4;
				};
				separator_set.erase(std::remove_if(separator_set.begin(), separator_set.end(), is_too_similar), separator_set.end());

				ret.push_back(std::move(c));
			}*/
			
			return ret;

		}
	private:
		flow_cutter::Config config;
	};


	template<class ComputeSeparatorSet, class ComputeOrder>
	ArrayIDIDFunc compute_order_by_decompose_along_all_separators(
		ArrayIDIDFunc tail, ArrayIDIDFunc head, ArrayIDIDFunc input_node_id,
		const ComputeSeparatorSet&compute_separator_set,
		const ComputeOrder&compute_order
	){
		const int node_count = tail.image_count();

		ArrayIDIDFunc order;
		int smallest_width = std::numeric_limits<int>::max();

		auto separator_set = compute_separator_set(tail, head);
		if(separator_set.empty())
			return chain(compute_greedy_min_shortcut_order(tail, head), input_node_id);
	
static int recurse_level = 0;
++recurse_level;
/*
{
int arc_count = tail.preimage_count(); 
ofstream out("foo.dot");
out << "graph G{\n";
for(int i=0; i<arc_count; ++i){
	if(tail(i) < head(i))
		out << tail(i) << " -- "<<head(i) << ";\n"; 
}
out << "}" << endl;
}*/
/*
		cout << "Separators:" << endl;
		for(auto x: separator_set){
			for(auto y:x)
				cout << y << ":";
			cout << "\b\n";
		}

{
cout << "min shortcut:" << endl;
		
auto foo = compute_greedy_min_shortcut_order(tail, head);
for(auto y:foo)
	cout << y << ":";
cout << "\b\n";
}*/
		//cin.get();

//cout << "Foo" << endl;
		for(auto sep:separator_set){

			auto sub_input_node_id = id_id_func(node_count, node_count, [](int x){return x;});
			ArrayIDIDFunc order_of_iteration = compute_order_by_upholding_separator(tail, head, sub_input_node_id, sep, compute_order);
			assert(is_permutation(order_of_iteration));
			int width_of_iteration = compute_tree_width(tail, head, order_of_iteration);

//if(recurse_level == 1){
//	for(auto x:sep)
//		cout << x << ":";
//	cout << "\b\n";
//	cout << "width "<< width_of_iteration << endl;
//}
			if(width_of_iteration < smallest_width){
				order = chain(order_of_iteration, input_node_id);
				smallest_width = width_of_iteration;
			}
		}
		
--recurse_level;

		return order;
		/*auto order = compute_greedy_min_shortcut_order(tail, head);
		assert(is_permutation(order));
		order.set_image_count(input_node_id.image_count());
		for(auto&x:order)
			x = input_node_id(x);
		assert(cch_order::is_valid_partial_order(order));
		return order;*/
	}


	template<class ComputeSeparatorSet>
	ArrayIDIDFunc compute_low_tree_width_order(
		ArrayIDIDFunc tail, ArrayIDIDFunc head,
		ArrayIDIDFunc input_node_id,
		const ComputeSeparatorSet&compute_separator_set
	){
		auto recurse = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_low_tree_width_order(
				std::move(a_tail), std::move(a_head), 
				std::move(a_input_node_id),
				compute_separator_set
			);
		};

		auto orderer4 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_decompose_along_all_separators(
				std::move(a_tail), std::move(a_head), std::move(a_input_node_id),
				compute_separator_set, recurse
			);
		};

		auto orderer3 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_applying_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer4);
		};

		auto orderer2 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id){
			return compute_order_by_decomposing_along_articulation_points(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer3);
		};

		auto orderer1 = [&](ArrayIDIDFunc a_tail, ArrayIDIDFunc a_head, ArrayIDIDFunc a_input_node_id, ArrayIDFunc<int>a_weight){
			(void)a_weight;
			return compute_order_by_applying_reduction_rules(std::move(a_tail), std::move(a_head), std::move(a_input_node_id), orderer2);
		};

		auto weight = id_func(tail.preimage_count(), [](int){return 0;});
		auto order = cch_order::reorder_nodes_in_preorder_and_compute_unconnected_graph_order_if_component_is_non_trivial(
			std::move(tail), std::move(head), std::move(input_node_id), std::move(weight),
			orderer1, [](int){return false;}
		);

		return order; // NVRO
	}

	template<class ComputeSeparatorSet>
	ArrayIDIDFunc compute_low_tree_width_order(
		ArrayIDIDFunc tail, ArrayIDIDFunc head,
		const ComputeSeparatorSet&compute_separator_set
	){
		ArrayIDFunc<int>weight = id_func(tail.preimage_count(), [](int){return 0;});
		cch_order::make_graph_simple(tail, head, weight);
		return compute_low_tree_width_order(std::move(tail), std::move(head), identity_permutation(tail.image_count()), compute_separator_set);
	}

}

#endif


