#include "greedy_order.h"
#include "array_id_func.h"
#include "tiny_id_func.h"
#include "permutation.h"
#include "heap.h"
#include "min_max.h"
#include <vector>

ArrayIDFunc<std::vector<int>> build_dyn_array(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head){
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	ArrayIDFunc<std::vector<int>> neighbors(node_count);

	for(int i=0; i<arc_count; ++i)
		neighbors[tail(i)].push_back(head(i));

	for(int i=0; i<node_count; ++i)
		std::sort(neighbors[i].begin(), neighbors[i].end());

	return neighbors; // NVRO
}

template<class T>
struct NullAssign{
public:
	NullAssign&operator=(const T&){
		return *this;
	}
};

template<class T>
struct CountOutputIterator{
	typedef T value_type;
	typedef int difference_type;
	typedef T*pointer;
	typedef T&reference;
	typedef std::output_iterator_tag iterator_category;

	NullAssign<T> operator*()const{
		return {};
	}

	CountOutputIterator(int&n):
		n(&n){};

	CountOutputIterator&operator++(){
		++*n;
		return *this;
	}

	CountOutputIterator&operator++(int){
		++*n;
		return *this;
	}

	int*n;
};


template<class Iter1, class Iter2, class Iter3, class T>
Iter3 set_union_and_remove_element(
	Iter1 a, Iter1 a_end,
	Iter2 b, Iter2 b_end,
	Iter3 out,
	const T&remove_element1, const T&remove_element2
){
	while(a != a_end && b != b_end){
		if(*a < *b){
			if(*a != remove_element1 && *a != remove_element2)
				*out++ = *a;
			++a;
		}else if(*a > *b){
			if(*b != remove_element1 && *b != remove_element2)
				*out++ = *b;
			++b;
		}else if(*a == *b){
			if(*a != remove_element1 && *a != remove_element2)
				*out++ = *a;
			++b;
			++a;
		}
	}

	while(a != a_end){
		if(*a != remove_element1 && *a != remove_element2)
			*out++ = *a;
		++a;
	}

	while(b != b_end){
		if(*b != remove_element1 && *b != remove_element2)
			*out++ = *b;
		++b;
	}

	return out;
}

std::vector<int> contract_node(ArrayIDFunc<std::vector<int>>&graph, int node){
	std::vector<int>tmp;
	for(int x:graph(node)){
		tmp.clear();
		set_union_and_remove_element(
			graph(node).begin(), graph(node).end(),
			graph(x).begin(), graph(x).end(),
			std::back_inserter(tmp),
			node, 
			x
		);
		graph[x].swap(tmp);
	}

	return std::move(graph[node]);
}

int compute_number_of_shortcuts_added_if_contracted(const ArrayIDFunc<std::vector<int>>&graph, int node){
	int added = 0;
	for(int x:graph(node)){
		std::set_difference(
			graph(node).begin(), graph(node).end(),
			graph(x).begin(), graph(x).end(),
			CountOutputIterator<int>(added)
		);
		--added;
	}

	added /= 2;

	return added;
}

ArrayIDIDFunc compute_greedy_min_degree_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head){
	const int node_count = tail.image_count();
	
	auto g = build_dyn_array(tail, head);

	min_id_heap<int> q(node_count);

	for(int x=0; x<node_count; ++x)
		q.push(x, g(x).size());

	ArrayIDIDFunc order(node_count, node_count);
	int next_pos = 0;

	while(!q.empty()){
		auto x = q.pop();

		order[next_pos++] = x;

		for(auto y:contract_node(g, x)){
			q.push_or_set_key(y, g(y).size());
		}
	}

	return order; // NVRO
}

ArrayIDIDFunc compute_greedy_min_shortcut_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head){
	const int node_count = tail.image_count();

	auto g = build_dyn_array(tail, head);

	min_id_heap<int> q(node_count);

	for(int x=0; x<node_count; ++x)
		q.push(x, 100*compute_number_of_shortcuts_added_if_contracted(g,x) +  g(x).size());

	ArrayIDIDFunc order(node_count, node_count);
	int next_pos = 0;

	while(!q.empty()){
		auto x = q.pop();

		order[next_pos++] = x;

		for(auto y:contract_node(g, x)){
			q.push_or_set_key(y, 100*compute_number_of_shortcuts_added_if_contracted(g,y) + g(y).size());
		}
	}

	return order; // NVRO
}

ArrayIDIDFunc compute_greedy_min_shortcut_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, int hash_factor, int hash_modulo){
	const int node_count = tail.image_count();

	auto g = build_dyn_array(tail, head);

	min_id_heap<int> q(node_count);

	auto compute_key = [&](int x){
		int key = 100*compute_number_of_shortcuts_added_if_contracted(g,x) + g(x).size();
		if(key < 10)
			return key;
		else
			return key + (hash_factor * x) % (2*key/3) - key/3;
	};

	for(int x=0; x<node_count; ++x)
		q.push(x, compute_key(x));

	ArrayIDIDFunc order(node_count, node_count);
	int next_pos = 0;

	while(!q.empty()){
		auto x = q.pop();

		order[next_pos++] = x;

		for(auto y:contract_node(g, x))
			q.push_or_set_key(y, compute_key(y));
	}

	return order; // NVRO
}

ArrayIDIDFunc compute_greedy_min_shortcut_and_level_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head){
	const int node_count = tail.image_count();

	auto g = build_dyn_array(tail, head);

	min_id_heap<int> q(node_count);

	ArrayIDFunc<int> level(node_count);
	level.fill(0);

	auto compute_key = [&](int x){
		return compute_number_of_shortcuts_added_if_contracted(g,x) + level(x);
	};

	for(int x=0; x<node_count; ++x)
		q.push(x, compute_key(x));

	ArrayIDIDFunc order(node_count, node_count);
	int next_pos = 0;

	while(!q.empty()){
		auto x = q.pop();

		order[next_pos++] = x;

		for(auto y:contract_node(g, x)){
			max_to(level[y], level[x]+1);
			q.push_or_set_key(y, compute_key(y));	
		}
	}

	return order; // NVRO
}


ArrayIDIDFunc compute_greedy_independent_set_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, bool degree_guided){
	const int node_count = tail.image_count();
	
	auto g = build_dyn_array(tail, head);

	ArrayIDIDFunc order(node_count, node_count);
	int next_order_pos = 0;

	BitIDFunc in_independent_set(node_count);
	in_independent_set.fill(false);

	std::vector<int>not_in_order(node_count);
	for(int i=0; i<node_count; ++i)
		not_in_order[i] = i;

	while(next_order_pos != node_count){
		if(degree_guided)
			std::sort(not_in_order.begin(), not_in_order.end(), [&](int l, int r){return g(l).size() < g(r).size();});

		for(auto x:not_in_order){
			bool put_x_in_set = true;
			for(auto y:g(x)){
				if(in_independent_set(y)){
					put_x_in_set = false;
					break;
				}
			}
			if(put_x_in_set)
				in_independent_set.set(x, true);
		}

		auto i = std::partition(not_in_order.begin(), not_in_order.end(), [&](int x){return !in_independent_set(x);});

		for(auto j=i; j!=not_in_order.end(); ++j){
			order[next_order_pos++] = *j;
			in_independent_set.set(*j, false);

			contract_node(g, *j);
		}

		not_in_order.erase(i, not_in_order.end());

	}
	
	return order; // NVRO
}


bool is_node_simplicial(const ArrayIDFunc<std::vector<int>>&graph, int node){
	int node_degree = std::end(graph(node)) - std::begin(graph(node));

	#ifndef NDEBUG
	bool correct_answer = compute_number_of_shortcuts_added_if_contracted(graph, node)==0;
	#endif

	if(node_degree <= 1){
		assert(correct_answer == true);
		return true;
	}

	for(auto neighbor:graph(node)){
		int neighbor_degree = std::end(graph(neighbor)) - std::begin(graph(neighbor));
		if(node_degree > neighbor_degree){
			assert(correct_answer == false);
			return false;
		}
	}

	for(int neighbor:graph(node)){

		auto neighbor_iter = graph(neighbor).begin(), neighbor_end = graph(neighbor).end();
		auto node_iter = graph(node).begin(), node_end = graph(node).end();

		for(;;){
			if(node_iter == node_end)
				break;	
			
			if(*node_iter == neighbor){
				++node_iter;
				continue;
			}
			if(neighbor_iter != neighbor_end && *neighbor_iter == node){
				++neighbor_iter;
				continue;
			}
			if(neighbor_iter == neighbor_end){
				assert(correct_answer == false);
				return false;
			}
			if(*node_iter == *neighbor_iter){
				++node_iter;
				++neighbor_iter;
				continue;
			}
			if(*neighbor_iter < *node_iter){
				++neighbor_iter;
				continue;
			}
			assert(correct_answer == false);
			return false;
		}

	}

	assert(correct_answer == true);
	return true;
}

std::vector<int> contract_simplicial_node(ArrayIDFunc<std::vector<int>>&graph, int node){
	for(int x:graph(node))
		graph[x].erase(std::lower_bound(std::begin(graph[x]), std::end(graph[x]), node));
	return std::move(graph[node]);
}


ArrayIDIDFunc compute_minimum_elimination_tree_height_order_from_chordal_graph(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head){
	int node_count = tail.image_count();
	ArrayIDIDFunc order(node_count, node_count);
	int next_order_pos = 0;

	auto g = build_dyn_array(tail, head);

	ArrayIDFunc<int> node_height(node_count);
	node_height.fill(0);

	min_id_heap<int>simplicial_node_queue(node_count);
	for(int x=0; x<node_count; ++x)
		if(is_node_simplicial(g, x)){
			simplicial_node_queue.push(x, 0);
		}

	while(!simplicial_node_queue.empty()){
		int x = simplicial_node_queue.pop();
		order[next_order_pos++] = x; 
		for(auto y:contract_node(g, x)){
			if(node_height[y] < node_height[x]+1)
				node_height[y] = node_height[x]+1;
			if(simplicial_node_queue.contains(y) || is_node_simplicial(g, y))
				simplicial_node_queue.push_or_increase_key(y, node_height[y]);
		}
		
	}
	
	if(next_order_pos < node_count)
		throw std::runtime_error("Graph was not chordal");

	return order; // NVRO

}

