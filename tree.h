#ifndef TREE_H
#define TREE_H

#include "array_id_func.h"
#include "tiny_id_func.h"
#include <cstdlib>
#include <utility>

template<class Neighbors>
bool is_tree(const Neighbors&neighbors){
	int node_count = neighbors.preimage_count();
	if(node_count <= 1)
		return true;

	int reachable_node_count = 0;

	ArrayIDFunc<int>parent(node_count);
	parent.fill(-1);
	ArrayIDFunc<int>stack(node_count);
	int stack_end = 0;
	stack[stack_end++] = 0;
	while(stack_end != 0){
		++reachable_node_count;
		auto x = stack[--stack_end];
		for(auto y:neighbors(x)){
			if(parent(x) != y){
				if(parent(y) == -1){
					parent[y] = x;
					stack[stack_end++] = y;
				}else{
					return false;
				}
			}
		}
	}
	return reachable_node_count == node_count;
}

std::pair<ArrayIDIDFunc, ArrayIDIDFunc> generate_random_tree(int node_count){
	const int arc_count = 2*(node_count-1);

	std::pair<ArrayIDIDFunc, ArrayIDIDFunc>p;
	auto&tail = p.first;
	auto&head = p.second;

	tail = ArrayIDIDFunc(arc_count, node_count);
	head = ArrayIDIDFunc(arc_count, node_count);

	int arc_id = 0;

	for(int i=1; i<node_count; ++i){
		int p = std::rand()%i;
		tail[arc_id] = i;
		head[arc_id] = p;
		++arc_id;
		tail[arc_id] = p;
		head[arc_id] = i;
		++arc_id;
	}

	return p; // NVRO
}

#endif

