#ifndef GREEDY_ORDER_H
#define GREEDY_ORDER_H

#include "array_id_func.h"

ArrayIDIDFunc compute_greedy_min_degree_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);
ArrayIDIDFunc compute_greedy_min_shortcut_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);
ArrayIDIDFunc compute_greedy_independent_set_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, bool degree_guided);
ArrayIDIDFunc compute_greedy_min_shortcut_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, int hash_factor, int hash_modulo);
ArrayIDIDFunc compute_minimum_elimination_tree_height_order_from_chordal_graph(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);
ArrayIDIDFunc compute_greedy_min_shortcut_and_level_order(const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head);

#endif
