#!/bin/sh

mpicxx -I. -std=c++0x -DUSE_KAHIP -O3 -DNEBUG console.cpp .fancy_input.o .greedy_order.o .permutation.o .list_graph.o .file_utility.o .geo_pos.o -lpthread -lreadline -L. -lkahip -fopenmp -lm -o console_with_kahip

