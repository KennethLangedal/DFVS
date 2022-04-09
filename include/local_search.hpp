#pragma once

#include "reductions.hpp"

void reducing_peeling(graph &g, bitvector &dfvs, graph_search &gs, std::vector<bitvector> &SCC);

void remove_redundant(const graph &g, bitvector &dfvs);

void two_one_swap(const graph &g, bitvector &dfvs);

void analysis(graph &g);