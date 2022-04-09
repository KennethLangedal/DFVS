#pragma once
#include "reductions.hpp"

bitvector solve_req(graph &g, graph_search &gs, size_t cost, size_t ub, size_t d, size_t &lb_counter);

bitvector solve(graph &g);