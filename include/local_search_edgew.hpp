#pragma once
#include "sparse_graph.hpp"
#include <random>

class local_search_edgew {
private:
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> _out_edges, _in_edges;
    std::vector<uint32_t> _edge_weights, _vertices;

    std::vector<uint32_t> _next, _prev;
    uint32_t _first, _last, _N;

    std::vector<uint64_t> _config_order;

    bitvector _config, _not_config, _best_dfvs;

    std::mt19937 _reng;

    void _insert(uint32_t u, uint32_t i);
    void _remove(uint32_t u);
    void _assign_label(uint32_t u);
    void _relable();

    std::tuple<uint32_t, int32_t> _move_score_first_out(uint32_t u);
    std::tuple<uint32_t, int32_t> _move_score_last_in(uint32_t u);

    void _apply_move(uint32_t u, uint32_t i, bool pos);

public:
    local_search_edgew(const sparse_graph &g);

    void search(size_t iterations);

    void shuffle_solution();

    const bitvector &get_best() const;

    const bitvector &get_current() const;
};