#pragma once
#include "sparse_graph.hpp"
#include <random>

class local_search {
private:
    std::vector<uint32_t> _next, _prev, _vertices;
    std::vector<std::vector<uint32_t>> _out_edges_in_config, _in_edges_in_config, _out_edges_not_config, _in_edges_not_config;
    uint32_t _first, _last, _N, _current_cost, _best_cost;

    std::vector<uint64_t> _config_order;

    bitvector _not_config, _best_dfvs;

    std::mt19937 _reng;
    std::uniform_int_distribution<uint32_t> _dist_int;
    std::uniform_real_distribution<double> _dist_real;

    double _T;

    void _insert(uint32_t u, uint32_t i); // insert before i
    void _remove(uint32_t u);
    void _assign_label(uint32_t u);
    void _relable();

    void _search_step(const sparse_graph &g);

    std::tuple<uint32_t, float> _move_score_first_out(const sparse_graph &g, uint32_t u) const;
    std::tuple<uint32_t, float> _move_score_last_in(const sparse_graph &g, uint32_t u) const;
    void _apply_move(const sparse_graph &g, uint32_t u, uint32_t i, bool pos);

public:
    local_search(const sparse_graph &g, double T = 0.4, size_t seed = 0);

    void search(const sparse_graph &g, size_t iterations);

    void set_temperature(double T);

    void set_solution(const sparse_graph &g, const bitvector &fvs);

    void greedy_one_zero_swaps(const sparse_graph &g);

    void shuffle_solution(const sparse_graph &g);

    const bitvector &get_best() const;

    const bitvector &get_current() const;

    const uint32_t &get_best_cost() const;

    const uint32_t &get_current_cost() const;
};