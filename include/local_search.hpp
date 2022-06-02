#pragma once
#include "sparse_graph.hpp"
#include <random>

class local_search {
private:
    std::vector<uint32_t> _next, _prev, _vertices;
    std::vector<std::vector<uint32_t>> _out_edges_in_config, _in_edges_in_config, _out_edges_not_config, _in_edges_not_config;
    uint32_t _first, _last, _N, _current_cost, _best_cost;

    std::vector<uint64_t> _config_order;

    std::vector<bool> _config, _best_config;

    std::mt19937 _reng;
    std::uniform_int_distribution<uint32_t> _dist_int;
    std::uniform_real_distribution<double> _dist_real;

    double _T;

    void _insert(uint32_t u, uint32_t i); // insert before i
    void _remove(uint32_t u);
    void _assign_label(uint32_t u);
    void _relable();

    void _search_step();

    std::tuple<uint32_t, float> _move_score_first_out(uint32_t u) const;
    std::tuple<uint32_t, float> _move_score_last_in(uint32_t u) const;

    void _apply_move(uint32_t u, uint32_t i, bool pos);

    void _greedy_strictly_improving_move(const sparse_graph &g, uint32_t u);

    void _one_one_candidates(uint32_t u, std::vector<std::vector<uint32_t>> &cand_list);

    bool _two_one_swap_test(const sparse_graph &g, uint32_t u, uint16_t v, uint32_t w);

public:
    local_search(const sparse_graph &g, double T = 0.4, size_t seed = 0);

    void search(const sparse_graph &g, size_t iterations);

    bool check_every_two_one_swap(const sparse_graph &g);

    void set_temperature(double T);

    void greedy_one_zero_swaps();

    void greedy_one_zero_swaps_dfs(const sparse_graph &g);

    void shuffle_solution(const sparse_graph &g);

    bitvector get_best() const;

    float get_average_move_cost(const sparse_graph &g) const;

    uint32_t get_best_cost() const;

    uint32_t get_current_cost() const;
};