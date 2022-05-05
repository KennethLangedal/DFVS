#pragma once
#include "bitvector.hpp"
#include "sparse_graph.hpp"
#include <random>

class sparse_local_search {
private:
    std::vector<uint32_t> next, prev, vertices;
    uint32_t first, last, N;

    std::vector<uint64_t> config_order;

    bitvector config, best_dfvs;

    std::mt19937 reng;
    std::uniform_int_distribution<uint32_t> dist_int;
    std::uniform_real_distribution<double> dist_real;

    double T;

    void insert(uint32_t u, uint32_t i);
    void remove(uint32_t u);
    void assign_label(uint32_t u);
    void relable();

    void search_step(const sparse_graph &g);
    std::pair<uint32_t, float> move_score(const sparse_graph &g, uint32_t u, bool pos) const;
    void apply_move(const sparse_graph &g, uint32_t u, uint32_t i, bool pos);

public:
    sparse_local_search(const sparse_graph &g, double T = 0.4, size_t seed = 0);

    void search(const sparse_graph &g, size_t iterations);

    void set_temperature(double T);

    void set_solution(const sparse_graph &g, const bitvector &fvs);

    void greedy_one_zero_swaps(const sparse_graph &g);

    void greedy_two_one_swaps(const sparse_graph &g);

    void shuffle_solution(const sparse_graph &g);

    const bitvector &get_best() const;

    size_t get_current() const;
};