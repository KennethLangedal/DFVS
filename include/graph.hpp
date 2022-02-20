#pragma once
#include <iostream>
#include <vector>

#include "bitvector.hpp"

constexpr size_t N = 64;

class graph {
   private:
    std::vector<bitvector<N>> out_edges, in_edges;
    bitvector<N> active;

   public:
    void deactive_edge(size_t u, size_t v);
    void deactive_vertex(size_t u);
    void deactive_single_out(size_t u);
    void deactive_single_in(size_t u);

    size_t degree(size_t u) const;
    size_t in_degree(size_t u) const;
    bool self_loop(size_t u) const;
    
    const bitvector<N>& active_vertices() const;
    const bitvector<N>& out(size_t u) const;
    const bitvector<N>& in(size_t u) const;

    void print_edgelist(std::ostream& os) const;
    friend std::istream& operator>>(std::istream& is, graph& g);
};

std::istream& operator>>(std::istream& is, graph& g);